// ============================================================================
//  dsr_encoder.cpp  -  DSR (Digital Satellite Radio) encoder
// ============================================================================
//
//  Up to 16 stereo programmes  ->  4-QPSK baseband (I/Q) or real IF signal.
//
//  Signal flow (one pass per 2 ms audio block = 64 main frames):
//
//   audio source ──ffmpeg──> 32 kHz s16 stereo ──> 64-sample blocks
//        │  scale factor per block and channel (look-ahead of 2 blocks = 4 ms)
//        ▼
//   16/14-bit coding ─> 77-bit blocks with BCH(63,44) ─> interleave ─> frames A/B
//   (+ Barker sync, S bit = special service SA, ZI bit = scale factors + PI)
//        ▼
//   scrambler ─> differential encoder ─> QPSK mapping ─> RRC pulse shaping
//        ▼
//   cs16 | cf32 (I/Q)  or  real16 (IF = fs/4)  ──> file / stdout / TCP port
//
//  Build:  make          Usage: ./dsr_encoder --help
// ============================================================================
#include <getopt.h>
#include <signal.h>

#include <atomic>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <memory>
#include <mutex>
#include <string>
#include <thread>

#include "dsr.hpp"
#include "dsr_dsp.hpp"
#include "dsr_io.hpp"

using namespace dsr;

static std::atomic<bool> g_stop{false};
static void on_signal(int) { g_stop = true; }

// ----------------------------------------------------------------------------
// AudioSource: delivers 32 kHz / 16 bit / stereo PCM from a background thread.
// Source specification:
//   tone:FL[,FR]   synthetic sine test tone (Hz)
//   raw:PATH       raw s16le, 32 kHz, stereo file
//   anything else  passed to ffmpeg (files: wav/mp3/flac..., URLs: internet radio)
// A ring buffer decouples the (bursty) source from the strict 2 ms encoder clock.
// ----------------------------------------------------------------------------
class AudioSource {
 public:
  AudioSource(std::string spec, std::string ffmpeg, bool loop = false)
      : spec_(std::move(spec)), ffmpeg_(std::move(ffmpeg)), loop_(loop) {
    ring_.assign(kCap * 2, 0);
  }
  ~AudioSource() { stop(); }

  void start() { th_ = std::thread(&AudioSource::worker, this); }

  void stop() {
    { std::lock_guard<std::mutex> l(m_); stop_ = true; }
    cv_.notify_all();
    if (pid_ > 0) kill(pid_, SIGTERM);
    if (th_.joinable()) th_.join();
    if (pid_ > 0) { waitpid(pid_, nullptr, 0); pid_ = -1; }
  }

  size_t available() { std::lock_guard<std::mutex> l(m_); return count_; }
  bool finished() { std::lock_guard<std::mutex> l(m_); return eof_ && count_ == 0; }

  // Copy up to `frames` stereo frames (interleaved L,R) into dst.  If `block`
  // is set, wait until enough data is there or the source has ended.
  size_t pull(int16_t* dst, size_t frames, bool block) {
    std::unique_lock<std::mutex> l(m_);
    if (block) cv_.wait(l, [&] { return count_ >= frames || eof_ || stop_; });
    size_t n = std::min(frames, count_);
    for (size_t i = 0; i < n; ++i) {
      size_t r = (rd_ + i) % kCap;
      dst[2 * i] = ring_[2 * r];
      dst[2 * i + 1] = ring_[2 * r + 1];
    }
    rd_ = (rd_ + n) % kCap;
    count_ -= n;
    l.unlock();
    cv_.notify_all();
    return n;
  }

 private:
  static constexpr size_t kCap = 3 * kAudioRate;   // 3 s of audio

  // Producer side: blocks while the ring is full.
  void push(const int16_t* d, size_t frames) {
    size_t done = 0;
    while (done < frames) {
      std::unique_lock<std::mutex> l(m_);
      cv_.wait(l, [&] { return count_ < kCap || stop_; });
      if (stop_) return;
      size_t n = std::min(frames - done, kCap - count_);
      for (size_t i = 0; i < n; ++i) {
        size_t w = (rd_ + count_ + i) % kCap;
        ring_[2 * w] = d[2 * (done + i)];
        ring_[2 * w + 1] = d[2 * (done + i) + 1];
      }
      count_ += n;
      done += n;
      l.unlock();
      cv_.notify_all();
    }
  }

  void worker() {
    if (spec_.rfind("tone:", 0) == 0) {
      double fl = 440, fr = 440;
      std::sscanf(spec_.c_str() + 5, "%lf,%lf", &fl, &fr);
      if (spec_.find(',') == std::string::npos) fr = fl;
      double pl = 0, pr = 0;
      std::vector<int16_t> buf(2 * 256);
      while (!stop_) {
        for (int i = 0; i < 256; ++i) {
          buf[2 * i] = static_cast<int16_t>(9000 * std::sin(pl));
          buf[2 * i + 1] = static_cast<int16_t>(9000 * std::sin(pr));
          pl += 2 * M_PI * fl / kAudioRate;
          pr += 2 * M_PI * fr / kAudioRate;
        }
        push(buf.data(), 256);
      }
      return;
    }
    if (spec_.rfind("raw:", 0) == 0) {
      fd_ = ::open(spec_.c_str() + 4, O_RDONLY);
    } else {
      std::vector<std::string> a = {ffmpeg_, "-nostdin", "-loglevel", "error"};
      if (spec_.rfind("http", 0) == 0)       // keep live internet streams alive
        for (const char* s : {"-reconnect", "1", "-reconnect_streamed", "1", "-reconnect_delay_max", "5"}) a.push_back(s);
      if (loop_ && spec_.rfind("http", 0) != 0) { a.push_back("-stream_loop"); a.push_back("-1"); }
      a.push_back("-i");
      a.push_back(spec_);
      for (const char* s : {"-vn", "-f", "s16le", "-acodec", "pcm_s16le", "-ac", "2", "-ar", "32000", "-"}) a.push_back(s);
      fd_ = io::spawn_reader(a, &pid_);
    }
    if (fd_ < 0) std::fprintf(stderr, "cannot open source '%s'\n", spec_.c_str());
    std::vector<uint8_t> raw(16384);
    size_t have = 0;
    while (fd_ >= 0 && !stop_) {
      ssize_t r = io::read_some(fd_, raw.data() + have, raw.size() - have);
      if (r <= 0) break;
      have += static_cast<size_t>(r);
      size_t frames = have / 4;                        // 4 bytes per stereo frame
      push(reinterpret_cast<const int16_t*>(raw.data()), frames);
      size_t used = frames * 4;
      std::memmove(raw.data(), raw.data() + used, have - used);
      have -= used;
    }
    if (fd_ >= 0) close(fd_);
    { std::lock_guard<std::mutex> l(m_); eof_ = true; }
    cv_.notify_all();
  }

  std::string spec_, ffmpeg_;
  bool loop_;
  std::vector<int16_t> ring_;
  size_t rd_ = 0, count_ = 0;
  std::mutex m_;
  std::condition_variable cv_;
  bool eof_ = false, stop_ = false;
  pid_t pid_ = -1;
  int fd_ = -1;
  std::thread th_;
};

// ----------------------------------------------------------------------------
// Output sample format
// ----------------------------------------------------------------------------
enum class Format { CS16, CF32, REAL16, CS8 };

// ----------------------------------------------------------------------------
// The encoder proper
// ----------------------------------------------------------------------------
struct Options {
  std::string src[kStations];
  StationInfo st[kStations];
  std::string out = "dsr.cs16";
  Format fmt = Format::CS16;
  int sps = 2;
  bool realtime = false, realtime_set = false;
  int prefill_ms = 300;
  double duration = 0;       // seconds, 0 = until sources end
  double amp = 0.18;         // rms of one I/Q component relative to full scale
  std::string ffmpeg = "ffmpeg";
  bool verbose = false;
  bool loop = false;         // repeat file sources forever
  bool conj = false;         // mirror the spectrum (negate Q)
};

// One 2 ms block of PCM for all stations plus the derived scale factors.
struct BlockData {
  int16_t pcm[kStations][2][kBlockSamples];
  uint8_t sk[kStations][2];
};
// Pre-computed per-frame content of one block (so the frame loop is trivial).
struct BlockPlan {
  u16 cw[kStations][2][kBlockSamples];
  u8  zi[kStations][kBlockSamples];
};

class Encoder {
 public:
  explicit Encoder(const Options& o) : o_(o), shaper_(o.sps) {
    for (int s = 0; s < kStations; ++s)
      if (o.st[s].active) src_[s] = std::make_unique<AudioSource>(o.src[s], o.ffmpeg, o.loop);
    // Normalise so that one I/Q component has rms = amp * full scale.
    scale_ = o.amp / std::sqrt(shaper_.unit_power());
  }

  int run() {
    // ---- open the sink ----
    if (o_.out == "-") fd_ = 1;
    else if (o_.out.rfind("tcp://", 0) == 0) {
      // Server that stays open: decoders may come and go while the encoder runs.
      const char* p = std::strrchr(o_.out.c_str(), ':');
      int port = p ? std::atoi(p + 1) : 5000;
      // queue limit ~1.5 s of stream data per client
      size_t bps = (o_.fmt == Format::CF32 ? 8 : o_.fmt == Format::CS16 ? 4 : 2) * size_t(o_.sps * kSymbolRate);   // cs8 and real16: 2 bytes/sample
      if (!bc_.start(port, bps * 3 / 2)) return 1;
      std::fprintf(stderr, "serving on 127.0.0.1:%d - decoders can connect/disconnect at any time\n", port);
      use_tcp_ = true;
    } else fd_ = ::open(o_.out.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd_ < 0 && !use_tcp_) { std::perror("output"); return 1; }

    for (auto& s : src_) if (s) s->start();
    if (o_.realtime) wait_prefill();

    // ---- prime the 3-block look-ahead queue (blocks 0,1,2) ----
    bool more = true;
    for (int i = 0; i < 3; ++i) more = read_block(queue_[i]) || more;
    std::memset(&plan_, 0, sizeof plan_);
    make_plan_silent();

    int ending = more ? -1 : 3;
    const int64_t t0 = io::now_ns();
    uint64_t blocks = 0;
    const int64_t block_ns = 2000000;      // 2 ms

    // The frame grid starts with the first 16 bits of a special-service frame
    // (Williard word); the audio block grid begins right after it (spec 3.3.2).
    // Frames 0..15 therefore belong to a silent "pre-roll" block (q = 48..63).
    emit_frames(48, 64);

    while (!g_stop) {
      make_plan(queue_[0], queue_[2]);     // audio of block m, scale factors of block m+2
      emit_frames(0, 64);
      ++blocks;

      // advance the queue and fetch block m+3
      queue_[0] = queue_[1];
      queue_[1] = queue_[2];
      bool got = read_block(queue_[2]);
      if (ending < 0 && !got) ending = 3;
      if (ending >= 0 && --ending < 0) break;
      if (o_.duration > 0 && blocks * 0.002 >= o_.duration) break;

      if (o_.realtime) io::sleep_until_ns(t0 + int64_t(blocks) * block_ns);
    }
    flush();
    if (use_tcp_) bc_.stop();
    if (fd_ > 2) close(fd_);
    std::fprintf(stderr, "done: %llu blocks (%.2f s)\n", (unsigned long long)blocks, blocks * 0.002);
    return 0;
  }

 private:
  // Wait until every active source has some audio buffered (live streams).
  void wait_prefill() {
    size_t need = size_t(o_.prefill_ms) * kAudioRate / 1000;
    int64_t limit = io::now_ns() + 15000000000ll;
    for (auto& s : src_) {
      if (!s) continue;
      while (s->available() < need && !s->finished() && io::now_ns() < limit && !g_stop) usleep(10000);
    }
  }

  // Read one block of 64 stereo samples for every station.  Returns true if any
  // station still delivered audio.  In real-time mode missing samples become silence.
  bool read_block(BlockData& b) {
    std::memset(&b, 0, sizeof b);
    bool any = false;
    int16_t tmp[2 * kBlockSamples];
    for (int s = 0; s < kStations; ++s) {
      b.sk[s][0] = b.sk[s][1] = 7;
      if (!src_[s]) continue;
      std::memset(tmp, 0, sizeof tmp);
      size_t n = src_[s]->pull(tmp, kBlockSamples, !o_.realtime);
      if (n > 0 || !src_[s]->finished()) any = true;
      for (int i = 0; i < kBlockSamples; ++i) {
        b.pcm[s][0][i] = tmp[2 * i];
        b.pcm[s][1][i] = tmp[2 * i + 1];
      }
      b.sk[s][0] = static_cast<uint8_t>(block_scale_factor(b.pcm[s][0]));
      b.sk[s][1] = static_cast<uint8_t>(block_scale_factor(b.pcm[s][1]));
    }
    return any;
  }

  // Frame content for the silent pre-roll block.
  void make_plan_silent() {
    for (int s = 0; s < kStations; ++s)
      for (int q = 0; q < kBlockSamples; ++q) {
        bool act = o_.st[s].active;
        plan_.cw[s][0][q] = plan_.cw[s][1][q] = act ? 0 : 0x3FFF;
        plan_.zi[s][q] = act ? 0 : 1;
      }
  }

  // Build the per-frame code words and ZI bits of one block.
  //   audio : block m (coded with its own scale factors)
  //   look  : block m+2 - its scale factors go into the ZI frame, i.e. they are
  //           transmitted 4 ms before the audio they belong to (spec 3.5.2).
  void make_plan(const BlockData& audio, const BlockData& look) {
    u32 pi = pi_.next_word();
    for (int s = 0; s < kStations; ++s) {
      if (!o_.st[s].active) {            // unused channel: "Dauer-Eins" (spec 3.4.2.1)
        for (int q = 0; q < kBlockSamples; ++q) {
          plan_.cw[s][0][q] = plan_.cw[s][1][q] = 0x3FFF;
          plan_.zi[s][q] = 1;
        }
        continue;
      }
      for (int c = 0; c < 2; ++c)
        for (int q = 0; q < kBlockSamples; ++q)
          plan_.cw[s][c][q] = encode_sample(audio.pcm[s][c][q], audio.sk[s][c]);
      // ZI frame: 3 x (SkL, SkR, BCH) = 42 bit, then 22 PI bits (Bild 5)
      u8 w[14];
      zi_encode(look.sk[s][0], look.sk[s][1], w);
      for (int r = 0; r < 3; ++r)
        for (int i = 0; i < 14; ++i) plan_.zi[s][r * 14 + i] = w[i];
      for (int i = 0; i < 22; ++i) plan_.zi[s][42 + i] = (pi >> (21 - i)) & 1u;
    }
  }

  // Emit main frames q0 .. q1-1 of the current plan.
  void emit_frames(int q0, int q1) {
    for (int q = q0; q < q1; ++q) {
      const uint64_t n = frame_no_++;
      // ---- special-service bit of frame A (SA frame = 64 frames) ----
      int sa_bit = static_cast<int>(n % 64);
      uint64_t sa_no = n / 64;
      if (sa_bit == 0) {
        make_sa_bytes(o_.st, sa_no, sa_bytes_);
        sa_sync_ = (sa_no % 8 == 0) ? kSync1 : kSync2;
      }
      u8 sbit;
      if (sa_bit < 16) sbit = (sa_sync_ >> (15 - sa_bit)) & 1u;
      else sbit = (sa_bytes_[(sa_bit - 16) / 8] >> (7 - (sa_bit - 16) % 8)) & 1u;

      u8 fr[2][kFrameBits];
      for (int half = 0; half < 2; ++half) {
        StationSlot slot[8];
        for (int i = 0; i < 8; ++i) {
          int s = half * 8 + i;
          slot[i].cw[0] = plan_.cw[s][0][q];
          slot[i].cw[1] = plan_.cw[s][1][q];
          slot[i].zi = plan_.zi[s][q];
        }
        build_frame(half == 1, slot, half == 0 ? sbit : 0, fr[half]);   // S bit of B is 0
        scramble_frame(half == 1, fr[half]);
      }
      // ---- differential encoding + QPSK mapping + pulse shaping ----
      for (int i = 0; i < kFrameBits; ++i) {
        u8 A2, B2;
        diff_.step(fr[0][i], fr[1][i], A2, B2);
        // Bild 10: passband = A'' sin(wt) + B'' cos(wt)  ==  Re{(B'' - jA'') e^{jwt}}
        // bit 0 -> +1, bit 1 -> -1.
        float I = B2 ? -1.f : 1.f;
        float Q = (A2 ? 1.f : -1.f) * (o_.conj ? -1.f : 1.f);
        float oi[8], oq[8];
        shaper_.push(I, Q, oi, oq);
        for (int p = 0; p < o_.sps; ++p) put_sample(oi[p] * scale_, oq[p] * scale_);
      }
      if (out_i16_.size() * 2 + out_f32_.size() * 4 + out_i8_.size() >= (1u << 18)) flush();
    }
  }

  // Convert one baseband sample to the selected output format.
  void put_sample(float i, float q) {
    switch (o_.fmt) {
      case Format::CS16:
        out_i16_.push_back(clip(i)); out_i16_.push_back(clip(q)); break;
      case Format::CS8:                  // signed 8-bit I/Q (HackRF / SoapySDR "CS8"), full scale = +-127
        out_i8_.push_back(clip8(i)); out_i8_.push_back(clip8(q)); break;
      case Format::CF32:
        out_f32_.push_back(i); out_f32_.push_back(q); break;
      case Format::REAL16: {             // IF = fs/4: I*cos(pi n/2) - Q*sin(pi n/2)
        float v;
        switch (sample_no_ & 3) { case 0: v = i; break; case 1: v = -q; break; case 2: v = -i; break; default: v = q; }
        out_i16_.push_back(clip(v));
        break;
      }
    }
    ++sample_no_;
  }
  // Quantise with rounding (truncation would add 6 dB more quantisation noise,
  // which matters for 8 bit) and hard-limit at full scale.
  static int16_t clip(float v) {
    float x = v * 32767.f;
    return static_cast<int16_t>(std::lrintf(std::max(-32767.f, std::min(32767.f, x))));
  }
  static int8_t clip8(float v) {
    float x = v * 127.f;
    return static_cast<int8_t>(std::lrintf(std::max(-127.f, std::min(127.f, x))));
  }

  void flush() {
    if (!out_i16_.empty()) { put(out_i16_.data(), out_i16_.size() * 2); out_i16_.clear(); }
    if (!out_i8_.empty()) { put(out_i8_.data(), out_i8_.size()); out_i8_.clear(); }
    if (!out_f32_.empty()) { put(out_f32_.data(), out_f32_.size() * 4); out_f32_.clear(); }
  }
  void put(const void* d, size_t n) {
    if (use_tcp_) bc_.send(d, n);          // no client connected -> data is simply discarded
    else io::write_all(fd_, d, n);
  }

  Options o_;
  std::unique_ptr<AudioSource> src_[kStations];
  BlockData queue_[3];
  BlockPlan plan_;
  PiGenerator pi_;
  PulseShaper shaper_;
  DiffEncoder diff_;
  float scale_;
  int fd_ = -1;
  bool use_tcp_ = false;
  io::TcpBroadcaster bc_;
  uint64_t frame_no_ = 0, sample_no_ = 0;
  u8 sa_bytes_[6] = {0};
  u16 sa_sync_ = kSync1;
  std::vector<int16_t> out_i16_;
  std::vector<float> out_f32_;
  std::vector<int8_t> out_i8_;
};

static void usage() {
  std::printf(
      "dsr_encoder - Digital Satellite Radio (DSR) encoder\n\n"
      "  -i, --input N=SRC     station N (1..16): file, URL (anything ffmpeg reads),\n"
      "                        tone:FL[,FR] or raw:PATH (s16le 32 kHz stereo)\n"
      "  -n, --name N=TEXT     8-character station name (SK), default \"STATION n\"\n"
      "  -p, --pty N=T[/S]     programme type 0..15 (and optional sub type)\n"
      "  -m, --kind N=music|speech   speech/music flag (default music)\n"
      "  -o, --output OUT      file, '-' (stdout) or tcp://127.0.0.1:PORT (server)\n"
      "  -f, --format F        cs16 (default) | cs8 | cf32 | real16 (IF = fs/4)\n"
      "      --sps K           samples per symbol: 2 (20.48 MS/s, default), 4, 8;\n"
      "                        real16 needs 4 or 8 (default 4)\n"
      "      --realtime / --no-realtime   pace output at real time (default: on for\n"
      "                        tcp and stdout, off for files)\n"
      "      --prefill MS      real-time start-up buffering (default 300)\n"
      "  -t, --duration SEC    stop after SEC seconds\n"
      "      --amp A           rms of one I/Q component, 0..1 of full scale (0.18)\n"
      "      --loop            repeat file sources forever (ffmpeg -stream_loop)\n"
      "      --conj            mirror the spectrum (negate Q). Use if your RF chain/receiver\n"
      "                        needs the opposite rotation sense (see README)\n"
      "      --ffmpeg PATH     ffmpeg binary\n"
      "  -h, --help\n");
}

// Split "N=value".  Returns station index 0..15 or -1.
static int split_n(const char* arg, std::string& val) {
  const char* eq = std::strchr(arg, '=');
  if (!eq) return -1;
  int n = std::atoi(arg);
  if (n < 1 || n > kStations) return -1;
  val = eq + 1;
  return n - 1;
}

int main(int argc, char** argv) {
  signal(SIGPIPE, SIG_IGN);
  signal(SIGINT, on_signal);
  signal(SIGTERM, on_signal);

  Options o;
  std::string names[kStations];
  static option lopts[] = {{"input", 1, 0, 'i'},  {"name", 1, 0, 'n'},   {"pty", 1, 0, 'p'},
                           {"kind", 1, 0, 'm'},   {"output", 1, 0, 'o'}, {"format", 1, 0, 'f'},
                           {"sps", 1, 0, 1000},   {"realtime", 0, 0, 1001}, {"no-realtime", 0, 0, 1002},
                           {"prefill", 1, 0, 1003}, {"duration", 1, 0, 't'}, {"amp", 1, 0, 1004},
                           {"ffmpeg", 1, 0, 1005}, {"conj", 0, 0, 1006}, {"loop", 0, 0, 1007}, {"help", 0, 0, 'h'},  {0, 0, 0, 0}};
  int c, sps_set = 0;
  std::string v;
  while ((c = getopt_long(argc, argv, "i:n:p:m:o:f:t:h", lopts, nullptr)) != -1) {
    int n;
    switch (c) {
      case 'i': if ((n = split_n(optarg, v)) < 0) { usage(); return 1; } o.src[n] = v; o.st[n].active = true; break;
      case 'n': if ((n = split_n(optarg, v)) < 0) { usage(); return 1; } names[n] = v; break;
      case 'p': if ((n = split_n(optarg, v)) < 0) { usage(); return 1; }
                o.st[n].pty = std::atoi(v.c_str()) & 15;
                if (v.find('/') != std::string::npos) o.st[n].sub = std::atoi(v.c_str() + v.find('/') + 1) & 15;
                break;
      case 'm': if ((n = split_n(optarg, v)) < 0) { usage(); return 1; } o.st[n].music = (v != "speech"); break;
      case 'o': o.out = optarg; break;
      case 'f': {
        std::string f = optarg;
        o.fmt = f == "cf32" ? Format::CF32 : f == "real16" ? Format::REAL16 : f == "cs8" ? Format::CS8 : Format::CS16;
        break;
      }
      case 1000: o.sps = std::atoi(optarg); sps_set = 1; break;
      case 1001: o.realtime = o.realtime_set = true; break;
      case 1002: o.realtime = false; o.realtime_set = true; break;
      case 1003: o.prefill_ms = std::atoi(optarg); break;
      case 't': o.duration = std::atof(optarg); break;
      case 1004: o.amp = std::atof(optarg); break;
      case 1005: o.ffmpeg = optarg; break;
      case 1006: o.conj = true; break;
      case 1007: o.loop = true; break;
      default: usage(); return c == 'h' ? 0 : 1;
    }
  }
  int active = 0;
  for (int s = 0; s < kStations; ++s) {
    if (!o.st[s].active) continue;
    ++active;
    std::string nm = names[s].empty() ? "STATION " + std::to_string(s + 1) : names[s];
    sk_from_utf8(nm, o.st[s].name);
  }
  if (!active) { usage(); return 1; }
  if (o.fmt == Format::REAL16 && !sps_set) o.sps = 4;
  if (o.sps != 2 && o.sps != 4 && o.sps != 8) { std::fprintf(stderr, "--sps must be 2, 4 or 8\n"); return 1; }
  if (o.fmt == Format::REAL16 && o.sps < 4) { std::fprintf(stderr, "real16 needs --sps 4 or 8\n"); return 1; }
  if (!o.realtime_set) o.realtime = (o.out == "-" || o.out.rfind("tcp://", 0) == 0);

  std::fprintf(stderr, "DSR encoder: %d station(s), %.2f MS/s, %s, %s\n", active,
               o.sps * kSymbolRate / 1e6, o.fmt == Format::CS16 ? "cs16" : o.fmt == Format::CF32 ? "cf32" : o.fmt == Format::CS8 ? "cs8" : "real16",
               o.realtime ? "real time" : "as fast as possible");
  Encoder enc(o);
  return enc.run();
}
