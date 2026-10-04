// ============================================================================
//  dsr_decoder.cpp  -  DSR (Digital Satellite Radio) decoder / receiver
// ============================================================================
//
//  Baseband I/Q (or real IF) samples  ->  one selected stereo programme.
//
//   samples ─> QPSK demodulator (timing + carrier recovery, dsr_dsp.hpp)
//           ─> hard dibits (A'',B'') ─> differential decoder ─> bit streams A', B'
//           ─> Barker frame synchroniser (also detects swapped A/B = inverted spectrum)
//           ─> descrambler ─> deinterleaver ─> BCH(63,44) error correction
//           ─> Williard super-frame sync from the S bits (64-frame grid)
//           ─> scale factors from the ZI frames (BCH + 3x repetition)
//           ─> 14 -> 16 bit expansion ─> PCM of the selected station
//           ─> sound card (via aplay) and/or WAV file
//
//  While running, type on the keyboard:  1..16 + Enter = tune station,
//  l = list the programme offer, q = quit.
// ============================================================================
#include <getopt.h>
#include <cctype>
#include <signal.h>

#include <atomic>
#include <cstdio>
#include <cstring>
#include <functional>
#include <memory>
#include <string>
#include <thread>

#include "dsr.hpp"
#include "dsr_dsp.hpp"
#include "dsr_io.hpp"

using namespace dsr;

static std::atomic<bool> g_stop{false};
static void on_signal(int) { g_stop = true; }

// ----------------------------------------------------------------------------
// PCM sinks
// ----------------------------------------------------------------------------
class PcmSink {
 public:
  virtual ~PcmSink() {}
  virtual void write(const int16_t* stereo, size_t frames) = 0;
};

// Pipe into an external player; default "aplay" talks to ALSA / PipeWire.
class PlayerSink : public PcmSink {
 public:
  explicit PlayerSink(const std::string& cmd) {
    fd_ = io::spawn_writer(cmd, &pid_);
    if (fd_ < 0) std::fprintf(stderr, "cannot start player\n");
  }
  ~PlayerSink() override {
    if (fd_ >= 0) close(fd_);
    if (pid_ > 0) waitpid(pid_, nullptr, 0);
  }
  void write(const int16_t* s, size_t n) override {
    if (fd_ >= 0 && !io::write_all(fd_, s, n * 4)) { close(fd_); fd_ = -1; g_stop = true; }
  }
 private:
  int fd_ = -1;
  pid_t pid_ = -1;
};

// 16 bit stereo 32 kHz RIFF/WAVE file; the header is patched on close.
class WavSink : public PcmSink {
 public:
  explicit WavSink(const std::string& path) {
    f_ = std::fopen(path.c_str(), "wb");
    if (!f_) { std::perror("wav"); return; }
    uint8_t h[44] = {0};
    std::memcpy(h, "RIFF", 4); std::memcpy(h + 8, "WAVEfmt ", 8);
    put32(h + 16, 16); put16(h + 20, 1); put16(h + 22, 2);
    put32(h + 24, kAudioRate); put32(h + 28, kAudioRate * 4); put16(h + 32, 4); put16(h + 34, 16);
    std::memcpy(h + 36, "data", 4);
    std::fwrite(h, 1, 44, f_);
  }
  ~WavSink() override {
    if (!f_) return;
    uint8_t b[4];
    put32(b, uint32_t(36 + bytes_)); std::fseek(f_, 4, SEEK_SET); std::fwrite(b, 1, 4, f_);
    put32(b, uint32_t(bytes_));      std::fseek(f_, 40, SEEK_SET); std::fwrite(b, 1, 4, f_);
    std::fclose(f_);
  }
  void write(const int16_t* s, size_t n) override {
    if (f_) { std::fwrite(s, 4, n, f_); bytes_ += n * 4; }
  }
 private:
  static void put16(uint8_t* p, uint32_t v) { p[0] = v; p[1] = v >> 8; }
  static void put32(uint8_t* p, uint32_t v) { p[0] = v; p[1] = v >> 8; p[2] = v >> 16; p[3] = v >> 24; }
  FILE* f_ = nullptr;
  uint64_t bytes_ = 0;
};

// ----------------------------------------------------------------------------
// FrameSync: finds the 320-bit main frame boundary with the Barker words.
// Input: one bit of each of the two streams per call.  Frame A carries the
// Barker word, frame B its inverse; if the streams are swapped (spectrum in
// inverted position, spec 4.5 note 2) this is detected and corrected.
// Lock criterion: 3 consecutive frames with a perfect match at the same offset.
// Unlock: 8 consecutive frames with more than 4 wrong sync bits.
// ----------------------------------------------------------------------------
class FrameSync {
 public:
  using Callback = std::function<void(const u8* a, const u8* b)>;
  explicit FrameSync(Callback cb) : cb_(std::move(cb)) { reset_search(); }

  bool locked() const { return locked_; }
  bool swapped() const { return swap_; }

  void push(u8 x, u8 y) {
    if (!locked_) {
      sx_ = ((sx_ << 1) | x) & 0x7FF;
      sy_ = ((sy_ << 1) | y) & 0x7FF;
      if (t_ >= 10) {
        int h = -1;
        if (sx_ == kA11() && sy_ == kB11()) h = 0;
        else if (sx_ == kB11() && sy_ == kA11()) h = 1;
        if (h >= 0) {
          int o = static_cast<int>(t_ % kFrameBits);
          int64_t fr = static_cast<int64_t>(t_ / kFrameBits);
          cnt_[h][o] = (last_[h][o] + 1 == fr) ? cnt_[h][o] + 1 : 1;
          last_[h][o] = fr;
          if (cnt_[h][o] >= 3) {            // lock: current bit is bit 10 of the frame
            locked_ = true; swap_ = (h == 1); pos_ = 11; miss_ = 0;
            std::memcpy(fa_, kBarkerA, 11); std::memcpy(fb_, kBarkerB, 11);
          }
        }
      }
      ++t_;
      return;
    }
    fa_[pos_] = swap_ ? y : x;
    fb_[pos_] = swap_ ? x : y;
    if (++pos_ == kFrameBits) {
      int err = 0;
      for (int i = 0; i < kSyncBits; ++i) err += (fa_[i] != kBarkerA[i]) + (fb_[i] != kBarkerB[i]);
      miss_ = (err <= 4) ? 0 : miss_ + 1;
      cb_(fa_, fb_);
      pos_ = 0;
      if (miss_ >= 8) { locked_ = false; reset_search(); }
    }
    ++t_;
  }

 private:
  static u32 kA11() { u32 v = 0; for (int i = 0; i < 11; ++i) v = (v << 1) | kBarkerA[i]; return v; }
  static u32 kB11() { return ~kA11() & 0x7FF; }
  void reset_search() {
    sx_ = sy_ = 0;
    for (auto& r : cnt_) for (auto& c : r) c = 0;
    for (auto& r : last_) for (auto& c : r) c = -10;
  }
  Callback cb_;
  bool locked_ = false, swap_ = false;
  uint64_t t_ = 0;
  u32 sx_ = 0, sy_ = 0;
  int cnt_[2][kFrameBits];
  int64_t last_[2][kFrameBits];
  int pos_ = 0, miss_ = 0;
  u8 fa_[kFrameBits], fb_[kFrameBits];
};

// ----------------------------------------------------------------------------
// Receiver: everything above the bit layer.
// ----------------------------------------------------------------------------
class Receiver {
 public:
  Receiver(PcmSink* sink, int station, bool verbose)
      : sink_(sink), sel_(station), verbose_(verbose),
        sync_([this](const u8* a, const u8* b) { on_frame(a, b); }) {
    std::memset(skv_, 0, sizeof skv_);
    for (auto& s : info_) std::memset(s.name, ' ', 8);
    cnt_reset();
  }

  int  selected() const { return sel_; }

  // Feed hard dibits of the demodulator.
  void push_dibits(const u8* d, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      u8 A1, B1;
      diff_.step(d[i] >> 1, d[i] & 1, A1, B1);
      sync_.push(A1, B1);
    }
  }

  // Programme offer as text (used by the keyboard console and the control port).
  std::string offer_text() const {
    char buf[160];
    std::string t;
    std::snprintf(buf, sizeof buf, "Programme offer (frame lock %s%s, grid lock %s):\n", sync_.locked() ? "yes" : "no",
                  sync_.swapped() ? ", A/B swapped" : "", grid_ ? "yes" : "no");
    t += buf;
    int n = 0;
    for (int s = 0; s < kStations; ++s) {
      if (!pa_ok_[s]) continue;               // only channels announced as in use
      ++n;
      std::snprintf(buf, sizeof buf, "  %2d %c %-8s  %-8s  %s%s\n", s + 1, s == sel_ - 1 ? '*' : ' ',
                    name_of(s).c_str(), kPtyNames[info_[s].pty], info_[s].music ? "music" : "speech",
                    mono_[s] ? "  (2 x mono: L = programme A, R = programme B)" : "");
      t += buf;
    }
    if (!n) t += "  (no programme announced yet - wait a few seconds)\n";
    return t;
  }

  std::string name_of(int s) const {
    std::string nm;
    for (int i = 0; i < 8; ++i) nm += sk_to_utf8(info_[s].name[i]);
    return nm;
  }

  // Tune to station s (1..16).  Returns a one-line confirmation.
  std::string tune(int s) {
    if (s < 1 || s > kStations) return "station must be 1..16\n";
    sel_ = s;
    std::string r = "station " + std::to_string(s);
    if (pa_ok_[s - 1]) r += ": " + name_of(s - 1);
    else r += " (not announced in the programme offer)";
    return r + "\n";
  }

  // Step to the next/previous announced station.
  std::string step(int dir) {
    int s = sel_;
    for (int i = 0; i < kStations; ++i) {
      s = (s - 1 + dir + kStations) % kStations + 1;
      if (pa_ok_[s - 1]) return tune(s);
    }
    return "no announced stations yet\n";
  }

  void stats() const {
    std::fprintf(stderr, "[frames %llu  bch-corrected %llu  uncorrectable blocks %llu  audio blocks %llu  pi packets %llu]\n",
                 (unsigned long long)frames_, (unsigned long long)corrected_, (unsigned long long)uncorr_,
                 (unsigned long long)audio_blocks_, (unsigned long long)pi_packets_);
  }
  bool locked() const { return sync_.locked() && grid_; }

 private:
  // ---- called for every aligned pair of 320-bit frames ----
  void on_frame(const u8* a_in, const u8* b_in) {
    ++frames_;
    u8 fa[kFrameBits], fb[kFrameBits];
    std::memcpy(fa, a_in, kFrameBits);
    std::memcpy(fb, b_in, kFrameBits);
    scramble_frame(false, fa);                       // descramble = scramble
    scramble_frame(true, fb);
    u8 sbit = fa[kSbitIndex];

    // --- Williard grid search on the S bits (frame A only) ---
    s16_ = static_cast<u16>((s16_ << 1) | sbit);
    if (!grid_) {
      if (frames_ >= 16 && (s16_ == kSync1 || s16_ == kSync2)) {
        int o = static_cast<int>(frames_ % 64);
        cnt_[o] = (gl_[o] + 64 == (int64_t)frames_) ? cnt_[o] + 1 : 1;
        gl_[o] = static_cast<int64_t>(frames_);
        if (cnt_[o] >= 3) {            // this frame carries the last sync bit -> next frame starts a block
          grid_ = true; q_ = 0; blk_ = 0; sync_miss_ = 0; sa_in_sau_ = -1; sau_pos_ = -1;
          std::memset(skv_, 0, sizeof skv_);
        }
      }
      return;
    }

    // --- demultiplex this frame (grid is locked, q_ = position in the 64-frame block) ---
    StationSlot slotA[8], slotB[8];
    bool badA[8], badB[8];
    int ca = parse_frame(fa, slotA, badA);
    int cb = parse_frame(fb, slotB, badB);
    corrected_ += ca + cb;
    for (int i = 0; i < 8; ++i) {
      store(i, q_, slotA[i], badA[i]);
      store(8 + i, q_, slotB[i], badB[i]);
      uncorr_ += badA[i] + badB[i];
    }

    // --- special-service frame: bit index of this frame inside the SA frame ---
    int sa_bit = (q_ + 16) % 64;
    if (sa_bit < 16) {
      sync_acc_ = static_cast<u16>((sync_acc_ << 1) | sbit);
      if (sa_bit == 15) end_sync_word();
    } else {
      int by = (sa_bit - 16) / 8;
      sa_[by] = static_cast<u8>((sa_[by] << 1) | sbit);
      if (sa_bit == 63) end_sa_frame();
    }
    if (!grid_) return;

    if (q_ == 63) end_block();
    q_ = (q_ + 1) % 64;
  }

  void store(int s, int q, const StationSlot& sl, bool bad) {
    cw_[s][0][q] = sl.cw[0];
    cw_[s][1][q] = sl.cw[1];
    zi_[s][q] = sl.zi;
    bad_[s][q] = bad;
  }

  // 16 sync bits of the next SA frame are complete.
  void end_sync_word() {
    int e1 = __builtin_popcount(sync_acc_ ^ kSync1), e2 = __builtin_popcount(sync_acc_ ^ kSync2);
    if (e1 <= 2 && e1 <= e2) { sa_in_sau_ = 0; sync_miss_ = 0; }
    else if (e2 <= 2)        { sa_in_sau_ = (sa_in_sau_ < 0) ? -1 : (sa_in_sau_ + 1) % 8; sync_miss_ = 0; }
    else if (++sync_miss_ >= 4) { grid_ = false; cnt_reset(); }   // lost the block grid
    else sa_in_sau_ = (sa_in_sau_ < 0) ? -1 : (sa_in_sau_ + 1) % 8;
  }

  // 48 content bits (6 bytes) of an SA frame are complete.
  void end_sa_frame() {
    if (sa_in_sau_ < 0) return;
    std::memcpy(sau_[sa_in_sau_], sa_, 6);
    if (sa_in_sau_ == 7) end_sau();
  }

  // 8 SA frames (= one SAUE) complete.  Identify its type inside the 16-part SAUEUE.
  void end_sau() {
    bool zero = true;
    for (auto& f : sau_) for (int i = 0; i < 6; ++i) zero &= (f[i] == 0);
    if (zero) sau_pos_ = 7;                         // LB frame: reference point of the SAUEUE
    else if (sau_pos_ >= 0) sau_pos_ = (sau_pos_ + 1) % 16;
    if (sau_pos_ < 0) return;
    if (sau_pos_ < 7) {                             // SAUE/PA
      for (int k = 0; k < 8; ++k) parse_pa(2 * k, sau_[k]);
    } else if (sau_pos_ >= 8) {                     // SAUE/SK, character index sau_pos_-8
      int c = sau_pos_ - 8;
      for (int k = 0; k < 8; ++k) {
        info_[2 * k].name[c] = sau_[k][0];
        info_[2 * k + 1].name[c] = sau_[k][2];
        seen_[2 * k] = seen_[2 * k + 1] = true;
      }
    }
  }

  void parse_pa(int s0, const u8 b[6]) {
    for (int j = 0; j < 2; ++j) {
      u8 l = b[2 * j], r = b[2 * j + 1];
      int s = s0 + j;
      if (l == kPaUnused) { pa_ok_[s] = false; continue; }
      if (with_parity(l) != l || with_parity(r) != r) continue;       // parity check
      int ml = (l >> 1) & 3, mr = (r >> 1) & 3;
      if (ml == 1 && mr == 2) mono_[s] = false;                        // stereo pair
      else if (ml == 1 && mr == 1) mono_[s] = true;                    // two independent mono programmes
      else continue;
      info_[s].pty = l >> 4;
      info_[s].music = (l >> 3) & 1;
      info_[s].sub = r >> 4;
      pa_ok_[s] = true;
      seen_[s] = true;
    }
  }

  // One 2 ms block (64 frames) complete: output audio, then decode the ZI frames.
  void end_block() {
    // --- audio of the selected station; its scale factors arrived 4 ms ago ---
    int sel = sel_ - 1;
    if (sel >= 0 && sel < kStations) {
      const SkEntry& e = skv_[blk_ % 4][sel];
      if (e.valid) {
        int16_t out[2 * kBlockSamples];
        for (int i = 0; i < kBlockSamples; ++i)
          for (int c = 0; c < 2; ++c) {
            // concealment (spec 3.2.2): an uncorrectable BCH word -> hold the last sample
            if (bad_[sel][i]) held_[c] = held_[c];
            else held_[c] = decode_sample(cw_[sel][c][i], c == 0 ? e.l : e.r);
            out[2 * i + c] = held_[c];
          }
        sink_->write(out, kBlockSamples);
        ++audio_blocks_;
      }
    }
    // --- the ZI frame of this block carries the scale factors of block blk_+2 ---
    for (int s = 0; s < kStations; ++s) skv_[blk_ % 4][s].valid = false;
    for (int s = 0; s < kStations; ++s) {
      SkEntry& e = skv_[(blk_ + 2) % 4][s];
      int l = 0, r = 0;
      e.valid = decode_zi(zi_[s], l, r);
      e.l = static_cast<u8>(l);
      e.r = static_cast<u8>(r);
    }
    if (sel >= 0 && sel < kStations) feed_pi(zi_[sel] + 42);
    ++blk_;
    if (verbose_ && blk_ % 2500 == 0) stats();
  }

  struct SkEntry { bool valid; u8 l, r; };

  // 64 ZI bits -> two scale factors.  Three BCH-protected copies: take the
  // first decodable pair that at least two copies agree on, else any decodable one.
  static bool decode_zi(const u8* zi, int& l, int& r) {
    int L[3], R[3];
    bool ok[3];
    for (int i = 0; i < 3; ++i) ok[i] = zi_decode(zi + 14 * i, L[i], R[i]);
    for (int i = 0; i < 3; ++i)
      for (int j = i + 1; j < 3; ++j)
        if (ok[i] && ok[j] && L[i] == L[j] && R[i] == R[j]) { l = L[i]; r = R[i]; return true; }
    for (int i = 0; i < 3; ++i) if (ok[i]) { l = L[i]; r = R[i]; return true; }
    u8 m[14];                                                        // last resort: bit-wise majority
    for (int b = 0; b < 14; ++b) m[b] = (zi[b] + zi[14 + b] + zi[28 + b]) >= 2;
    return zi_decode(m, l, r);
  }

  // PI: scan the 22-bit words for packet headers (start word + Hamming-coded bytes).
  void feed_pi(const u8* w) {
    static const u8 start[12] = {0,0,0,0,0,0,1,1,1,1,1,1};
    if (pi_left_ > 0) { --pi_left_; return; }          // content words of the current packet
    if (pi_wait_) { pi_wait_ = false; finish_pi_header(w); return; }
    if (std::memcmp(w, start, 12) == 0) { pi_cache_.assign(w, w + 22); pi_wait_ = true; }
  }
  void finish_pi_header(const u8* second) {
    std::vector<u8> h = pi_cache_;
    h.insert(h.end(), second, second + 22);
    auto byte_at = [&](int off) { int v = 0; for (int i = 0; i < 8; ++i) v = (v << 1) | h[off + i]; return v; };
    int n1 = hamming84_decode(byte_at(12)), n2 = hamming84_decode(byte_at(20));
    int c1 = hamming84_decode(byte_at(28)), c2 = hamming84_decode(byte_at(36));
    if (n1 < 0 || n2 < 0 || c1 < 0 || c2 < 0) return;
    pi_left_ = (n1 << 4) | n2;
    ++pi_packets_;
    if (verbose_ && pi_packets_ <= 3)
      std::fprintf(stderr, "[PI packet: content id %d, %d word(s)]\n", (c1 << 4) | c2, pi_left_);
  }

  void cnt_reset() { for (auto& c : cnt_) c = 0; for (auto& l : gl_) l = -1000; }

  PcmSink* sink_;
  std::atomic<int> sel_;
  bool verbose_;
  DiffDecoder diff_;
  FrameSync sync_;

  // grid / block state
  bool grid_ = false;
  int q_ = 0;
  uint64_t blk_ = 0;
  u16 s16_ = 0, sync_acc_ = 0;
  int cnt_[64] = {0};
  int64_t gl_[64];
  int sync_miss_ = 0;
  uint64_t frames_ = 0, corrected_ = 0, uncorr_ = 0, audio_blocks_ = 0, pi_packets_ = 0;

  // per-block buffers (all 16 stations, so tuning is instantaneous)
  u16 cw_[kStations][2][64];
  u8 zi_[kStations][64];
  bool bad_[kStations][64];
  SkEntry skv_[4][kStations];
  int16_t held_[2] = {0, 0};

  // special service
  u8 sa_[6] = {0};
  u8 sau_[8][6];
  int sa_in_sau_ = -1, sau_pos_ = -1;
  StationInfo info_[kStations];
  bool seen_[kStations] = {false};
  bool pa_ok_[kStations] = {false};
  bool mono_[kStations] = {false};

  // PI
  std::vector<u8> pi_cache_;
  bool pi_wait_ = false;
  int pi_left_ = 0;
};

// ----------------------------------------------------------------------------
// Run-time control.  The same command set is available on the keyboard and on a
// local TCP "control port" (line based, e.g. `nc 127.0.0.1 6000`), so scripts or
// a GUI can switch programmes and fetch the programme list at any time.
//   <1..16> | s <n> : tune to station n        l | list : programme list
//   + / n : next    - / p : previous           h | ? : help      q : quit
// ----------------------------------------------------------------------------
static const char* kHelp =
    "commands:  1..16 = tune station | l = programme list | + / - = next / previous\n"
    "           h = help | q = quit\n";

static std::string handle_command(Receiver& rx, std::string line) {
  while (!line.empty() && (line.back() == '\n' || line.back() == '\r' || line.back() == ' ')) line.pop_back();
  size_t b = line.find_first_not_of(' ');
  if (b == std::string::npos) return "";
  line = line.substr(b);
  if (line == "q" || line == "quit") { g_stop = true; return "bye\n"; }
  if (line == "l" || line == "list") return rx.offer_text();
  if (line == "h" || line == "?" || line == "help") return kHelp;
  if (line == "+" || line == "n" || line == "next") return rx.step(+1);
  if (line == "-" || line == "p" || line == "prev") return rx.step(-1);
  const char* p = line.c_str();
  if (line[0] == 's' && line.size() > 1 && line.find_first_of("0123456789") != std::string::npos) p = line.c_str() + 1;
  if (std::isdigit(static_cast<unsigned char>(*p == ' ' ? *(p + 1) : *p))) return rx.tune(std::atoi(p));
  return std::string("unknown command '") + line + "'\n" + kHelp;
}

// Keyboard console (only meaningful if stdin is not the sample source).
static void console_loop(Receiver* rx) {
  char line[128];
  while (!g_stop && std::fgets(line, sizeof line, stdin)) {
    std::string r = handle_command(*rx, line);
    std::fputs(r.c_str(), stderr);
    std::fflush(stderr);
  }
}

// TCP control port: any number of clients, each served by its own thread.
static void control_server(Receiver* rx, int port) {
  int s = socket(AF_INET, SOCK_STREAM, 0);
  int one = 1;
  setsockopt(s, SOL_SOCKET, SO_REUSEADDR, &one, sizeof one);
  sockaddr_in a{};
  a.sin_family = AF_INET;
  a.sin_port = htons(static_cast<uint16_t>(port));
  inet_pton(AF_INET, "127.0.0.1", &a.sin_addr);
  if (bind(s, reinterpret_cast<sockaddr*>(&a), sizeof a) < 0 || listen(s, 4) < 0) { perror("control port"); return; }
  while (!g_stop) {
    int c = accept(s, nullptr, nullptr);
    if (c < 0) continue;
    std::thread([rx, c] {
      io::write_all(c, kHelp, std::strlen(kHelp));
      std::string buf;
      char tmp[256];
      ssize_t n;
      while (!g_stop && (n = io::read_some(c, tmp, sizeof tmp)) > 0) {
        buf.append(tmp, static_cast<size_t>(n));
        size_t nl;
        while ((nl = buf.find('\n')) != std::string::npos) {
          std::string r = handle_command(*rx, buf.substr(0, nl));
          io::write_all(c, r.c_str(), r.size());
          buf.erase(0, nl + 1);
        }
      }
      close(c);
    }).detach();
  }
}

static void usage() {
  std::printf(
      "dsr_decoder - Digital Satellite Radio (DSR) decoder\n\n"
      "  -i, --input SRC     file, '-' (stdin) or tcp://HOST:PORT\n"
      "  -f, --format F      cs16 (default) | cf32 | real16\n"
      "      --sps K         samples per symbol of the input: 2 (default), 4, 8\n"
      "  -s, --station N     programme to decode (1..16, default 1; changeable at run time)\n"
      "  -o, --wav FILE      write the programme to a 32 kHz stereo WAV file\n"
      "      --play [CMD]    play via CMD (default: aplay -q -t raw -f S16_LE -r 32000 -c 2 -)\n"
      "      --list          print the programme offer when the input ends\n"
      "      --control PORT  command port on 127.0.0.1 (line protocol, see below)\n"
      "\nRun-time commands (keyboard + Enter, or via the control port):\n"
      "  1..16 tune | l programme list | + / - next / previous | h help | q quit\n"
      "  -v, --verbose\n"
      "  -h, --help\n");
}

int main(int argc, char** argv) {
  signal(SIGPIPE, SIG_IGN);
  signal(SIGINT, on_signal);
  signal(SIGTERM, on_signal);

  std::string in = "dsr.cs16", wav, player;
  bool play = false, list = false, verbose = false;
  int control_port = 0;
  int fmt = 0, sps = 2, station = 1;       // fmt: 0 cs16, 1 cf32, 2 real16
  static option lo[] = {{"input", 1, 0, 'i'}, {"format", 1, 0, 'f'}, {"sps", 1, 0, 1000},
                        {"station", 1, 0, 's'}, {"wav", 1, 0, 'o'}, {"play", 2, 0, 1001},
                        {"list", 0, 0, 1002}, {"control", 1, 0, 1003}, {"verbose", 0, 0, 'v'}, {"help", 0, 0, 'h'}, {0, 0, 0, 0}};
  int c;
  bool sps_set = false;
  while ((c = getopt_long(argc, argv, "i:f:s:o:vh", lo, nullptr)) != -1) {
    switch (c) {
      case 'i': in = optarg; break;
      case 'f': fmt = !strcmp(optarg, "cf32") ? 1 : !strcmp(optarg, "real16") ? 2 : 0; break;
      case 1000: sps = std::atoi(optarg); sps_set = true; break;
      case 's': station = std::atoi(optarg); break;
      case 'o': wav = optarg; break;
      case 1001: play = true; if (optarg) player = optarg; break;
      case 1002: list = true; break;
      case 1003: control_port = std::atoi(optarg); break;
      case 'v': verbose = true; break;
      default: usage(); return c == 'h' ? 0 : 1;
    }
  }
  if (fmt == 2 && !sps_set) sps = 4;
  if (sps != 2 && sps != 4 && sps != 8) { std::fprintf(stderr, "--sps must be 2, 4 or 8\n"); return 1; }
  if (wav.empty()) play = true;
  if (player.empty()) player = "aplay -q -t raw -f S16_LE -r 32000 -c 2 -";

  // ---- input ----
  int fd;
  bool is_tcp = false;
  std::string tcp_host;
  int tcp_port = 0;
  if (in == "-") fd = 0;
  else if (in.rfind("tcp://", 0) == 0) {
    std::string hp = in.substr(6);
    size_t col = hp.rfind(':');
    is_tcp = true;
    tcp_host = hp.substr(0, col);
    tcp_port = std::atoi(hp.c_str() + col + 1);
    fd = io::tcp_connect(tcp_host, tcp_port);
  } else fd = ::open(in.c_str(), O_RDONLY);
  if (fd < 0) { std::perror("input"); return 1; }

  // ---- output ----
  struct Tee : PcmSink {
    std::vector<std::unique_ptr<PcmSink>> v;
    void write(const int16_t* s, size_t n) override { for (auto& k : v) k->write(s, n); }
  } tee;
  if (!wav.empty()) tee.v.emplace_back(new WavSink(wav));
  if (play) tee.v.emplace_back(new PlayerSink(player));

  Receiver rx(&tee, station, verbose);
  Demodulator demod(sps, fmt == 2);

  // ---- run-time control: keyboard (unless stdin carries samples) and optional port ----
  if (fd != 0) {
    std::thread(console_loop, &rx).detach();
    if (isatty(0)) std::fputs(kHelp, stderr);
  }
  if (control_port) {
    std::thread(control_server, &rx, control_port).detach();
    std::fprintf(stderr, "control port: 127.0.0.1:%d\n", control_port);
  }

  // ---- main loop ----
  const size_t kChunk = 1 << 16;                     // samples per read
  size_t bps = fmt == 1 ? 8 : fmt == 2 ? 2 : 4;      // bytes per sample
  std::vector<uint8_t> raw(kChunk * bps);
  std::vector<cf> cx(kChunk);
  std::vector<u8> dib;
  size_t carry = 0;
  while (!g_stop) {
    ssize_t r = io::read_some(fd, raw.data() + carry, raw.size() - carry);
    if (r <= 0) {
      if (!is_tcp || g_stop) break;
      // The encoder went away (or was restarted): keep the decoder alive and reconnect.
      std::fprintf(stderr, "connection lost - reconnecting ...\n");
      close(fd);
      fd = io::tcp_connect(tcp_host, tcp_port, 3600);
      if (fd < 0) break;
      std::fprintf(stderr, "reconnected\n");
      carry = 0;
      continue;
    }
    size_t have = carry + static_cast<size_t>(r), n = have / bps;
    for (size_t k = 0; k < n; ++k) {
      const uint8_t* p = &raw[k * bps];
      if (fmt == 0)      { int16_t a, b; std::memcpy(&a, p, 2); std::memcpy(&b, p + 2, 2); cx[k] = {a / 32768.f, b / 32768.f}; }
      else if (fmt == 1) { float a, b; std::memcpy(&a, p, 4); std::memcpy(&b, p + 4, 4); cx[k] = {a, b}; }
      else               { int16_t a; std::memcpy(&a, p, 2); cx[k] = {a / 32768.f, 0.f}; }
    }
    carry = have - n * bps;
    std::memmove(raw.data(), raw.data() + n * bps, carry);
    dib.clear();
    demod.process(cx.data(), n, dib);
    rx.push_dibits(dib.data(), dib.size());
  }
  if (list) std::fputs(rx.offer_text().c_str(), stdout);
  rx.stats();
  // The control thread may still sit in fgets(stdin) and hold the stdio lock that
  // exit() needs -> close the sinks explicitly (patches the WAV header) and _exit.
  tee.v.clear();
  std::fflush(stdout);
  std::fflush(stderr);
  _exit(0);
}
