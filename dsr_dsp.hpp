// ============================================================================
//  dsr_dsp.hpp  -  Baseband signal processing for the DSR QPSK signal
// ============================================================================
//
//  Transmit side : symbol (+-1, +-1) -> root-raised-cosine pulse shaping
//  Receive side  : [real IF -> complex mix] -> decimating low-pass -> matched
//                  filter -> AGC -> Gardner timing recovery (cubic interpolator)
//                  -> decision-directed Costas loop -> hard decisions
//
//  Sample rates: the signal is always processed with an integer number of
//  samples per symbol (sps).  Symbol rate = 10.24 MBd, so sps=2 -> 20.48 MS/s.
// ============================================================================
#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "dsr.hpp"

namespace dsr {

struct cf { float r, i; };   // complex float sample

// ----------------------------------------------------------------------------
// Root-raised-cosine impulse response (spec 4.4: 50 % cos roll-off).
// The spec gives the *frequency* response S(f) = sqrt(raised cosine), which is
// exactly the RRC filter.  Closed-form time-domain expression, T = 1 symbol.
// ----------------------------------------------------------------------------
inline std::vector<float> rrc_taps(int sps, int span_symbols, double beta = kRollOff) {
  const double pi = M_PI;
  int n = span_symbols * sps + 1;
  std::vector<float> h(n);
  for (int k = 0; k < n; ++k) {
    double t = (k - (n - 1) / 2.0) / sps;      // time in symbol periods
    double v;
    if (std::fabs(t) < 1e-9) {
      v = 1.0 - beta + 4.0 * beta / pi;
    } else if (std::fabs(std::fabs(t) - 1.0 / (4.0 * beta)) < 1e-9) {
      v = beta / std::sqrt(2.0) *
          ((1 + 2 / pi) * std::sin(pi / (4 * beta)) + (1 - 2 / pi) * std::cos(pi / (4 * beta)));
    } else {
      v = (std::sin(pi * t * (1 - beta)) + 4 * beta * t * std::cos(pi * t * (1 + beta))) /
          (pi * t * (1 - (4 * beta * t) * (4 * beta * t)));
    }
    h[k] = static_cast<float>(v);
  }
  return h;
}

// Hamming-windowed sinc low-pass, cutoff fc in cycles/sample (0..0.5), unity DC gain.
inline std::vector<float> lowpass_taps(int n, double fc) {
  std::vector<float> h(n);
  double sum = 0;
  for (int k = 0; k < n; ++k) {
    double m = k - (n - 1) / 2.0;
    double s = (std::fabs(m) < 1e-12) ? 2 * fc : std::sin(2 * M_PI * fc * m) / (M_PI * m);
    double w = 0.54 - 0.46 * std::cos(2 * M_PI * k / (n - 1));
    h[k] = static_cast<float>(s * w);
    sum += h[k];
  }
  for (auto& x : h) x = static_cast<float>(x / sum);
  return h;
}

// ----------------------------------------------------------------------------
// Pulse shaper: zero-stuffing + FIR implemented as a polyphase structure.
// For every symbol it produces `sps` output samples.
// ----------------------------------------------------------------------------
class PulseShaper {
 public:
  explicit PulseShaper(int sps, int span = 16) : sps_(sps) {
    std::vector<float> h = rrc_taps(sps, span);
    int n = static_cast<int>(h.size());
    L_ = (n + sps - 1) / sps;                           // taps per polyphase branch
    ph_.assign(static_cast<size_t>(sps) * L_, 0.f);
    double e = 0;
    for (int p = 0; p < sps; ++p)
      for (int j = 0; j < L_; ++j) {
        int idx = p + j * sps;
        if (idx < n) ph_[p * L_ + j] = h[idx];
      }
    for (float x : h) e += double(x) * x;
    power_ = e / sps;                                   // mean power per sample for unit symbols
    hi_.assign(2 * L_, 0.f);
    hq_.assign(2 * L_, 0.f);
  }
  double unit_power() const { return power_; }

  void push(float I, float Q, float* oi, float* oq) {
    pos_ = (pos_ == 0) ? L_ - 1 : pos_ - 1;             // newest symbol at hist[pos_]
    hi_[pos_] = hi_[pos_ + L_] = I;                     // duplicated -> contiguous window
    hq_[pos_] = hq_[pos_ + L_] = Q;
    for (int p = 0; p < sps_; ++p) {
      const float* h = &ph_[p * L_];
      float ai = 0, aq = 0;
      for (int j = 0; j < L_; ++j) { ai += h[j] * hi_[pos_ + j]; aq += h[j] * hq_[pos_ + j]; }
      oi[p] = ai;
      oq[p] = aq;
    }
  }

 private:
  int sps_, L_, pos_ = 0;
  double power_;
  std::vector<float> ph_, hi_, hq_;
};

// ----------------------------------------------------------------------------
// Streaming FIR filter with optional decimation (complex samples, real taps).
// ----------------------------------------------------------------------------
class FirFilter {
 public:
  FirFilter(std::vector<float> taps, int decim = 1) : h_(std::move(taps)), D_(decim) {}
  // Appends outputs to `out`.
  void process(const cf* x, size_t n, std::vector<cf>& out) {
    buf_.insert(buf_.end(), x, x + n);
    const size_t N = h_.size();
    size_t pos = 0;
    while (pos + N <= buf_.size()) {
      float ar = 0, ai = 0;
      const cf* b = &buf_[pos];
      for (size_t k = 0; k < N; ++k) { ar += b[k].r * h_[k]; ai += b[k].i * h_[k]; }
      out.push_back({ar, ai});
      pos += D_;
    }
    buf_.erase(buf_.begin(), buf_.begin() + std::min(pos, buf_.size()));
  }
 private:
  std::vector<float> h_;
  int D_;
  std::vector<cf> buf_;
};

// ----------------------------------------------------------------------------
// QPSK demodulator: complex/real samples in -> hard dibits out.
// Output byte = (A'' << 1) | B''  (bit definitions: see dsr_encoder modulator).
// ----------------------------------------------------------------------------
class Demodulator {
 public:
  Demodulator(int sps_in, bool real_input)
      : real_(real_input),
        D_(sps_in / 2),
        mf_(rrc_taps(2, 16)),
        dec_(std::vector<float>{1.f}, 1) {
    // Cut-off halfway between the signal edge (7.68 MHz) and the first alias
    // (fs_out - 7.68 MHz = 12.8 MHz) -> 10.24 MHz, normalised to the input rate.
    if (D_ > 1) dec_ = FirFilter(lowpass_taps(15 * D_ + 1, 10.24e6 / (sps_in * 10.24e6)), D_);
  }

  void process(const cf* in, size_t n, std::vector<u8>& out) {
    // --- 1. front end: real IF at fs/4 -> complex baseband (multiply by e^{-j pi n/2}) ---
    tmp_.clear();
    if (real_) {
      for (size_t k = 0; k < n; ++k) {
        float r = in[k].r;
        switch (mix_) {
          case 0: tmp_.push_back({ r, 0}); break;
          case 1: tmp_.push_back({ 0, -r}); break;
          case 2: tmp_.push_back({-r, 0}); break;
          default: tmp_.push_back({ 0,  r}); break;
        }
        mix_ = (mix_ + 1) & 3;
      }
    } else {
      tmp_.assign(in, in + n);
    }
    // --- 2. low-pass + decimation to 2 samples/symbol ---
    s2_.clear();
    if (D_ > 1) dec_.process(tmp_.data(), tmp_.size(), s2_); else s2_.swap(tmp_);
    // --- 3. matched filter (RRC) at 2 sps ---
    mfo_.clear();
    mf_.process(s2_.data(), s2_.size(), mfo_);
    // --- 4. AGC: leaky power estimate, target mean power kAgcTarget ---
    for (cf s : mfo_) {
      pwr_ += 1e-4f * ((s.r * s.r + s.i * s.i) - pwr_);
      float g = std::sqrt(kAgcTarget / (pwr_ + 1e-12f));
      y_.push_back({s.r * g, s.i * g});
    }
    // --- 5. symbol timing + carrier recovery ---
    run_loops(out);
    // drop consumed samples (keep a few for the interpolator)
    long keep_from = static_cast<long>(std::floor(pos_)) - 4;
    if (keep_from - base_ > 8192) {
      y_.erase(y_.begin(), y_.begin() + (keep_from - base_));
      base_ = keep_from;
    }
  }

 private:
  static constexpr float kAgcTarget = 1.5f;

  // 4-point Lagrange (cubic) interpolation at absolute position t.
  cf interp(double t) const {
    long i = static_cast<long>(std::floor(t));
    float mu = static_cast<float>(t - i);
    const cf* p = &y_[i - base_ - 1];
    float cm1 = -mu * (mu - 1) * (mu - 2) / 6.f;
    float c0  = (mu + 1) * (mu - 1) * (mu - 2) / 2.f;
    float c1  = -(mu + 1) * mu * (mu - 2) / 2.f;
    float c2  = (mu + 1) * mu * (mu - 1) / 6.f;
    return {cm1 * p[0].r + c0 * p[1].r + c1 * p[2].r + c2 * p[3].r,
            cm1 * p[0].i + c0 * p[1].i + c1 * p[2].i + c2 * p[3].i};
  }

  void run_loops(std::vector<u8>& out) {
    // Gardner timing-error detector, 2 samples per symbol:
    //   e = Re{ (y_n - y_{n-1}) * conj(y_{n-1/2}) }
    // positive e = sampling too late -> advance strobe later by less than 2 samples.
    constexpr float kTp = 0.02f, kTi = 0.0001f;       // timing PI loop
    constexpr float kCa = 0.03f, kCb = 2.2e-4f;       // Costas PI loop (phase, frequency)
    while (static_cast<long>(std::floor(pos_)) + 2 < base_ + static_cast<long>(y_.size())) {
      cf cur = interp(pos_);
      cf mid = interp(pos_ - 1.0);
      float et = (cur.r - prev_.r) * mid.r + (cur.i - prev_.i) * mid.i;
      et = std::max(-2.f, std::min(2.f, et));
      prev_ = cur;
      tint_ += kTi * et;
      pos_ += 2.0 + 0.0 - (kTp * et + tint_);

      // carrier derotation z = cur * conj(ph)
      cf z = {cur.r * ph_.r + cur.i * ph_.i, cur.i * ph_.r - cur.r * ph_.i};
      float sr = z.r >= 0 ? 1.f : -1.f, si = z.i >= 0 ? 1.f : -1.f;
      float mag = std::sqrt(z.r * z.r + z.i * z.i) + 1e-6f;
      float ep = (sr * z.i - si * z.r) / mag;         // decision-directed QPSK phase error
      freq_ += kCb * ep;
      float d = freq_ + kCa * ep;                     // phase increment (rad), small
      float c = 1 - d * d / 2 + d * d * d * d / 24, s = d - d * d * d / 6;
      cf np = {ph_.r * c - ph_.i * s, ph_.r * s + ph_.i * c};
      float m2 = np.r * np.r + np.i * np.i;
      float nrm = 1.5f - 0.5f * m2;                   // cheap re-normalisation
      ph_ = {np.r * nrm, np.i * nrm};

      // Hard decision.  Modulator: I = +-1 from B'' (bit 0 -> +1),
      // Q = -(+-1) from A''.  => B'' = (I<0), A'' = (Q>0).
      u8 a = z.i > 0, b = z.r < 0;
      out.push_back(static_cast<u8>((a << 1) | b));
    }
  }

  bool real_;
  int D_;
  FirFilter mf_, dec_;
  std::vector<cf> tmp_, s2_, mfo_, y_;
  int mix_ = 0;
  float pwr_ = 1.f;
  long base_ = 0;
  double pos_ = 3.0;
  cf prev_{0, 0};
  float tint_ = 0;
  cf ph_{1, 0};
  float freq_ = 0;
};

}  // namespace dsr
