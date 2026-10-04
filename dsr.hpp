// ============================================================================
//  dsr.hpp  -  Core definitions of the DSR (Digital Satellite Radio) system
// ============================================================================
//
//  Reference: "Technische Richtlinie 3 R 1, Digitaler Satelliten Rundfunk
//  (DSR), Spezifikation des Uebertragungsverfahrens", IRT, Ausgabe 3, Nov. 1989.
//  Section numbers in the comments below ("spec 3.2.2" ...) refer to it.
//
//  This header contains everything that is *bit-exact protocol logic* and is
//  shared by the encoder and the decoder:
//
//    1. system constants
//    2. generic systematic BCH encoder / syndrome decoder
//    3. audio coding: 16 bit -> 14 bit floating point with a 3 bit scale factor
//    4. 77-bit block, double-block interleaving, main frame (320 bit)
//    5. scrambler and differential encoder / decoder
//    6. special service channel (programme offer PA, station name SK)
//    7. programme-associated information (PI) packets
//
//  Conventions
//  -----------
//  * A "bit" is stored in one byte (value 0 or 1).  This wastes memory but keeps
//    the code close to the figures of the specification.
//  * In every bit/byte sequence the first element is the first one that is
//    transmitted (spec chapter 1: "Leseweise von links nach rechts").
//  * A "stereo channel" (station) is a pair of audio channels (L, R).  The
//    system carries 16 stereo channels: 1..8 in main frame A, 9..16 in frame B.
//    In the code the stations are numbered 0..15.
// ============================================================================
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

namespace dsr {

using u8  = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

// ----------------------------------------------------------------------------
// 1. System constants
// ----------------------------------------------------------------------------
constexpr int    kAudioRate      = 32000;   // audio sampling rate = main frame rate (Hz)
constexpr int    kBlockSamples   = 64;      // samples per scale-factor block = 2 ms
constexpr int    kStations       = 16;      // stereo programmes per transponder
constexpr int    kFrameBits      = 320;     // bits per main frame (spec 3.2)
constexpr int    kSyncBits       = 11;      // Barker word length
constexpr int    kSbitIndex      = 11;      // position of the special-service bit
constexpr int    kFirstScrambled = 12;      // scrambling starts at the 13th bit
constexpr int    kScrambledBits  = 308;     // 320 - 12
constexpr int    kBlockBits      = 77;      // one 77-bit block (spec 3.2.2)
constexpr double kSymbolRate     = 10.24e6; // QPSK symbols per second (= bits/s per frame stream)
constexpr double kRollOff        = 0.5;     // cosine roll-off factor (spec 4.4)

// Barker-11 synchronisation words.  Frame A carries the word, frame B the
// inverted word (spec 3.2.1).
constexpr u8 kBarkerA[kSyncBits] = {1,1,1,0,0,0,1,0,0,1,0};
constexpr u8 kBarkerB[kSyncBits] = {0,0,0,1,1,1,0,1,1,0,1};

// 16-bit Williard words carried in the special-service bits of frame A
// (spec 3.3.2 and 3.4.3).  SYNC1 marks the start of a special-service
// super-frame (SAUE), SYNC2 starts the 7 following special-service frames.
constexpr u16 kSync1 = 0x05CF;  // 0000 0101 1100 1111
constexpr u16 kSync2 = 0x05FF;  // 0000 0101 1111 1111

// ----------------------------------------------------------------------------
// 2. Systematic BCH code with table-driven error correction
// ----------------------------------------------------------------------------
// Textbook background:  a cyclic code of length n is defined by a generator
// polynomial g(x) of degree r.  Systematic encoding appends the r parity bits
//   p(x) = d(x) * x^r  mod g(x)
// to the k data bits d.  The received word is a codeword iff its remainder
// modulo g(x) (the "syndrome") is zero.  For small codes it is cheap to
// pre-compute the syndromes of ALL error patterns up to weight t and look the
// received syndrome up in a hash table: that is a complete bounded-distance
// decoder without any Galois-field arithmetic.
//
// The remainder is computed with the classical bit-serial LFSR division
// circuit; bits are fed most significant (= first transmitted) bit first.
class Bch {
 public:
  // gen  : generator polynomial including the x^deg term, bit i = coefficient of x^i
  // deg  : degree of the generator = number of parity bits
  // n    : code word length (shortened codes: length actually transmitted)
  // tmax : correct all error patterns up to this weight
  Bch(u32 gen, int deg, int n, int tmax)
      : deg_(deg), n_(n), mask_((1u << deg) - 1), low_(gen & ((1u << deg) - 1)) {
    build_table(0, 0, 0, tmax);
  }

  // Remainder of the (zero-extended) bit sequence, see class comment.
  u32 syndrome(const u8* bits, int count) const {
    u32 r = 0;
    for (int i = 0; i < count; ++i) {
      u32 fb = ((r >> (deg_ - 1)) & 1u) ^ bits[i];
      r = (r << 1) & mask_;
      if (fb) r ^= low_;
    }
    return r;
  }

  // Compute `deg` parity bits for `k` data bits (highest power first).
  void parity(const u8* data, int k, u8* par) const {
    u32 r = syndrome(data, k);
    for (int i = 0; i < deg_; ++i) par[i] = (r >> (deg_ - 1 - i)) & 1u;
  }

  // Correct `w` (n bits) in place.  Returns the number of corrected bits or
  // -1 if the error pattern is not correctable.
  int correct(u8* w) const {
    u32 s = syndrome(w, n_);
    if (s == 0) return 0;
    auto it = table_.find(s);
    if (it == table_.end()) return -1;
    int c = 0;
    for (int i = 0; i < n_; ++i)
      if ((it->second >> i) & 1u) { w[i] ^= 1; ++c; }
    return c;
  }

 private:
  // Enumerate all error patterns of weight 1..tmax (positions in increasing order).
  void build_table(int start, int depth, u64 pat, int tmax) {
    for (int p = start; p < n_; ++p) {
      u64 np = pat | (1ull << p);
      u8 bits[64] = {0};
      for (int i = 0; i < n_; ++i) bits[i] = (np >> i) & 1u;
      table_.emplace(syndrome(bits, n_), np);
      if (depth + 1 < tmax) build_table(p + 1, depth + 1, np, tmax);
    }
  }
  int deg_, n_;
  u32 mask_, low_;
  std::unordered_map<u32, u64> table_;
};

// (63,44) BCH code protecting the 11 MSBs of four audio channels (spec 3.2.2):
//   g(x) = x^19 + x^15 + x^10 + x^9 + x^8 + x^6 + x^4 + 1
// Minimum distance 8 -> corrects 3 errors and reliably detects 4.
inline const Bch& bch63_44() {
  static const Bch code((1u << 19) | (1u << 15) | (1u << 10) | (1u << 9) |
                            (1u << 8) | (1u << 6) | (1u << 4) | 1u,
                        19, 63, 3);
  return code;
}

// Shortened (14,6) BCH code derived from BCH(15,7), protecting the two scale
// factors of a stereo pair (spec 3.5.2):  g(x) = x^8 + x^7 + x^6 + x^4 + 1.
inline const Bch& bch15_7() {
  static const Bch code((1u << 8) | (1u << 7) | (1u << 6) | (1u << 4) | 1u, 8, 15, 2);
  return code;
}

// ----------------------------------------------------------------------------
// 3. Audio coding: 16 bit PCM <-> 14 bit "floating point" (spec 2.2 / 2.3)
// ----------------------------------------------------------------------------
// A block of 64 samples shares one 3-bit scale factor k (0..7).  k is the
// number of bits directly after the sign bit Y1 that are identical to the sign
// bit in ALL 64 samples of the block.  Those k redundant bits are removed and
// the remaining bits are shifted towards the sign bit (Bild 1).  The 14-bit
// code word is therefore
//
//     Y1, Y(k+2), Y(k+3), ..., Y16, then zero padding up to 14 bits
//
// (Y1 = sign / MSB, Y16 = LSB).  With k = 0 the two LSBs Y15/Y16 are lost
// (14 bit resolution), with large k all 16 bits survive for quiet signals.

// Number of redundant sign bits of one sample, limited to 7.
inline int sample_redundancy(std::int16_t v) {
  u16 x = static_cast<u16>(v ^ (v >> 15));      // sign cancelled: leading zeros = redundancy+1
  if (x == 0) return 7;
  int lz = __builtin_clz(static_cast<unsigned>(x)) - 16;   // leading zeros inside 16 bits
  return std::min(7, lz - 1);
}

// Scale factor of a whole block = minimum redundancy of all samples.
inline int block_scale_factor(const std::int16_t* s, int n = kBlockSamples) {
  int k = 7;
  for (int i = 0; i < n; ++i) k = std::min(k, sample_redundancy(s[i]));
  return k;
}

// 16-bit sample -> 14-bit code word (first transmitted bit = bit 13).
inline u16 encode_sample(std::int16_t v, int k) {
  u16 u = static_cast<u16>(v);
  auto Y = [&](int i) { return (u >> (16 - i)) & 1u; };   // Y(i), i = 1..16
  u16 w = static_cast<u16>(Y(1) << 13);
  for (int i = 1; i <= 13; ++i) {
    int idx = k + 1 + i;
    if (idx <= 16) w |= static_cast<u16>(Y(idx) << (13 - i));
  }
  return w;
}

// 14-bit code word -> 16-bit sample.  Missing LSBs are filled with zeros.
inline std::int16_t decode_sample(u16 w, int k) {
  auto T = [&](int i) { return (w >> (13 - i)) & 1u; };   // transmitted bit i = 0..13
  u16 u = 0;
  u32 sign = T(0);
  u |= static_cast<u16>(sign << 15);
  for (int i = 2; i <= k + 1; ++i) u |= static_cast<u16>(sign << (16 - i));   // restore redundancy
  for (int j = 0; j <= 12 && k + 2 + j <= 16; ++j)
    u |= static_cast<u16>(T(1 + j) << (16 - (k + 2 + j)));
  return static_cast<std::int16_t>(u);
}

// ----------------------------------------------------------------------------
// 3b. Scale factor protection inside the ZI frame (spec 3.5.2)
// ----------------------------------------------------------------------------
// Data: a leading 0, then SkL (3 bit, MSB first), then SkR (3 bit).  8 BCH
// parity bits are appended; the leading 0 is not transmitted -> 14 bits.
inline void zi_encode(int skL, int skR, u8 out14[14]) {
  u8 d[7] = {0, u8((skL >> 2) & 1), u8((skL >> 1) & 1), u8(skL & 1),
             u8((skR >> 2) & 1), u8((skR >> 1) & 1), u8(skR & 1)};
  u8 par[8];
  bch15_7().parity(d, 7, par);
  for (int i = 0; i < 6; ++i) out14[i] = d[i + 1];
  for (int i = 0; i < 8; ++i) out14[6 + i] = par[i];
}

// Returns false if the 14 bits are not correctable.
inline bool zi_decode(const u8 in14[14], int& skL, int& skR) {
  u8 w[15];
  w[0] = 0;
  std::memcpy(w + 1, in14, 14);
  if (bch15_7().correct(w) < 0 || w[0] != 0) return false;
  skL = (w[1] << 2) | (w[2] << 1) | w[3];
  skR = (w[4] << 2) | (w[5] << 1) | w[6];
  return true;
}

// ----------------------------------------------------------------------------
// 4. The 77-bit block, double blocks and the 320-bit main frame
// ----------------------------------------------------------------------------
// One 77-bit block serves TWO stereo channels (I, II), see Bild 3:
//
//   | MSB(LI) 11 | MSB(RI) 11 | MSB(LII) 11 | MSB(RII) 11 | BCH 19 | ZI(I) ZI(II) | LSB 4x3 |
//   '------------------- BCH(63,44) code word ---------------------'
//
// cw[0..3] = 14-bit code words of LI, RI, LII, RII; zi[0..1] = ZI bits of I, II.
inline void build_block(const u16 cw[4], const u8 zi[2], u8 out[kBlockBits]) {
  u8 msb[44];
  for (int c = 0; c < 4; ++c)
    for (int b = 0; b < 11; ++b) msb[c * 11 + b] = (cw[c] >> (13 - b)) & 1u;
  u8 par[19];
  bch63_44().parity(msb, 44, par);
  std::memcpy(out, msb, 44);
  std::memcpy(out + 44, par, 19);
  out[63] = zi[0];
  out[64] = zi[1];
  for (int c = 0; c < 4; ++c)
    for (int b = 0; b < 3; ++b) out[65 + c * 3 + b] = (cw[c] >> (2 - b)) & 1u;
}

// Inverse operation.  Returns the number of corrected bits, or -1 if the BCH
// word was uncorrectable (the caller then conceals the audio).
inline int parse_block(const u8 in[kBlockBits], u16 cw[4], u8 zi[2]) {
  u8 w[63];
  std::memcpy(w, in, 63);
  int e = bch63_44().correct(w);
  if (e < 0) std::memcpy(w, in, 63);          // keep the raw bits for the caller
  for (int c = 0; c < 4; ++c) {
    u16 v = 0;
    for (int b = 0; b < 11; ++b) v = static_cast<u16>((v << 1) | w[c * 11 + b]);
    for (int b = 0; b < 3; ++b) v = static_cast<u16>((v << 1) | in[65 + c * 3 + b]);
    cw[c] = v;
  }
  zi[0] = in[63];
  zi[1] = in[64];
  return e;
}

// Bit-wise interleaving ("verkaemmen") of two consecutive 77-bit blocks into a
// 154-bit double block: x0 y0 x1 y1 ... (spec 3.2).  Purpose: a double bit
// error produced by a differential demodulator hits two *different* blocks.
inline void interleave(const u8 x[kBlockBits], const u8 y[kBlockBits], u8 out[2 * kBlockBits]) {
  for (int i = 0; i < kBlockBits; ++i) { out[2 * i] = x[i]; out[2 * i + 1] = y[i]; }
}
inline void deinterleave(const u8 in[2 * kBlockBits], u8 x[kBlockBits], u8 y[kBlockBits]) {
  for (int i = 0; i < kBlockBits; ++i) { x[i] = in[2 * i]; y[i] = in[2 * i + 1]; }
}

// Everything a main frame carries for one stereo channel.
struct StationSlot {
  u16 cw[2];   // 14-bit code words of the left / right channel
  u8  zi;      // one bit of the additional-information (ZI) frame
};

// Main frame layout (Bild 2 / Bild 6):
//   [ Barker 11 | S 1 | double block 1 (154) | double block 2 (154) ] = 320 bit
//   double block 1 = blocks of stereo channels (1,2) and (3,4)   (A)  or (9,10),(11,12) (B)
//   double block 2 = blocks of stereo channels (5,6) and (7,8)   (A)  or (13,14),(15,16) (B)
inline void build_frame(bool is_b, const StationSlot slot[8], u8 sbit, u8 out[kFrameBits]) {
  std::memcpy(out, is_b ? kBarkerB : kBarkerA, kSyncBits);
  out[kSbitIndex] = sbit;
  u8 blk[4][kBlockBits];
  for (int b = 0; b < 4; ++b) {
    u16 cw[4] = {slot[2 * b].cw[0], slot[2 * b].cw[1], slot[2 * b + 1].cw[0], slot[2 * b + 1].cw[1]};
    u8 zi[2] = {slot[2 * b].zi, slot[2 * b + 1].zi};
    build_block(cw, zi, blk[b]);
  }
  interleave(blk[0], blk[1], out + 12);
  interleave(blk[2], blk[3], out + 12 + 2 * kBlockBits);
}

// Parse a (descrambled) main frame.  bad[i] is set for stereo channel i when the
// BCH word of its block was uncorrectable.  Returns the number of corrected bits.
inline int parse_frame(const u8 in[kFrameBits], StationSlot slot[8], bool bad[8]) {
  u8 blk[4][kBlockBits];
  deinterleave(in + 12, blk[0], blk[1]);
  deinterleave(in + 12 + 2 * kBlockBits, blk[2], blk[3]);
  int corrected = 0;
  for (int b = 0; b < 4; ++b) {
    u16 cw[4];
    u8 zi[2];
    int e = parse_block(blk[b], cw, zi);
    slot[2 * b].cw[0] = cw[0];     slot[2 * b].cw[1] = cw[1];
    slot[2 * b + 1].cw[0] = cw[2]; slot[2 * b + 1].cw[1] = cw[3];
    slot[2 * b].zi = zi[0];        slot[2 * b + 1].zi = zi[1];
    bad[2 * b] = bad[2 * b + 1] = (e < 0);
    if (e > 0) corrected += e;
  }
  return corrected;
}

// ----------------------------------------------------------------------------
// 5. Scrambler and differential encoder (spec 4.2, 4.3)
// ----------------------------------------------------------------------------
// 9-stage LFSR, g(x) = x^9 + x^4 + 1, clocked at the bit rate.  It is reset to
// 0 1 0 1 1 1 1 0 1 (r8 ... r0) at the 13th bit of EVERY main frame, so the
// 308 mask bits are identical for every frame and can be pre-computed.
//   frame A is XORed with r0,   frame B with r3 XOR r0.
// The Barker word and the S bit (first 12 bits) are NOT scrambled.
struct ScramblerMasks {
  u8 a[kScrambledBits];
  u8 b[kScrambledBits];
  ScramblerMasks() {
    u8 r[9] = {1, 0, 1, 1, 1, 1, 0, 1, 0};   // r[0] = r0 ... r[8] = r8  (init 0 1 0 1 1 1 1 0 1 read r8..r0)
    for (int i = 0; i < kScrambledBits; ++i) {
      a[i] = r[0];
      b[i] = r[3] ^ r[0];
      u8 fb = r[0] ^ r[4];                    // taps x^0 and x^4 feed back into r8 (x^9)
      for (int k = 0; k < 8; ++k) r[k] = r[k + 1];
      r[8] = fb;
    }
  }
};
inline const ScramblerMasks& scrambler_masks() {
  static const ScramblerMasks m;
  return m;
}
// Scrambling and descrambling are the same XOR operation.
inline void scramble_frame(bool is_b, u8 f[kFrameBits]) {
  const u8* m = is_b ? scrambler_masks().b : scrambler_masks().a;
  for (int i = 0; i < kScrambledBits; ++i) f[kFirstScrambled + i] ^= m[i];
}

// Differential encoder (spec 4.3).  (a,b) = previous output pair (A''_{n-1}, B''_{n-1}).
//   if A'^B' == 0 :  A'' = a ^ A'   B'' = b ^ B'
//   if A'^B' == 1 :  A'' = b ^ A'   B'' = a ^ B'      (quadrant swap = +-90 deg)
// The mapping makes the receiver immune to the 4-fold QPSK carrier ambiguity.
struct DiffEncoder {
  u8 a = 0, b = 0;
  void step(u8 A1, u8 B1, u8& A2, u8& B2) {
    if ((A1 ^ B1) == 0) { A2 = a ^ A1; B2 = b ^ B1; }
    else                { A2 = b ^ A1; B2 = a ^ B1; }
    a = A2; b = B2;
  }
};
// Inverse.  S = A'^B' is recovered from the parity change of the pair.
struct DiffDecoder {
  u8 a = 0, b = 0;
  void step(u8 A2, u8 B2, u8& A1, u8& B1) {
    u8 s = (A2 ^ B2) ^ (a ^ b);
    if (s == 0) { A1 = A2 ^ a; B1 = B2 ^ b; }
    else        { A1 = A2 ^ b; B1 = B2 ^ a; }
    a = A2; b = B2;
  }
};

// ----------------------------------------------------------------------------
// 6. Special-service channel: programme offer (PA) and station name (SK)
// ----------------------------------------------------------------------------
// The S bits of 64 consecutive frames A form a 64-bit special-service frame SA:
//     16 bit sync | 6 bytes
// Bytes 1..4 hold PA or SK data of four mono channels (= 2 stereo channels),
// bytes 5/6 are the D/E bytes (spec 3.4.2).  8 SA frames (16 ms) form a
// special-service super-frame SAUE covering all 16 stereo channels, and 16 of
// those (256 ms) form SAUEUE: 7 x SAUE/PA, 1 x SAUE/LB (empty), 8 x SAUE/SK
// (one station-name character per SAUE/SK) (spec 3.4.3, 3.4.4, Bild 4).

// Programme types (spec 3.4.2.1 table, identical to RDS PTY of EBU Tech 3244).
inline const char* const kPtyNames[16] = {
    "KEINE",   "NACHRICH", "POLITIK", "SPEZWORT", "SPORT",   "LERNEN",
    "HOER+LIT","KULTUR",   "WISSEN",  "UNTERHAL", "POP",     "ROCK",
    "U-MUSIK", "L-KLASS",  "E-KLASS", "SPEZ MUS"};

struct StationInfo {
  bool active = false;       // is the channel used at all?
  u8   name[8] = {' ',' ',' ',' ',' ',' ',' ',' '};   // SK characters (RDS code)
  int  pty = 0;              // programme type 0..15 (main)
  int  sub = -1;             // programme sub-type 0..15 (-1: repeat main type)
  bool music = true;         // speech / music flag K
};

// Set the parity bit (bit 0) so that the whole byte has even parity (spec: P=0
// for an even number of ones in bits 1..7).
inline u8 with_parity(u8 b) {
  b &= 0xFE;
  return static_cast<u8>(b | (__builtin_parity(b) & 1));
}
constexpr u8 kPaUnused = 0x09;   // 0000 1001: channel not in use (spec 3.4.2.1)

// PA-L : pppp K 0 1 P      PA-R : ssss 0 1 0 P     (stereo pair; bits 6/7 = 01 and 10)
inline u8 pa_left(const StationInfo& s)  { return with_parity(u8((s.pty << 4) | (s.music << 3) | 0x02)); }
inline u8 pa_right(const StationInfo& s) {
  int sub = s.sub >= 0 ? s.sub : s.pty;    // "repeat the main type if there is no sub-type"
  return with_parity(u8((sub << 4) | 0x04));
}

// Map a UTF-8 string to 8 SK characters.  RDS/EBU Tech 3244 character set; code
// words 0x7F and 0xE0..0xFF are excluded by the spec (3.4.2.2) because they could
// imitate the sync words.  ASCII 0x20..0x7E is identical in RDS; a few German
// letters are mapped explicitly (see table 1 of the spec).  Others become '?'.
inline u8 sk_map_codepoint(u32 cp) {
  // Where the RDS table differs from ASCII: '$' sits at 0xAB (0x24 is the currency
  // sign), and ^ ` ~ do not exist at all.
  if (cp == '$') return 0xAB;
  if (cp == '^' || cp == '`' || cp == '~') return '?';
  if (cp >= 0x20 && cp <= 0x7E) return static_cast<u8>(cp);
  switch (cp) {
    case 0xE4: return 0x91;  // ae
    case 0xF6: return 0x97;  // oe
    case 0xFC: return 0x99;  // ue
    case 0xDF: return 0x8D;  // sz
    case 0xC4: return 0xD1;  // AE
    case 0xD6: return 0xD7;  // OE
    case 0xDC: return 0xD9;  // UE
    default:   return '?';
  }
}
inline void sk_from_utf8(const std::string& s, u8 out[8]) {
  std::memset(out, ' ', 8);
  int n = 0;
  for (size_t i = 0; i < s.size() && n < 8;) {
    u8 c = static_cast<u8>(s[i]);
    u32 cp; int len;
    if (c < 0x80)            { cp = c;        len = 1; }
    else if ((c >> 5) == 6)  { cp = c & 0x1F; len = 2; }
    else if ((c >> 4) == 14) { cp = c & 0x0F; len = 3; }
    else                     { cp = c & 0x07; len = 4; }
    for (int k = 1; k < len && i + k < s.size(); ++k) cp = (cp << 6) | (static_cast<u8>(s[i + k]) & 0x3F);
    i += len;
    out[n++] = sk_map_codepoint(cp);
  }
}
inline std::string sk_to_utf8(u8 c) {
  switch (c) {
    case 0x91: return "\xC3\xA4"; case 0x97: return "\xC3\xB6"; case 0x99: return "\xC3\xBC";
    case 0x8D: return "\xC3\x9F"; case 0xD1: return "\xC3\x84"; case 0xD7: return "\xC3\x96";
    case 0xD9: return "\xC3\x9C"; case 0xAB: return "$";
    default: return (c >= 0x20 && c <= 0x7E) ? std::string(1, char(c)) : std::string("?");
  }
}

// Content (6 bytes) of special-service frame number `sa_no` (counted from the
// stream start, one per 2 ms).  Position inside the 256 ms SAUEUE:
//     sau = (sa_no / 8) % 16 :  0..6 PA, 7 LB (all zero), 8..15 SK character 0..7
//     k   =  sa_no % 8       :  stereo channels 2k, 2k+1
inline void make_sa_bytes(const StationInfo st[kStations], u64 sa_no, u8 out[6]) {
  int sau = static_cast<int>((sa_no / 8) % 16);
  int k   = static_cast<int>(sa_no % 8);
  std::memset(out, 0, 6);                        // D/E bytes and the LB frame stay zero
  if (sau == 7) return;
  const StationInfo& s0 = st[2 * k];
  const StationInfo& s1 = st[2 * k + 1];
  if (sau < 7) {                                 // programme offer
    out[0] = s0.active ? pa_left(s0)  : kPaUnused;
    out[1] = s0.active ? pa_right(s0) : kPaUnused;
    out[2] = s1.active ? pa_left(s1)  : kPaUnused;
    out[3] = s1.active ? pa_right(s1) : kPaUnused;
  } else {                                       // station name: same character in L and R
    int c = sau - 8;
    out[0] = out[1] = s0.name[c];
    out[2] = out[3] = s1.name[c];
  }
}

// ----------------------------------------------------------------------------
// 7. Programme-associated information (PI), Anhang 2
// ----------------------------------------------------------------------------
// Each ZI frame (2 ms) carries one 22-bit PI word per stereo channel.  PI data is
// organised in packets: 2 header words (44 bit) + n content words.  Header:
//   12 bit start word 000000111111 | length (2 bytes) | content id (2 bytes)
// where every byte is two Hamming(8,4) coded nibbles (Tabelle 2).

// Hamming(8,4): information bits b8,b6,b4,b2, protection bits b7,b5,b3,b1.
inline u8 hamming84_encode(u8 nib) {
  u8 b8 = (nib >> 3) & 1, b6 = (nib >> 2) & 1, b4 = (nib >> 1) & 1, b2 = nib & 1;
  u8 b7 = b8 ^ b6 ^ b4;
  u8 b5 = b6 ^ b4 ^ (b2 ^ 1);
  u8 b3 = (b8 ^ 1) ^ b4 ^ b2;
  u8 b1 = (b2 ^ 1) ^ b8 ^ b6;
  return static_cast<u8>((b8 << 7) | (b7 << 6) | (b6 << 5) | (b5 << 4) | (b4 << 3) | (b3 << 2) | (b2 << 1) | b1);
}
// Returns the nibble (0..15) or -1 if not correctable.  Distance-4 code:
// 1 error is corrected, 2 are detected.
inline int hamming84_decode(u8 byte) {
  for (int n = 0; n < 16; ++n)
    if (__builtin_popcount(hamming84_encode(static_cast<u8>(n)) ^ byte) <= 1) return n;
  return -1;
}

// Emits an endless sequence of "dummy" PI packets (length 0, content id 0): the
// minimal valid use of the PI channel.  A real service would insert its payload here.
class PiGenerator {
 public:
  PiGenerator() {
    int n = 0;
    for (int i = 11; i >= 0; --i) bits_[n++] = (0x03F >> i) & 1;           // 000000111111
    for (int rep = 0; rep < 2; ++rep) {                                      // length=0, id=0
      u8 hi = hamming84_encode(0), lo = hamming84_encode(0);
      for (int i = 7; i >= 0; --i) bits_[n++] = (hi >> i) & 1;
      for (int i = 7; i >= 0; --i) bits_[n++] = (lo >> i) & 1;
    }
  }
  // Next 22-bit word (MSB = first transmitted bit).
  u32 next_word() {
    u32 w = 0;
    for (int i = 0; i < 22; ++i) w = (w << 1) | bits_[pos_ + i];
    pos_ = (pos_ + 22) % 44;
    return w;
  }
 private:
  u8 bits_[44];
  int pos_ = 0;
};

}  // namespace dsr
