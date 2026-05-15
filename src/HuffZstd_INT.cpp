#include "HuffZstd_INT.h"
#include "sperr_helper.h"

// Vendored SZ3 lossless pipeline, with SZo's speed-optimized Huffman encoder.
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "zstd.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Payload format version. byte 9 of m_payload carries this.
//   1 = SZo Huffman + ZigZag bins, raster scan order (Sprint 4 layout)
//   2 = adds subband-order permutation before Huffman (Sprint 5)
// The permutation is fully determined by chunk dims (via can_use_dyadic and
// calc_approx_detail_len), so no metadata travels with the payload.
constexpr uint8_t kFormatVersion = 2;

// Axis-aligned box [xa, xb) × [ya, yb) × [za, zb) within the raster-ordered
// coefficient buffer. Sprint 5 uses these to permute bins by subband.
struct SubbandBox {
  size_t xa, xb, ya, yb, za, zb;
};

// Enumerate subbands of SPERR's 3D wavelet transform, ordered lowest-freq
// first. Handles BOTH transform paths in CDF97::dwt3d:
//   - Dyadic Mallat (can_use_dyadic returns L): 1 + 7L subbands. LLL at the
//     deepest level (smallest box, lowest 3D freq) first, then 7 detail
//     subbands at each lev = L down to 1.
//   - Wavelet packet (non-dyadic dims, e.g. Uf48 500×500×100 chunked to
//     256×256×100 — z has fewer levels than xy): Z-axis is 1D-Mallat'd with
//     Lz levels independently, then each Z slice is 2D-Mallat'd with Lxy
//     levels. Total (Lz+1)·(1+3·Lxy) subbands. We enumerate Z-bands lowest
//     freq first, and within each Z-band the 2D-LL first then 2D detail
//     subbands from deepest level outward.
auto enumerate_subbands(sperr::dims_type dims) -> std::vector<SubbandBox>
{
  std::vector<SubbandBox> out;
  auto dy = sperr::can_use_dyadic(dims);

  if (dy) {
    // 3D dyadic path.
    const size_t L = *dy;
    std::vector<size_t> ax(L + 1), ay(L + 1), az(L + 1);
    for (size_t k = 0; k <= L; ++k) {
      ax[k] = sperr::calc_approx_detail_len(dims[0], k)[0];
      ay[k] = sperr::calc_approx_detail_len(dims[1], k)[0];
      az[k] = sperr::calc_approx_detail_len(dims[2], k)[0];
    }
    out.reserve(1 + 7 * L);
    out.push_back({0, ax[L], 0, ay[L], 0, az[L]});  // LLL
    for (size_t lev = L; lev > 0; --lev) {
      const size_t xS = ax[lev], xL = ax[lev - 1];
      const size_t yS = ay[lev], yL = ay[lev - 1];
      const size_t zS = az[lev], zL = az[lev - 1];
      out.push_back({xS, xL,  0, yS,  0, zS});  // HLL
      out.push_back({ 0, xS, yS, yL,  0, zS});  // LHL
      out.push_back({xS, xL, yS, yL,  0, zS});  // HHL
      out.push_back({ 0, xS,  0, yS, zS, zL});  // LLH
      out.push_back({xS, xL,  0, yS, zS, zL});  // HLH
      out.push_back({ 0, xS, yS, yL, zS, zL});  // LHH
      out.push_back({xS, xL, yS, yL, zS, zL});  // HHH
    }
    return out;
  }

  // Wavelet-packet path: Z-1D then XY-2D, independent level counts.
  const size_t Lz  = sperr::num_of_xforms(dims[2]);
  const size_t Lxy = sperr::num_of_xforms(std::min(dims[0], dims[1]));
  if (Lz == 0 && Lxy == 0) {
    // Trivially small dims — no decomposition happens at all.
    out.push_back({0, dims[0], 0, dims[1], 0, dims[2]});
    return out;
  }

  std::vector<size_t> ax(Lxy + 1), ay(Lxy + 1), az(Lz + 1);
  for (size_t k = 0; k <= Lxy; ++k) {
    ax[k] = sperr::calc_approx_detail_len(dims[0], k)[0];
    ay[k] = sperr::calc_approx_detail_len(dims[1], k)[0];
  }
  for (size_t k = 0; k <= Lz; ++k)
    az[k] = sperr::calc_approx_detail_len(dims[2], k)[0];

  // Z-band ranges, lowest freq first: [0, az[Lz]) then [az[lev], az[lev-1]) for
  // lev = Lz..1.
  std::vector<std::pair<size_t, size_t>> z_bands;
  z_bands.reserve(Lz + 1);
  z_bands.emplace_back(0, az[Lz]);
  for (size_t lev = Lz; lev > 0; --lev)
    z_bands.emplace_back(az[lev], az[lev - 1]);

  out.reserve((Lz + 1) * (1 + 3 * Lxy));
  for (const auto& [za, zb] : z_bands) {
    // 2D-LL at the deepest XY level (lowest 2D freq within this Z-band).
    out.push_back({0, ax[Lxy], 0, ay[Lxy], za, zb});
    // 3 detail subbands per 2D level, deepest first.
    for (size_t lev = Lxy; lev > 0; --lev) {
      const size_t xS = ax[lev], xL = ax[lev - 1];
      const size_t yS = ay[lev], yL = ay[lev - 1];
      out.push_back({xS, xL,  0, yS, za, zb});  // HL
      out.push_back({ 0, xS, yS, yL, za, zb});  // LH
      out.push_back({xS, xL, yS, yL, za, zb});  // HH
    }
  }
  return out;
}

// ZigZag map signed int32 to non-negative int32 so the Huffman alphabet starts
// at 0. Most-common signed value (0, dominant in SPERR's quantized output) maps
// to bin=0, which is exactly what SZo's `valid_len1` decoder fast path expects.
inline int32_t zigzag_encode(int32_t v) noexcept {
  return (v << 1) ^ (v >> 31);
}
inline int32_t zigzag_decode(int32_t b) noexcept {
  return (static_cast<uint32_t>(b) >> 1) ^ -(b & 1);
}

// Encoder workhorse. Inputs are already non-negative ZigZag bins; we know the
// alphabet bounds (stateNum) and pass a pre-built dense freq array so SZo's
// 65 K-cap init scan is skipped entirely.
auto encode_huffzstd(const std::vector<int32_t>& bins,
                     const std::vector<size_t>& freq,
                     int stateNum) -> sperr::vec8_type
{
  const size_t N = bins.size();

  // Scratch budget: tree header is O(stateNum) bytes (SZo writes 1 byte per
  // slot); bitstream is at most ~8 bytes per symbol (very loose upper bound,
  // never approached in practice).
  std::vector<SZ3::uchar> huff_scratch(8 * N + size_t(stateNum) + 4096);

  const auto t0 = std::chrono::steady_clock::now();
  SZ3::HuffmanEncoder<int32_t> enc;
  enc.preprocess_encode_dense(N, const_cast<size_t*>(freq.data()), stateNum);
  SZ3::uchar* hp = huff_scratch.data();
  enc.save(hp);
  const size_t tree_len = static_cast<size_t>(hp - huff_scratch.data());
  enc.encode(bins.data(), N, hp);
  enc.postprocess_encode();
  const size_t huff_len = static_cast<size_t>(hp - huff_scratch.data());
  const size_t stream_len = huff_len - tree_len;
  const auto t1 = std::chrono::steady_clock::now();

  std::vector<SZ3::uchar> zst_buf(ZSTD_compressBound(huff_len) + sizeof(size_t));
  SZ3::Lossless_zstd zstd;
  const size_t zst_len = zstd.compress(huff_scratch.data(),
                                       huff_len,
                                       zst_buf.data(),
                                       zst_buf.size());
  const auto t2 = std::chrono::steady_clock::now();

  if (const char* env = std::getenv("SPERR_HZSTATS"); env && env[0] != '0') {
    const double t_huff =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    const double t_zst =
        std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
    std::fprintf(stderr,
                 "[HZSTAT] N=%zu stateNum=%d tree=%zu stream=%zu huff_total=%zu "
                 "zst_out=%zu zst_ratio=%.4f t_huff=%.4fs t_zst=%.4fs\n",
                 N, stateNum, tree_len, stream_len, huff_len, zst_len,
                 double(zst_len) / double(huff_len), t_huff, t_zst);
  }

  sperr::vec8_type out(zst_len);
  std::memcpy(out.data(), zst_buf.data(), zst_len);
  return out;
}

auto decode_huffzstd(const uint8_t* zst_src, size_t zst_src_len, size_t N)
    -> std::unique_ptr<int32_t[]>
{
  SZ3::uchar* huff_buf = nullptr;
  size_t huff_len = 0;
  SZ3::Lossless_zstd zstd;
  zstd.decompress(const_cast<SZ3::uchar*>(zst_src), zst_src_len, huff_buf, huff_len);

  SZ3::HuffmanEncoder<int32_t> dec;
  const SZ3::uchar* rp = huff_buf;
  size_t remaining = huff_len;
  dec.load(rp, remaining);
  // SZo's decode returns a raw new[] pointer (with SIMD tail padding).
  std::unique_ptr<int32_t[]> bins{dec.decode(rp, N)};
  dec.postprocess_decode();
  std::free(huff_buf);
  return bins;
}

}  // namespace

template <typename T>
void sperr::HuffZstd_INT<T>::encode()
{
  const auto& mags = this->m_coeff_buf;
  auto& signs = this->m_sign_array;
  const size_t N = mags.size();

  // Largest magnitude → num_bitplanes_eq (kept identical to SPECK_INT so
  // SPECK_FLT's uint-width selection still works).
  uint64_t max_mag = 0;
  for (size_t i = 0; i < N; ++i) {
    const uint64_t v = static_cast<uint64_t>(mags[i]);
    if (v > max_mag) max_mag = v;
  }
  uint8_t num_bp_eq = 0;
  for (uint64_t v = max_mag; v != 0; v >>= 1) ++num_bp_eq;
  this->m_num_bitplanes = num_bp_eq;

  constexpr size_t hdr_size = 10;  // [num_bp:1] + [m_total_bits:8] + [fmt_ver:1]

  // All-zero short-circuit: just the header, no Huffman/ZSTD payload.
  if (max_mag == 0) {
    m_payload.assign(hdr_size, 0);
    m_payload[9] = kFormatVersion;
    this->m_total_bits = 0;
    return;
  }

  // ZigZag with int32 requires max_mag < 2^30 (so 2*max_mag < 2^31 stays in
  // signed int32 range). SPERR's realistic alphabets are << 2^30; this is a
  // defensive bound, not expected to fire.
  if (max_mag >= (uint64_t(1) << 30))
    throw std::runtime_error("HuffZstd: max_mag exceeds int32 ZigZag bound");

  // Build ZigZag bins + dense frequency array. Iterate subband-by-subband
  // (lowest frequency first) so that all-zero high-frequency subbands cluster
  // contiguously — gives the outer zstd long runs of the bin=0 code to LZ-match.
  const int32_t state_num = static_cast<int32_t>(2 * max_mag + 2);
  std::vector<int32_t> bins(N);
  std::vector<size_t> freq(static_cast<size_t>(state_num), 0);
  const auto boxes = enumerate_subbands(this->m_dims);
  const size_t Nx = this->m_dims[0];
  const size_t plane = Nx * this->m_dims[1];
  size_t pos = 0;
  for (const auto& bx : boxes) {
    for (size_t z = bx.za; z < bx.zb; ++z) {
      for (size_t y = bx.ya; y < bx.yb; ++y) {
        const size_t base = z * plane + y * Nx;
        for (size_t x = bx.xa; x < bx.xb; ++x) {
          const size_t idx = base + x;
          const int32_t mag = static_cast<int32_t>(static_cast<uint32_t>(mags[idx]));
          const int32_t s = signs.rbit(idx) ? -mag : mag;
          const int32_t b = zigzag_encode(s);
          bins[pos] = b;
          ++freq[static_cast<size_t>(b)];
          ++pos;
        }
      }
    }
  }
  assert(pos == N);

  // Optional: dump int32 ZigZag bins for offline analysis (qcat huffmanZstd
  // baseline comparison etc). Env-gated; per-chunk file suffix .cN where N is
  // a process-local atomic counter so concurrent OMP chunks don't collide.
  if (const char* prefix = std::getenv("SPERR_HZ_DUMP_BINS"); prefix && prefix[0]) {
    static std::atomic<int> g_dump_counter{0};
    const int idx = g_dump_counter.fetch_add(1);
    const std::string base = std::string(prefix) + ".c" + std::to_string(idx);
    std::ofstream bf(base + ".bin", std::ios::binary);
    bf.write(reinterpret_cast<const char*>(bins.data()),
             static_cast<std::streamsize>(bins.size() * sizeof(int32_t)));
    std::ofstream mf(base + ".meta");
    mf << "N=" << N << "\nstateNum=" << state_num << "\nmax_mag=" << max_mag << "\n";
  }

  auto zst_payload = encode_huffzstd(bins, freq, state_num);
  const uint64_t bits_for_full_len =
      static_cast<uint64_t>(zst_payload.size() + 1u) * 8u;  // +1 for fmt_ver byte
  this->m_total_bits = bits_for_full_len;

  m_payload.resize(hdr_size + zst_payload.size());
  m_payload[0] = num_bp_eq;
  std::memcpy(m_payload.data() + 1, &bits_for_full_len, sizeof(bits_for_full_len));
  m_payload[9] = kFormatVersion;
  std::memcpy(m_payload.data() + hdr_size, zst_payload.data(), zst_payload.size());
}

template <typename T>
void sperr::HuffZstd_INT<T>::use_bitstream(const void* p, size_t len)
{
  assert(len >= 10);
  const auto* p8 = static_cast<const uint8_t*>(p);

  std::memcpy(&this->m_num_bitplanes, p8, sizeof(this->m_num_bitplanes));
  std::memcpy(&this->m_total_bits, p8 + 1, sizeof(this->m_total_bits));
  this->m_avail_bits = this->m_total_bits;

  // Format-version gate. Reject anything we don't know how to parse.
  const uint8_t fmt = p8[9];
  if (fmt != kFormatVersion)
    throw std::runtime_error("HuffZstd: unknown payload format_version");

  m_payload.assign(p8, p8 + len);
}

template <typename T>
void sperr::HuffZstd_INT<T>::decode()
{
  const size_t N = this->m_dims[0] * this->m_dims[1] * this->m_dims[2];

  this->m_coeff_buf.assign(N, uint_type{0});
  this->m_sign_array.resize(N);
  this->m_sign_array.reset();

  if (this->m_num_bitplanes == 0)
    return;

  constexpr size_t hdr_size = 10;
  assert(m_payload.size() >= hdr_size);
  const uint8_t* zst_src = m_payload.data() + hdr_size;
  const size_t   zst_src_len = m_payload.size() - hdr_size;

  auto bins = decode_huffzstd(zst_src, zst_src_len, N);

  // Inverse subband permutation: bins[] are in subband-scan order; place each
  // value at its original raster position.
  const auto boxes = enumerate_subbands(this->m_dims);
  const size_t Nx = this->m_dims[0];
  const size_t plane = Nx * this->m_dims[1];
  size_t pos = 0;
  for (const auto& bx : boxes) {
    for (size_t z = bx.za; z < bx.zb; ++z) {
      for (size_t y = bx.ya; y < bx.yb; ++y) {
        const size_t base = z * plane + y * Nx;
        for (size_t x = bx.xa; x < bx.xb; ++x) {
          const size_t idx = base + x;
          const int32_t s = zigzag_decode(bins[pos]);
          if (s < 0) {
            this->m_sign_array.wtrue(idx);
            this->m_coeff_buf[idx] = static_cast<uint_type>(static_cast<uint32_t>(-s));
          }
          else {
            this->m_coeff_buf[idx] = static_cast<uint_type>(static_cast<uint32_t>(s));
          }
          ++pos;
        }
      }
    }
  }
  assert(pos == N);
}

template <typename T>
auto sperr::HuffZstd_INT<T>::encoded_bitstream_len() const -> size_t
{
  return m_payload.size();
}

template <typename T>
void sperr::HuffZstd_INT<T>::append_encoded_bitstream(vec8_type& buf) const
{
  buf.insert(buf.end(), m_payload.cbegin(), m_payload.cend());
}

// Explicit instantiations — match SPECK_FLT.
template class sperr::HuffZstd_INT<uint8_t>;
template class sperr::HuffZstd_INT<uint16_t>;
template class sperr::HuffZstd_INT<uint32_t>;
template class sperr::HuffZstd_INT<uint64_t>;
