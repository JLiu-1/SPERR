#include "LC_INT.h"
#include "sperr_helper.h"

// Vendored LC pipeline (BIT_4 → RZE_4 → RZE_1, in-memory wrapper).
#include "../third_party/lc/lc_pipeline.hpp"

// ZSTD for the outer wrap.
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "zstd.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace {

// Payload format version. byte 9 of m_payload carries this.
//   1 = ZigZag bins, subband-reorder permutation, LC pipeline BIT_4 RZE_4 RZE_1
//       wrapped in zstd. Same subband layout as HuffZstd_INT Sprint 5.
constexpr uint8_t kFormatVersion = 1;

struct SubbandBox {
  size_t xa, xb;
  size_t ya, yb;
  size_t za, zb;
};

// Same subband enumeration as HuffZstd_INT. We duplicate (rather than share)
// to keep this backend self-contained; the deterministic layout from chunk
// dims means encoder and decoder produce the same boxes.
auto enumerate_subbands(sperr::dims_type dims) -> std::vector<SubbandBox>
{
  std::vector<SubbandBox> out;
  auto dy = sperr::can_use_dyadic(dims);

  if (dy) {
    const size_t L = *dy;
    std::vector<size_t> ax(L + 1), ay(L + 1), az(L + 1);
    for (size_t k = 0; k <= L; ++k) {
      ax[k] = sperr::calc_approx_detail_len(dims[0], k)[0];
      ay[k] = sperr::calc_approx_detail_len(dims[1], k)[0];
      az[k] = sperr::calc_approx_detail_len(dims[2], k)[0];
    }
    out.reserve(1 + 7 * L);
    out.push_back({0, ax[L], 0, ay[L], 0, az[L]});
    for (size_t lev = L; lev > 0; --lev) {
      const size_t xS = ax[lev], xL = ax[lev - 1];
      const size_t yS = ay[lev], yL = ay[lev - 1];
      const size_t zS = az[lev], zL = az[lev - 1];
      out.push_back({xS, xL,  0, yS,  0, zS});
      out.push_back({ 0, xS, yS, yL,  0, zS});
      out.push_back({xS, xL, yS, yL,  0, zS});
      out.push_back({ 0, xS,  0, yS, zS, zL});
      out.push_back({xS, xL,  0, yS, zS, zL});
      out.push_back({ 0, xS, yS, yL, zS, zL});
      out.push_back({xS, xL, yS, yL, zS, zL});
    }
    return out;
  }

  const size_t Lz  = sperr::num_of_xforms(dims[2]);
  const size_t Lxy = sperr::num_of_xforms(std::min(dims[0], dims[1]));
  if (Lz == 0 && Lxy == 0) {
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

  std::vector<std::pair<size_t, size_t>> z_bands;
  z_bands.reserve(Lz + 1);
  z_bands.emplace_back(0, az[Lz]);
  for (size_t lev = Lz; lev > 0; --lev)
    z_bands.emplace_back(az[lev], az[lev - 1]);

  out.reserve((Lz + 1) * (1 + 3 * Lxy));
  for (const auto& [za, zb] : z_bands) {
    out.push_back({0, ax[Lxy], 0, ay[Lxy], za, zb});
    for (size_t lev = Lxy; lev > 0; --lev) {
      const size_t xS = ax[lev], xL = ax[lev - 1];
      const size_t yS = ay[lev], yL = ay[lev - 1];
      out.push_back({xS, xL,  0, yS, za, zb});
      out.push_back({ 0, xS, yS, yL, za, zb});
      out.push_back({xS, xL, yS, yL, za, zb});
    }
  }
  return out;
}

inline int32_t zigzag_encode(int32_t v) noexcept {
  return (v << 1) ^ (v >> 31);
}
inline int32_t zigzag_decode(int32_t b) noexcept {
  return (static_cast<uint32_t>(b) >> 1) ^ -(b & 1);
}

// Encode bins[] (ZigZag-encoded, non-negative int32) via LC pipeline, then zstd.
auto encode_lc_zstd(const std::vector<int32_t>& bins) -> sperr::vec8_type
{
  const size_t in_bytes = bins.size() * sizeof(int32_t);
  const long long lc_cap = sperr_lc::max_encoded_size(static_cast<long long>(in_bytes));
  std::vector<uint8_t> lc_buf(static_cast<size_t>(lc_cap));

  const long long lc_size = sperr_lc::encode(
      reinterpret_cast<const sperr_lc::byte*>(bins.data()),
      static_cast<long long>(in_bytes),
      lc_buf.data(),
      lc_cap);
  if (lc_size < 0)
    throw std::runtime_error("LC_INT: encode buffer too small");

  std::vector<SZ3::uchar> zst_buf(ZSTD_compressBound(static_cast<size_t>(lc_size)) + sizeof(size_t));
  SZ3::Lossless_zstd zstd;
  const size_t zst_len = zstd.compress(lc_buf.data(),
                                       static_cast<size_t>(lc_size),
                                       zst_buf.data(),
                                       zst_buf.size());

  sperr::vec8_type out(zst_len);
  std::memcpy(out.data(), zst_buf.data(), zst_len);
  return out;
}

// Decode zstd then LC. Returns the int32 bins (length N).
auto decode_lc_zstd(const uint8_t* zst_src, size_t zst_src_len, size_t N) -> std::vector<int32_t>
{
  SZ3::uchar* lc_buf = nullptr;
  size_t lc_len = 0;
  SZ3::Lossless_zstd zstd;
  zstd.decompress(const_cast<SZ3::uchar*>(zst_src), zst_src_len, lc_buf, lc_len);

  std::vector<int32_t> bins(N);
  const long long got = sperr_lc::decode(
      lc_buf,
      static_cast<long long>(lc_len),
      reinterpret_cast<sperr_lc::byte*>(bins.data()),
      static_cast<long long>(N * sizeof(int32_t)));
  std::free(lc_buf);
  if (got != static_cast<long long>(N * sizeof(int32_t)))
    throw std::runtime_error("LC_INT: decoded length mismatch");
  return bins;
}

}  // namespace

template <typename T>
void sperr::LC_INT<T>::encode()
{
  const auto& mags = this->m_coeff_buf;
  auto& signs = this->m_sign_array;
  const size_t N = mags.size();

  uint64_t max_mag = 0;
  for (size_t i = 0; i < N; ++i) {
    const uint64_t v = static_cast<uint64_t>(mags[i]);
    if (v > max_mag) max_mag = v;
  }
  uint8_t num_bp_eq = 0;
  for (uint64_t v = max_mag; v != 0; v >>= 1) ++num_bp_eq;
  this->m_num_bitplanes = num_bp_eq;

  constexpr size_t hdr_size = 10;

  if (max_mag == 0) {
    m_payload.assign(hdr_size, 0);
    m_payload[9] = kFormatVersion;
    this->m_total_bits = 0;
    return;
  }

  if (max_mag >= (uint64_t(1) << 30))
    throw std::runtime_error("LC_INT: max_mag exceeds int32 ZigZag bound");

  std::vector<int32_t> bins(N);
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
          bins[pos++] = zigzag_encode(s);
        }
      }
    }
  }
  assert(pos == N);

  auto zst_payload = encode_lc_zstd(bins);
  const uint64_t bits_for_full_len =
      static_cast<uint64_t>(zst_payload.size() + 1u) * 8u;
  this->m_total_bits = bits_for_full_len;

  m_payload.resize(hdr_size + zst_payload.size());
  m_payload[0] = num_bp_eq;
  std::memcpy(m_payload.data() + 1, &bits_for_full_len, sizeof(bits_for_full_len));
  m_payload[9] = kFormatVersion;
  std::memcpy(m_payload.data() + hdr_size, zst_payload.data(), zst_payload.size());
}

template <typename T>
void sperr::LC_INT<T>::use_bitstream(const void* p, size_t len)
{
  assert(len >= 10);
  const auto* p8 = static_cast<const uint8_t*>(p);

  std::memcpy(&this->m_num_bitplanes, p8, sizeof(this->m_num_bitplanes));
  std::memcpy(&this->m_total_bits, p8 + 1, sizeof(this->m_total_bits));
  this->m_avail_bits = this->m_total_bits;

  const uint8_t fmt = p8[9];
  if (fmt != kFormatVersion)
    throw std::runtime_error("LC_INT: unknown payload format_version");

  m_payload.assign(p8, p8 + len);
}

template <typename T>
void sperr::LC_INT<T>::decode()
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

  auto bins = decode_lc_zstd(zst_src, zst_src_len, N);

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
auto sperr::LC_INT<T>::encoded_bitstream_len() const -> size_t
{
  return m_payload.size();
}

template <typename T>
void sperr::LC_INT<T>::append_encoded_bitstream(vec8_type& buf) const
{
  buf.insert(buf.end(), m_payload.cbegin(), m_payload.cend());
}

template class sperr::LC_INT<uint8_t>;
template class sperr::LC_INT<uint16_t>;
template class sperr::LC_INT<uint32_t>;
template class sperr::LC_INT<uint64_t>;
