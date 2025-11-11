#include "CDF97_F.h"

#include <algorithm>
#include <cassert>
#include <numeric>  // std::accumulate()
#include <type_traits>

#ifdef __AVX2__
#include <immintrin.h>
#endif

// Destructor
sperr::CDF97_F::~CDF97_F()
{
  if (m_aligned_buf)
    std::free(m_aligned_buf);
}

template <typename T>
auto sperr::CDF97::copy_data(const T* data, size_t len, dims_type dims) -> RTNType
{
  static_assert(std::is_floating_point<T>::value, "!! Only floating point values are supported !!");
  if (len != dims[0] * dims[1] * dims[2])
    return RTNType::WrongLength;

  m_data_buf.resize(len);
  std::copy(data, data + len, m_data_buf.begin());

  m_dims = dims;

  auto max_col = std::max(std::max(dims[0], dims[1]), dims[2]);
  if (max_col * sizeof(float) > m_aligned_buf_bytes) {
    if (m_aligned_buf)
      std::free(m_aligned_buf);
    size_t alignment = 64;  // 256 bits
    size_t alloc_chunks = (max_col * 16 + 31) / alignment;
    m_aligned_buf_bytes = alignment * alloc_chunks;
    m_aligned_buf = static_cast<float*>(std::aligned_alloc(alignment, m_aligned_buf_bytes));
  }

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[0] * dims[2]), dims[1] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(max_slice);

  return RTNType::Good;
}
template auto sperr::CDF97::copy_data(const float*, size_t, dims_type) -> RTNType;
template auto sperr::CDF97::copy_data(const float*, size_t, dims_type) -> RTNType;

auto sperr::CDF97::take_data(vecd_type&& buf, dims_type dims) -> RTNType
{
  if (buf.size() != dims[0] * dims[1] * dims[2])
    return RTNType::WrongLength;

  m_data_buf.resize(buf.size());
  for (size_t i = 0; i < buf.size(); ++i) {
      m_data_buf[i] = static_cast<float>(buf[i]);
  }
  buf.clear();
  m_dims = dims;

  auto max_col = std::max(std::max(dims[0], dims[1]), dims[2]);
  if (max_col * sizeof(float) > m_aligned_buf_bytes) {
    if (m_aligned_buf)
      std::free(m_aligned_buf);
    size_t alignment = 64;  // 256 bits
    size_t alloc_chunks = (max_col * 16 + 31) / alignment;
    m_aligned_buf_bytes = alignment * alloc_chunks;
    m_aligned_buf = static_cast<float*>(std::aligned_alloc(alignment, m_aligned_buf_bytes));
  }

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[0] * dims[2]), dims[1] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(max_slice);

  return RTNType::Good;
}

auto sperr::CDF97::view_data() const -> const vecf_type&
{
  return m_data_buf;
}

auto sperr::CDF97::release_data() -> vecd_type
{
  vecd_type data(m_data_buf.size());
  for (size_t i = 0; i < m_data_buf.size(); ++i) {
      data[i] = static_cast<float>(m_data_buf[i]);
  }
  m_data_buf.clear();
  return data;
}

auto sperr::CDF97::get_dims() const -> std::array<size_t, 3>
{
  return m_dims;
}

void sperr::CDF97::dwt1d()
{
  auto num_xforms = sperr::num_of_xforms(m_dims[0]);
  m_dwt1d(m_data_buf.data(), m_data_buf.size(), num_xforms);
}

void sperr::CDF97::idwt1d()
{
  auto num_xforms = sperr::num_of_xforms(m_dims[0]);
  m_idwt1d(m_data_buf.data(), m_data_buf.size(), num_xforms);
}

void sperr::CDF97::dwt2d()
{
  auto xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  m_dwt2d(m_data_buf.data(), {m_dims[0], m_dims[1]}, xy);
}

void sperr::CDF97::idwt2d()
{
  auto xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  m_idwt2d(m_data_buf.data(), {m_dims[0], m_dims[1]}, xy);
}

auto sperr::CDF97_F::idwt2d_multi_res() -> std::vector<vecd_type>
{
  const auto xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  auto ret = std::vector<vecf_type>();

  if (xy > 0) {
    ret.reserve(xy);
    for (size_t lev = xy; lev > 0; lev--) {
      auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev);
      auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev);
      ret.emplace_back(m_sub_slice({x, y}));
      m_idwt2d_one_level(m_data_buf.copy_data(), {x + xd, y + yd});
    }
  }

  auto ret_d = std::vector<vecd_type>(ret.size());
  for(size_t i=0;i<ret.size();i++){
    ret_d[i].resize(ret[i].size());
    for(size_t j=0;j<ret[i].size();j++){
      ret_d[i][j]=ret[i][j];
    }
    ret[i].clear();
  }

  return ret_d;
}

void sperr::CDF97::dwt3d()
{
  auto dyadic = sperr::can_use_dyadic(m_dims);
  if (dyadic)
    m_dwt3d_dyadic(*dyadic);
  else
    m_dwt3d_wavelet_packet();
}

void sperr::CDF97::idwt3d()
{
  auto dyadic = sperr::can_use_dyadic(m_dims);
  if (dyadic)
    m_idwt3d_dyadic(*dyadic);
  else
    m_idwt3d_wavelet_packet();
}

void sperr::CDF97::idwt3d_multi_res(std::vector<vecd_type>& h_d)
{
  auto dyadic = sperr::can_use_dyadic(m_dims);
  std::vector<vecf_type> h;
  if (dyadic) {
    h.resize(*dyadic);
    for (size_t lev = *dyadic; lev > 0; lev--) {
      auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev);
      auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev);
      auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev);
      auto& buf = h[*dyadic - lev];
      buf.resize(x * y * z);
      m_sub_volume({x, y, z}, buf.begin());
      m_idwt3d_one_level(m_data_buf.begin(), {x + xd, y + yd, z + zd});
    }
  }
  else
    m_idwt3d_wavelet_packet();
  
  h_d.resize(h.size());
  for(size_t i=0;i<h.size();i++){
    h_d[i].resize(h[i].size());
    for(size_t j=0;j<h[i].size();j++){
      h_d[i][j]=h[i][j];
    }
    h[i].clear();
  }
}

void sperr::CDF97::m_dwt3d_wavelet_packet()
{
  /*
   *             Z
   *            /
   *           /
   *          /________
   *         /       /|
   *        /       / |
   *     0 |-------/-------> X
   *       |       |  |
   *       |       |  /
   *       |       | /
   *       |_______|/
   *       |
   *       |
   *       Y
   */

  const size_t plane_size_xy = m_dims[0] * m_dims[1];

  // First transform along the Z dimension
  //
  const auto num_xforms_z = sperr::num_of_xforms(m_dims[2]);

  for (size_t y = 0; y < m_dims[1]; y++) {
    const auto y_offset = y * m_dims[0];

    // Re-arrange values of one XZ slice so that they form many z_columns
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_slice_buf[z + x * m_dims[2]] = m_data_buf[cube_start_idx + x];
    }

    // DWT1D on every z_column
    for (size_t x = 0; x < m_dims[0]; x++)
      m_dwt1d(m_slice_buf.data() + x * m_dims[2], m_dims[2], num_xforms_z);

    // Put back values of the z_columns to the cube
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_data_buf[cube_start_idx + x] = m_slice_buf[z + x * m_dims[2]];
    }
  }

  // Second transform each plane
  //
  const auto num_xforms_xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));

  for (size_t z = 0; z < m_dims[2]; z++) {
    const size_t offset = plane_size_xy * z;
    m_dwt2d(m_data_buf.data() + offset, {m_dims[0], m_dims[1]}, num_xforms_xy);
  }
}

void sperr::CDF97::m_idwt3d_wavelet_packet()
{
  const size_t plane_size_xy = m_dims[0] * m_dims[1];

  // First, inverse transform each plane
  //
  auto num_xforms_xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  for (size_t i = 0; i < m_dims[2]; i++) {
    const size_t offset = plane_size_xy * i;
    m_idwt2d(m_data_buf.data() + offset, {m_dims[0], m_dims[1]}, num_xforms_xy);
  }

  /*
   * Second, inverse transform along the Z dimension
   *
   *             Z
   *            /
   *           /
   *          /________
   *         /       /|
   *        /       / |
   *     0 |-------/-------> X
   *       |       |  |
   *       |       |  /
   *       |       | /
   *       |_______|/
   *       |
   *       |
   *       Y
   */

  // Process one XZ slice at a time
  //
  const auto num_xforms_z = sperr::num_of_xforms(m_dims[2]);
  for (size_t y = 0; y < m_dims[1]; y++) {
    const auto y_offset = y * m_dims[0];

    // Re-arrange values on one slice so that they form many z_columns
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_slice_buf[z + x * m_dims[2]] = m_data_buf[cube_start_idx + x];
    }

    // IDWT1D on every z_column
    for (size_t x = 0; x < m_dims[0]; x++)
      m_idwt1d(m_slice_buf.data() + x * m_dims[2], m_dims[2], num_xforms_z);

    // Put back values from the z_columns to the cube
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_data_buf[cube_start_idx + x] = m_slice_buf[z + x * m_dims[2]];
    }
  }
}

void sperr::CDF97::m_dwt3d_dyadic(size_t num_xforms)
{
  for (size_t lev = 0; lev < num_xforms; lev++) {
    auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev);
    auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev);
    auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev);
    m_dwt3d_one_level({x, y, z});
  }
}

void sperr::CDF97::m_idwt3d_dyadic(size_t num_xforms)
{
  for (size_t lev = num_xforms; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev - 1);
    auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev - 1);
    auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev - 1);
    m_idwt3d_one_level({x, y, z});
  }
}

//
// Private Methods
//
void sperr::CDF97::m_dwt1d(float* array, size_t array_len, size_t num_of_lev)
{
  for (size_t lev = 0; lev < num_of_lev; lev++) {
    m_gather(array, array_len, m_aligned_buf);
    this->QccWAVCDF97AnalysisSymmetric(m_aligned_buf, array_len);
    std::copy(m_aligned_buf, m_aligned_buf + array_len, array);
    array_len -= array_len / 2;
  }
}

void sperr::CDF97::m_idwt1d(float* array, size_t array_len, size_t num_of_lev)
{
  for (size_t lev = num_of_lev; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(array_len, lev - 1);
    this->QccWAVCDF97SynthesisSymmetric(array, x);
    m_scatter(array, x, m_aligned_buf);
    std::copy(m_aligned_buf, m_aligned_buf + x, array);
  }
}

void sperr::CDF97::m_dwt2d(float* plane, std::array<size_t, 2> len_xy, size_t num_of_lev)
{
  for (size_t lev = 0; lev < num_of_lev; lev++) {
    auto [x, xd] = sperr::calc_approx_detail_len(len_xy[0], lev);
    auto [y, yd] = sperr::calc_approx_detail_len(len_xy[1], lev);
    m_dwt2d_one_level(plane, {x, y});
  }
}

void sperr::CDF97::m_idwt2d(float* plane, std::array<size_t, 2> len_xy, size_t num_of_lev)
{
  for (size_t lev = num_of_lev; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(len_xy[0], lev - 1);
    auto [y, yd] = sperr::calc_approx_detail_len(len_xy[1], lev - 1);
    m_idwt2d_one_level(plane, {x, y});
  }
}

void sperr::CDF97::m_dwt2d_one_level(float* plane, std::array<size_t, 2> len_xy)
{
  // First, perform DWT along X for every row
  for (size_t i = 0; i < len_xy[1]; i++) {
    auto* pos = plane + i * m_dims[0];
    m_gather(pos, len_xy[0], m_aligned_buf);
    this->QccWAVCDF97AnalysisSymmetric(m_aligned_buf, len_xy[0]);
    std::copy(m_aligned_buf, m_aligned_buf + len_xy[0], pos);
  }

  // Second, perform DWT along Y for every column
  for (size_t x = 0; x < len_xy[0]; x++) {
    for (size_t y = 0; y < len_xy[1]; y++)
      m_slice_buf[y] = plane[y * m_dims[0] + x];
    m_gather(m_slice_buf.data(), len_xy[1], m_aligned_buf);
    this->QccWAVCDF97AnalysisSymmetric(m_aligned_buf, len_xy[1]);
    for (size_t y = 0; y < len_xy[1]; y++)
      plane[y * m_dims[0] + x] = m_aligned_buf[y];
  }
}

void sperr::CDF97::m_idwt2d_one_level(float* plane, std::array<size_t, 2> len_xy)
{
  // First, perform IDWT along Y for every column
  for (size_t x = 0; x < len_xy[0]; x++) {
    for (size_t y = 0; y < len_xy[1]; y++)
      m_slice_buf[y] = plane[y * m_dims[0] + x];
    this->QccWAVCDF97SynthesisSymmetric(m_slice_buf.data(), len_xy[1]);
    m_scatter(m_slice_buf.data(), len_xy[1], m_aligned_buf);
    for (size_t y = 0; y < len_xy[1]; y++)
      plane[y * m_dims[0] + x] = m_aligned_buf[y];
  }

  // Second, perform IDWT along X for every row
  for (size_t i = 0; i < len_xy[1]; i++) {
    auto* pos = plane + i * m_dims[0];
    this->QccWAVCDF97SynthesisSymmetric(pos, len_xy[0]);
    m_scatter(pos, len_xy[0], m_aligned_buf);
    std::copy(m_aligned_buf, m_aligned_buf + len_xy[0], pos);
  }
}

void sperr::CDF97::m_dwt3d_one_level(std::array<size_t, 3> len_xyz)
{
  // First, do one level of transform on all XY planes.
  const auto plane_size_xy = m_dims[0] * m_dims[1];
  const auto col_len = len_xyz[2];
  for (size_t z = 0; z < col_len; z++) {
    const size_t offset = plane_size_xy * z;
    m_dwt2d_one_level(m_data_buf.data() + offset, {len_xyz[0], len_xyz[1]});
  }

  // Second, do one level of transform on all Z columns.  Strategy:
  // 1) extract eight Z columns to buffer space `m_slice_buf`
  // 2) transform these eight columns
  // 3) put the Z columns back to their locations in the volume.
  //
  // Note: the reason to process 16 columns at a time is that a cache line
  // is usually 64 bytes, or 16 floats. That means when you pay the cost to retrieve
  // one value from the Z column, its neighboring 15 values are available for free!

  for (size_t y = 0; y < len_xyz[1]; y++) {
    for (size_t x = 0; x < len_xyz[0]; x += 8) {
      const size_t xy_offset = y * m_dims[0] + x;
      const auto stride = std::min(16ul, len_xyz[0] - x);

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_slice_buf[z + i * col_len] = m_data_buf[z * plane_size_xy + xy_offset + i];
      }

      for (size_t i = 0; i < stride; i++) {
        auto* itr = m_slice_buf.data() + i * col_len;
        m_gather(itr, col_len, m_aligned_buf);
        this->QccWAVCDF97AnalysisSymmetric(m_aligned_buf, col_len);
        std::copy(m_aligned_buf, m_aligned_buf + col_len, itr);
      }

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_data_buf[z * plane_size_xy + xy_offset + i] = m_slice_buf[z + i * col_len];
      }
    }
  }
}

void sperr::CDF97::m_idwt3d_one_level(std::array<size_t, 3> len_xyz)
{
  const auto plane_size_xy = m_dims[0] * m_dims[1];
  const auto col_len = len_xyz[2];

  // First, do one level of inverse transform on all Z columns.  Strategy:
  // 1) extract eight Z columns to buffer space `m_slice_buf`
  // 2) transform these eight columns
  // 3) put the Z columns back to their appropriate locations in the volume.
  //
  // Note: the reason to process 16 columns at a time is that a cache line
  // is usually 64 bytes, or 16 floats. That means when you pay the cost to retrieve
  // one value from the Z column, its neighboring 15 values are available for free!

  for (size_t y = 0; y < len_xyz[1]; y++) {
    for (size_t x = 0; x < len_xyz[0]; x += 8) {
      const size_t xy_offset = y * m_dims[0] + x;
      const auto stride = std::min(16ul, len_xyz[0] - x);

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_slice_buf[z + i * col_len] = m_data_buf[z * plane_size_xy + xy_offset + i];
      }

      for (size_t i = 0; i < stride; i++) {
        auto* itr = m_slice_buf.data() + i * col_len;
        this->QccWAVCDF97SynthesisSymmetric(itr, col_len);
        m_scatter(itr, col_len, m_aligned_buf);
        std::copy(m_aligned_buf, m_aligned_buf + col_len, itr);
      }

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_data_buf[z * plane_size_xy + xy_offset + i] = m_slice_buf[z + i * col_len];
      }
    }
  }

  // Second, do one level of inverse transform on all XY planes.
  for (size_t z = 0; z < len_xyz[2]; z++) {
    const size_t offset = plane_size_xy * z;
    m_idwt2d_one_level(m_data_buf.data() + offset, {len_xyz[0], len_xyz[1]});
  }
}

void sperr::CDF97::m_gather(const float* src, size_t len, float* dst) const
{
#ifdef __AVX2__
  const float* src_end = src + len;
  float* dst_evens = dst;
  float* dst_odds = dst + len - len / 2;

  // Process 16 elements at a time
  for (; src + 16 <= src_end; src += 16) {
    __m256 v0 = _mm256_loadu_ps(src);      //  0 1 2 3 4 5 6 7
    __m256 v1 = _mm256_loadu_ps(src + 8);  //  8 9 10 11 12 13 14 15

    __m256 lo = _mm256_unpacklo_ps(v0, v1);  // 0 8 1 9 2 10 3 11
    __m256 hi = _mm256_unpackhi_ps(v0, v1);  // 4 12 5 13 6 14 7 15


    __m256 even_part1 = _mm256_shuffle_ps(lo, hi, 0b10001000); // 0 2 4 6
    __m256 even_part2 = _mm256_shuffle_ps(lo, hi, 0b11011101); // 8 10 12 14
    __m256 evens = _mm256_permute2f128_ps(even_part1, even_part2, 0x20); // 拼接成 0 2 4 6 8 10 12 14

    __m256 odd_part1 = _mm256_shuffle_ps(lo, hi, 0b11011101);  // 1 3 5 7
    __m256 odd_part2 = _mm256_shuffle_ps(lo, hi, 0b11111111);  // 9 11 13 15
    __m256 odds = _mm256_permute2f128_ps(odd_part1, odd_part2, 0x20); // 拼接成 1 3 5 7 9 11 13 15

    _mm256_storeu_ps(dst_evens, evens);
    _mm256_storeu_ps(dst_odds, odds);

    dst_evens += 8;
    dst_odds += 8;
  }
  for (; src < src_end - 1; src += 2) {
    *(dst_evens++) = *src;
    *(dst_odds++) = *(src + 1);
  }

  if (src < src_end)
    *dst_evens = *src;
#else
  size_t low_count = len - len / 2, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *dst = *(src + i * 2);
    ++dst;
  }
  for (size_t i = 0; i < high_count; i++) {
    *dst = *(src + i * 2 + 1);
    ++dst;
  }
#endif
}

void sperr::CDF97::m_scatter(const float* begin, size_t len, float* dst) const
{
#ifdef __AVX2__
  const float* even_end = begin + len - len / 2;
  const float* odd_beg = even_end;
  const float* dst_end = dst + len;

  // Process 16 elements at a time
  for (; begin + 8 < even_end; begin += 8) {
    __m256 v0 = _mm256_loadu_ps(begin);    // e0, e1, e2, e3, e4, e5, e6, e7
    __m256 v1 = _mm256_loadu_ps(odd_beg);  // o0, o1, o2, o3, o4, o5, o6, o7

    __m256 lo = _mm256_unpacklo_ps(v0, v1);
    __m256 hi = _mm256_unpackhi_ps(v0, v1);

    __m256 result1 = _mm256_permute2f128_ps(lo, hi, 0x20);
    __m256 result2 = _mm256_permute2f128_ps(lo, hi, 0x31);

    _mm256_store_ps(dst,     result1);
    _mm256_store_ps(dst + 8, result2);

    dst     += 16;
    odd_beg += 8;
  }


  for (; dst < dst_end - 1; dst += 2) {
    *dst = *(begin++);
    *(dst + 1) = *(odd_beg++);
  }

  if (dst < dst_end)
    *dst = *begin;
#else
  size_t low_count = len - len / 2, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *(dst + i * 2) = *begin;
    ++begin;
  }
  for (size_t i = 0; i < high_count; i++) {
    *(dst + i * 2 + 1) = *begin;
    ++begin;
  }
#endif
}

auto sperr::CDF97::m_sub_slice(std::array<size_t, 2> subdims) const -> vecf_type
{
  assert(subdims[0] <= m_dims[0] && subdims[1] <= m_dims[1]);

  auto ret = vecf_type(subdims[0] * subdims[1]);
  auto dst = ret.begin();
  for (size_t y = 0; y < subdims[1]; y++) {
    auto beg = m_data_buf.begin() + y * m_dims[0];
    std::copy(beg, beg + subdims[0], dst);
    dst += subdims[0];
  }

  return ret;
}

void sperr::CDF97::m_sub_volume(dims_type subdims, float* dst) const
{
  assert(subdims[0] <= m_dims[0] && subdims[1] <= m_dims[1] && subdims[2] <= m_dims[2]);

  const auto slice_len = m_dims[0] * m_dims[1];
  for (size_t z = 0; z < subdims[2]; z++) {
    for (size_t y = 0; y < subdims[1]; y++) {
      auto beg = m_data_buf.begin() + z * slice_len + y * m_dims[0];
      std::copy(beg, beg + subdims[0], dst);
      dst += subdims[0];
    }
  }
}

//
// Methods from QccPack
//
void sperr::CDF97::QccWAVCDF97AnalysisSymmetric(float* signal, size_t len)
{
  size_t even_len = len - len / 2;
  size_t odd_len = len / 2;
  float* even = signal;
  float* odd = signal + even_len;

  // Process all the odd elements
  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] += ALPHA * (even[i] + even[i + 1]);
  odd[odd_len - 1] += ALPHA * (even[odd_len - 1] + even[even_len - 1]);

  // Process all the even elements
  even[0] += 2.0 * BETA * odd[0];
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] += BETA * (odd[i - 1] + odd[i]);
  even[even_len - 1] += BETA * (odd[even_len - 2] + odd[odd_len - 1]);

  // Process all the odd elements
  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] += GAMMA * (even[i] + even[i + 1]);
  odd[odd_len - 1] += GAMMA * (even[odd_len - 1] + even[even_len - 1]);

  // Process even elements
  even[0] = EPSILON * (even[0] + 2.0 * DELTA * odd[0]);
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] = EPSILON * (even[i] + DELTA * (odd[i - 1] + odd[i]));
  even[even_len - 1] =
      EPSILON * (even[even_len - 1] + DELTA * (odd[even_len - 2] + odd[odd_len - 1]));

  // Process odd elements
  for (size_t i = 0; i < odd_len; i++)
    odd[i] *= -INV_EPSILON;
}

void sperr::CDF97::QccWAVCDF97SynthesisSymmetric(float* signal, size_t len)
{
  size_t even_len = len - len / 2;
  size_t odd_len = len / 2;
  float* even = signal;
  float* odd = signal + even_len;

  // Process odd elements
  for (size_t i = 0; i < odd_len; i++)
    odd[i] *= (-EPSILON);

  // Process even elements
  even[0] = even[0] * INV_EPSILON - 2.0 * DELTA * odd[0];
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] = even[i] * INV_EPSILON - DELTA * (odd[i - 1] + odd[i]);
  even[even_len - 1] =
      even[even_len - 1] * INV_EPSILON - DELTA * (odd[even_len - 2] + odd[odd_len - 1]);

  // Process odd elements
  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] -= GAMMA * (even[i] + even[i + 1]);
  odd[odd_len - 1] -= GAMMA * (even[odd_len - 1] + even[even_len - 1]);

  // Process even elements
  even[0] -= 2.0 * BETA * odd[0];
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] -= BETA * (odd[i - 1] + odd[i]);
  even[even_len - 1] -= BETA * (odd[even_len - 2] + odd[odd_len - 1]);

  // Process odd elements
  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] -= ALPHA * (even[i] + even[i + 1]);
  odd[odd_len - 1] -= ALPHA * (even[odd_len - 1] + even[even_len - 1]);
}