#include "CDF97.h"

#include <algorithm>
#include <cassert>
#include <numeric>  // std::accumulate()
#include <type_traits>

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
  if (max_col * 2 > m_qcc_buf.size())
    m_qcc_buf.resize(std::max(m_qcc_buf.size(), max_col) * 2);

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[1] * dims[2]), dims[0] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(std::max(m_slice_buf.size(), max_slice));

  return RTNType::Good;
}

template auto sperr::CDF97::copy_data(const float*, size_t, dims_type) -> RTNType;
template auto sperr::CDF97::copy_data(const double*, size_t, dims_type) -> RTNType;

auto sperr::CDF97::take_data(vecd_type&& buf, dims_type dims) -> RTNType
{
  if (buf.size() != dims[0] * dims[1] * dims[2])
    return RTNType::WrongLength;

  m_data_buf = std::move(buf);
  m_dims = dims;

  auto max_col = std::max(std::max(dims[0], dims[1]), dims[2]);
  if (max_col * 2 > m_qcc_buf.size())
    m_qcc_buf.resize(std::max(m_qcc_buf.size(), max_col) * 2);

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[1] * dims[2]), dims[0] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(std::max(m_slice_buf.size(), max_slice));

  return RTNType::Good;
}

auto sperr::CDF97::release_data() -> vecd_type&&
{
  return std::move(m_data_buf);
}


auto sperr::CDF97::view_data() const -> const vecd_type&
{
  return m_data_buf;
}

auto sperr::CDF97::get_dims() const -> dims_type
{
  return m_dims;
}

void sperr::CDF97::dwt1d()
{
  auto n = std::min(m_dims[0], std::min(m_dims[1], m_dims[2]));
  auto num_xforms = sperr::num_of_xforms(n);
  m_dwt1d(m_data_buf.begin(), m_data_buf.size(), num_xforms);
}

void sperr::CDF97::idwt1d()
{
  auto n = std::min(m_dims[0], std::min(m_dims[1], m_dims[2]));
  auto num_xforms = sperr::num_of_xforms(n);
  m_idwt1d(m_data_buf.begin(), m_data_buf.size(), num_xforms);
}

void sperr::CDF97::dwt2d()
{
  auto n = std::min(m_dims[0], m_dims[1]);
  auto num_xforms = sperr::num_of_xforms(n);
  m_dwt2d(m_data_buf.begin(), {m_dims[0], m_dims[1]}, num_xforms);
}

void sperr::CDF97::idwt2d()
{
  auto n = std::min(m_dims[0], m_dims[1]);
  auto num_xforms = sperr::num_of_xforms(n);
  m_idwt2d(m_data_buf.begin(), {m_dims[0], m_dims[1]}, num_xforms);
}

auto sperr::CDF97::idwt2d_multi_res() -> std::vector<vecd_type>
{
  const auto xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  auto ret = std::vector<vecd_type>();

  if (xy > 0) {
    ret.reserve(xy);
    for (size_t lev = xy; lev > 0; lev--) {
      auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev);
      auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev);
      ret.emplace_back(m_sub_slice({x, y}));
      m_idwt2d_one_level(m_data_buf.begin(), {x + xd, y + yd});
    }
  }

  return ret;
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

void sperr::CDF97::idwt3d_multi_res(std::vector<vecd_type>& h)
{
  auto dyadic = sperr::can_use_dyadic(m_dims);

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
}

void sperr::CDF97::m_dwt1d(itd_type array, size_t array_len, size_t num_of_lev)
{
  for (size_t lev = 0; lev < num_of_lev; lev++) {
    auto [x, xd] = sperr::calc_approx_detail_len(array_len, lev);
    m_dwt1d_one_level(array, x);
  }
}

void sperr::CDF97::m_idwt1d(itd_type array, size_t array_len, size_t num_of_lev)
{
  for (size_t lev = num_of_lev; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(array_len, lev - 1);
    m_idwt1d_one_level(array, x);
  }
}

void sperr::CDF97::m_dwt2d(itd_type plane, std::array<size_t, 2> len_xy, size_t num_of_lev)
{
  for (size_t lev = 0; lev < num_of_lev; lev++) {
    auto [x, xd] = sperr::calc_approx_detail_len(len_xy[0], lev);
    auto [y, yd] = sperr::calc_approx_detail_len(len_xy[1], lev);
    m_dwt2d_one_level(plane, {x, y});
  }
}

void sperr::CDF97::m_idwt2d(itd_type plane, std::array<size_t, 2> len_xy, size_t num_of_lev)
{
  for (size_t lev = num_of_lev; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(len_xy[0], lev - 1);
    auto [y, yd] = sperr::calc_approx_detail_len(len_xy[1], lev - 1);
    m_idwt2d_one_level(plane, {x, y});
  }
}

void sperr::CDF97::m_dwt1d_one_level(itd_type array, size_t array_len)
{
  if (array_len <= 1)
    return;

  if (m_qcc_buf.size() < array_len)
    m_qcc_buf.resize(array_len);

  // Copy into contiguous buffer for lifting
  std::copy(array, array + array_len, m_qcc_buf.begin());

  if (array_len % 2 == 0) {
    QccWAVCDF97AnalysisSymmetricEvenEven(m_qcc_buf.data(), array_len);
    m_gather_even(m_qcc_buf.cbegin(), m_qcc_buf.cbegin() + array_len, array);
  }
  else {
    QccWAVCDF97AnalysisSymmetricOddEven(m_qcc_buf.data(), array_len);
    m_gather_odd(m_qcc_buf.cbegin(), m_qcc_buf.cbegin() + array_len, array);
  }
}

void sperr::CDF97::m_dwt1d_one_level_strided(itd_type base, size_t len, ptrdiff_t stride)
{
  if (len <= 1)
    return;

  if (m_qcc_buf.size() < len)
    m_qcc_buf.resize(len);

  if (len % 2 == 0) {
    QccWAVCDF97AnalysisSymmetricEvenEvenStrided(base->data(), len, stride);
  }
  else {
    QccWAVCDF97AnalysisSymmetricOddEvenStrided(base->data(), len, stride);
  }

  // Pack [even, odd] into [L | H] layout along the same strided line.
  const size_t low_count  = (len + 1) / 2;
  const size_t high_count = len / 2;
  double* tmp = m_qcc_buf.data();

  for (size_t i = 0; i < low_count; ++i)
    tmp[i] = base[(2 * i) * stride];

  for (size_t i = 0; i < high_count; ++i)
    tmp[low_count + i] = base[(2 * i + 1) * stride];

  for (size_t i = 0; i < len; ++i)
    base[i * stride] = tmp[i];
}

void sperr::CDF97::m_idwt1d_one_level(itd_type array, size_t array_len)
{
  if (array_len <= 1)
    return;

  if (m_qcc_buf.size() < array_len)
    m_qcc_buf.resize(array_len);

  if (array_len % 2 == 0) {
    m_scatter_even(array, array + array_len, m_qcc_buf.begin());
    QccWAVCDF97SynthesisSymmetricEvenEven(m_qcc_buf.data(), array_len);
  }
  else {
    m_scatter_odd(array, array + array_len, m_qcc_buf.begin());
    QccWAVCDF97SynthesisSymmetricOddEven(m_qcc_buf.data(), array_len);
  }

  std::copy(m_qcc_buf.cbegin(), m_qcc_buf.cbegin() + array_len, array);
}

void sperr::CDF97::m_idwt1d_one_level_strided(itd_type base, size_t len, ptrdiff_t stride)
{
  if (len <= 1)
    return;

  if (m_qcc_buf.size() < len)
    m_qcc_buf.resize(len);

  const size_t low_count  = (len + 1) / 2;
  const size_t high_count = len / 2;

  // Read [L | H] from strided storage into contiguous buffer
  for (size_t i = 0; i < len; ++i)
    m_qcc_buf[i] = base[i * stride];

  // Scatter into even/odd interleaving expected by lifting synthesis,
  // then apply the appropriate symmetric extension kernel in-place.
  if (len % 2 == 0) {
    // even-even
    for (size_t i = 0; i < low_count; ++i)
      base[(2 * i) * stride] = m_qcc_buf[i];
    for (size_t i = 0; i < high_count; ++i)
      base[(2 * i + 1) * stride] = m_qcc_buf[low_count + i];

    QccWAVCDF97SynthesisSymmetricEvenEvenStrided(base, len, stride);
  }
  else {
    // odd-even
    for (size_t i = 0; i < low_count; ++i)
      base[(2 * i) * stride] = m_qcc_buf[i];
    for (size_t i = 0; i < high_count; ++i)
      base[(2 * i + 1) * stride] = m_qcc_buf[low_count + i];

    QccWAVCDF97SynthesisSymmetricOddEvenStrided(base, len, stride);
  }
}

void sperr::CDF97::m_dwt2d_one_level(itd_type plane, std::array<size_t, 2> len_xy)
{
  // X is contiguous (stride = 1), Y uses stride = m_dims[0].
  if (len_xy[0] <= 1 && len_xy[1] <= 1)
    return;

  const auto max_len = std::max(len_xy[0], len_xy[1]);
  if (m_qcc_buf.size() < max_len)
    m_qcc_buf.resize(max_len);

  // 1) DWT along X for each row (contiguous).
  for (size_t y = 0; y < len_xy[1]; ++y) {
    double* row = plane + y * m_dims[0];
    m_dwt1d_one_level(row, len_xy[0]);
  }

  // 2) DWT along Y for each column (stride-aware).
  const ptrdiff_t stride_y = static_cast<ptrdiff_t>(m_dims[0]);
  for (size_t x = 0; x < len_xy[0]; ++x) {
    double* col0 = plane + x;  // (x,0)
    m_dwt1d_one_level_strided(col0, len_xy[1], stride_y);
  }
}

void sperr::CDF97::m_idwt2d_one_level(itd_type plane, std::array<size_t, 2> len_xy)
{
  if (len_xy[0] <= 1 && len_xy[1] <= 1)
    return;

  const auto max_len = std::max(len_xy[0], len_xy[1]);
  if (m_qcc_buf.size() < max_len)
    m_qcc_buf.resize(max_len);

  // 1) IDWT along Y for each column (stride-aware).
  const ptrdiff_t stride_y = static_cast<ptrdiff_t>(m_dims[0]);
  for (size_t x = 0; x < len_xy[0]; ++x) {
    double* col0 = plane + x;  // (x,0)
    m_idwt1d_one_level_strided(col0, len_xy[1], stride_y);
  }

  // 2) IDWT along X for each row (contiguous).
  for (size_t y = 0; y < len_xy[1]; ++y) {
    double* row = plane + y * m_dims[0];
    m_idwt1d_one_level(row, len_xy[0]);
  }
}

void sperr::CDF97::m_dwt3d_dyadic(size_t num_xforms)
{
  const auto plane_size_xy = m_dims[0] * m_dims[1];
  const auto plane_size_xz = m_dims[0] * m_dims[2];
  const auto plane_size_yz = m_dims[1] * m_dims[2];

  for (size_t lev = 0; lev < num_xforms; lev++) {
    auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev);
    auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev);
    auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev);

    // XY planes
    for (size_t iz = 0; iz < z; iz++) {
      auto offset = iz * plane_size_xy;
      m_dwt2d_one_level(m_data_buf.begin() + offset, {x, y});
    }

    // XZ planes
    for (size_t iy = 0; iy < y; iy++) {
      auto offset = iy * m_dims[0];
      auto slice = m_slice_buf.begin();
      for (size_t iz = 0; iz < z; iz++) {
        auto beg = m_data_buf.begin() + iz * plane_size_xy + offset;
        std::copy(beg, beg + x, slice);
        slice += x;
      }
      m_dwt2d_one_level(m_slice_buf.begin(), {x, z});
      slice = m_slice_buf.begin();
      for (size_t iz = 0; iz < z; iz++) {
        auto beg = m_data_buf.begin() + iz * plane_size_xy + offset;
        std::copy(slice, slice + x, beg);
        slice += x;
      }
    }

    // YZ planes
    for (size_t ix = 0; ix < x; ix++) {
      auto slice = m_slice_buf.begin();
      for (size_t iz = 0; iz < z; iz++) {
        auto beg = m_data_buf.begin() + iz * plane_size_xy + ix;
        for (size_t iy = 0; iy < y; iy++)
          *(slice++) = *(beg + iy * m_dims[0]);
      }
      m_dwt2d_one_level(m_slice_buf.begin(), {y, z});
      slice = m_slice_buf.begin();
      for (size_t iz = 0; iz < z; iz++) {
        auto beg = m_data_buf.begin() + iz * plane_size_xy + ix;
        for (size_t iy = 0; iy < y; iy++)
          *(beg + iy * m_dims[0]) = *(slice++);
      }
    }
  }
}

void sperr::CDF97::m_dwt3d_wavelet_packet()
{
  const size_t plane_size_xy = m_dims[0] * m_dims[1];

  // First, transform each plane (2D)
  auto num_xforms_xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  for (size_t i = 0; i < m_dims[2]; i++) {
    const size_t offset = plane_size_xy * i;
    m_dwt2d(m_data_buf.begin() + offset, {m_dims[0], m_dims[1]}, num_xforms_xy);
  }

  // Second, transform along Z dimension (wavelet packet style)
  const auto num_xforms_z = sperr::num_of_xforms(m_dims[2]);

  for (size_t lev = 0; lev < num_xforms_z; lev++) {
    auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev);

    // Process one XZ slice at a time
    auto slice_begin = m_slice_buf.begin();
    for (size_t y = 0; y < m_dims[1]; y++) {
      // Copy one XZ slice
      auto dst = slice_begin;
      for (size_t zz = 0; zz < z; zz++) {
        auto src = m_data_buf.begin() + zz * plane_size_xy + y * m_dims[0];
        std::copy(src, src + m_dims[0], dst);
        dst += m_dims[0];
      }

      // 2D DWT on that slice (X-Z)
      m_dwt2d(slice_begin, {m_dims[0], z}, 1);

      // Copy back
      dst = slice_begin;
      for (size_t zz = 0; zz < z; zz++) {
        auto src = m_data_buf.begin() + zz * plane_size_xy + y * m_dims[0];
        std::copy(dst, dst + m_dims[0], src);
        dst += m_dims[0];
      }
    }
  }
}

void sperr::CDF97::m_dwt3d_one_level(itd_type vol, std::array<size_t, 3> len_xyz)
{
  const size_t plane_size_xy = m_dims[0] * m_dims[1];

  // Ensure workspace for the longest line.
  const auto max_len = std::max({len_xyz[0], len_xyz[1], len_xyz[2]});
  if (m_qcc_buf.size() < max_len)
    m_qcc_buf.resize(max_len);

  // 1) One level of 2D transform on all XY planes.
  for (size_t z = 0; z < len_xyz[2]; ++z) {
    const size_t offset = plane_size_xy * z;
    m_dwt2d_one_level(vol + offset, {len_xyz[0], len_xyz[1]});
  }

  // 2) One level along Z for every (x,y) column using stride-aware lifting.
  const ptrdiff_t stride_z = static_cast<ptrdiff_t>(plane_size_xy);
  for (size_t y = 0; y < len_xyz[1]; ++y) {
    for (size_t x = 0; x < len_xyz[0]; ++x) {
      double* col0 = vol + y * m_dims[0] + x;  // (x,y,0)
      m_dwt1d_one_level_strided(col0, len_xyz[2], stride_z);
    }
  }
}

void sperr::CDF97::m_idwt3d_dyadic(size_t num_xforms)
{
  const auto plane_size_xy = m_dims[0] * m_dims[1];
  const auto plane_size_xz = m_dims[0] * m_dims[2];
  const auto plane_size_yz = m_dims[1] * m_dims[2];

  for (size_t lev = num_xforms; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev - 1);
    auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev - 1);
    auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev - 1);

    // YZ planes
    for (size_t ix = 0; ix < x; ix++) {
      auto slice = m_slice_buf.begin();
      for (size_t iz = 0; iz < z; iz++) {
        auto beg = m_data_buf.begin() + iz * plane_size_xy + ix;
        for (size_t iy = 0; iy < y; iy++)
          *(slice++) = *(beg + iy * m_dims[0]);
      }
      m_idwt2d_one_level(m_slice_buf.begin(), {y, z});
      slice = m_slice_buf.begin();
      for (size_t iz = 0; iz < z; iz++) {
        auto beg = m_data_buf.begin() + iz * plane_size_xy + ix;
        for (size_t iy = 0; iy < y; iy++)
          *(beg + iy * m_dims[0]) = *(slice++);
      }
    }

    // XZ planes
    for (size_t iy = 0; iy < y; iy++) {
      auto offset = iy * m_dims[0];
      auto slice = m_slice_buf.begin();
      for (size_t iz = 0; iz < z; iz++) {
        auto beg = m_data_buf.begin() + iz * plane_size_xy + offset;
        std::copy(beg, beg + x, slice);
        slice += x;
      }
      m_idwt2d_one_level(m_slice_buf.begin(), {x, z});
      slice = m_slice_buf.begin();
      for (size_t iz = 0; iz < z; iz++) {
        auto beg = m_data_buf.begin() + iz * plane_size_xy + offset;
        std::copy(slice, slice + x, beg);
        slice += x;
      }
    }

    // XY planes
    for (size_t iz = 0; iz < z; iz++) {
      auto offset = iz * plane_size_xy;
      m_idwt2d_one_level(m_data_buf.begin() + offset, {x, y});
    }
  }
}

void sperr::CDF97::m_idwt3d_wavelet_packet()
{
  const size_t plane_size_xy = m_dims[0] * m_dims[1];

  // First, inverse transform each plane
  auto num_xforms_xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  for (size_t i = 0; i < m_dims[2]; i++) {
    const size_t offset = plane_size_xy * i;
    m_idwt2d(m_data_buf.begin() + offset, {m_dims[0], m_dims[1]}, num_xforms_xy);
  }

  /*
   * Second, inverse transform along the Z dimension
   */

  const auto num_xforms_z = sperr::num_of_xforms(m_dims[2]);

  for (size_t lev = num_xforms_z; lev > 0; lev--) {
    auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev - 1);

    // Process one XZ slice at a time
    auto slice_begin = m_slice_buf.begin();
    for (size_t y = 0; y < m_dims[1]; y++) {
      // Copy one XZ slice
      auto dst = slice_begin;
      for (size_t zz = 0; zz < z; zz++) {
        auto src = m_data_buf.begin() + zz * plane_size_xy + y * m_dims[0];
        std::copy(src, src + m_dims[0], dst);
        dst += m_dims[0];
      }

      // 2D IDWT on that slice (X-Z)
      m_idwt2d(slice_begin, {m_dims[0], z}, 1);

      // Copy back
      dst = slice_begin;
      for (size_t zz = 0; zz < z; zz++) {
        auto src = m_data_buf.begin() + zz * plane_size_xy + y * m_dims[0];
        std::copy(dst, dst + m_dims[0], src);
        dst += m_dims[0];
      }
    }
  }
}

void sperr::CDF97::m_idwt3d_one_level(itd_type vol, std::array<size_t, 3> len_xyz)
{
  const size_t plane_size_xy = m_dims[0] * m_dims[1];

  const auto max_len = std::max({len_xyz[0], len_xyz[1], len_xyz[2]});
  if (m_qcc_buf.size() < max_len)
    m_qcc_buf.resize(max_len);

  // 1) Inverse along Z (stride-aware) for every (x,y) column.
  const ptrdiff_t stride_z = static_cast<ptrdiff_t>(plane_size_xy);
  for (size_t y = 0; y < len_xyz[1]; ++y) {
    for (size_t x = 0; x < len_xyz[0]; ++x) {
      double* col0 = vol + y * m_dims[0] + x;  // (x,y,0)
      m_idwt1d_one_level_strided(col0, len_xyz[2], stride_z);
    }
  }

  // 2) Inverse 2D transform on all XY planes.
  for (size_t z = 0; z < len_xyz[2]; ++z) {
    const size_t offset = plane_size_xy * z;
    m_idwt2d_one_level(vol + offset, {len_xyz[0], len_xyz[1]});
  }
}

void sperr::CDF97::m_gather_even(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 0);  // This function specifically for even length input
  auto half = len / 2;

  // Low-pass: even samples
  for (size_t i = 0; i < half; i++)
    dest[i] = begin[i * 2];

  // High-pass: odd samples
  for (size_t i = 0; i < half; i++)
    dest[half + i] = begin[i * 2 + 1];
}

void sperr::CDF97::m_gather_odd(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 1);  // This function specifically for odd length input
  size_t low_count = len / 2 + 1, high_count = len / 2;

  for (size_t i = 0; i < low_count; i++)
    dest[i] = begin[i * 2];

  for (size_t i = 0; i < high_count; i++)
    dest[low_count + i] = begin[i * 2 + 1];
}

void sperr::CDF97::m_scatter_even(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 0);  // This function specifically for even length input
  auto half = len / 2;

  for (size_t i = 0; i < half; i++)
    dest[i * 2] = begin[i];

  for (size_t i = 0; i < half; i++)
    dest[i * 2 + 1] = begin[half + i];
}

void sperr::CDF97::m_scatter_odd(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 1);  // This function specifically for odd length input
  size_t low_count = len / 2 + 1, high_count = len / 2;

  for (size_t i = 0; i < low_count; i++)
    dest[i * 2] = begin[i];

  for (size_t i = 0; i < high_count; i++)
    dest[i * 2 + 1] = begin[low_count + i];
}

auto sperr::CDF97::m_sub_slice(std::array<size_t, 2> subdims) const -> vecd_type
{
  assert(subdims[0] <= m_dims[0] && subdims[1] <= m_dims[1]);

  auto ret = vecd_type(subdims[0] * subdims[1]);
  auto dst = ret.begin();
  for (size_t y = 0; y < subdims[1]; y++) {
    auto beg = m_data_buf.begin() + y * m_dims[0];
    std::copy(beg, beg + subdims[0], dst);
    dst += subdims[0];
  }
  return ret;
}

void sperr::CDF97::m_sub_slice(std::array<size_t, 2> subdims, itd_type dest) const
{
  assert(subdims[0] <= m_dims[0] && subdims[1] <= m_dims[1]);

  for (size_t y = 0; y < subdims[1]; y++) {
    auto beg = m_data_buf.begin() + y * m_dims[0];
    std::copy(beg, beg + subdims[0], dest);
    dest += subdims[0];
  }
}

auto sperr::CDF97::m_sub_volume(std::array<size_t, 3> subdims) const -> vecd_type
{
  assert(subdims[0] <= m_dims[0] && subdims[1] <= m_dims[1] && subdims[2] <= m_dims[2]);

  auto ret = vecd_type(subdims[0] * subdims[1] * subdims[2]);
  auto dst = ret.begin();
  const auto plane_size_xy = m_dims[0] * m_dims[1];

  for (size_t z = 0; z < subdims[2]; z++)
    for (size_t y = 0; y < subdims[1]; y++) {
      auto beg = m_data_buf.begin() + z * plane_size_xy + y * m_dims[0];
      std::copy(beg, beg + subdims[0], dst);
      dst += subdims[0];
    }

  return ret;
}

void sperr::CDF97::m_sub_volume(std::array<size_t, 3> subdims, itd_type dest) const
{
  assert(subdims[0] <= m_dims[0] && subdims[1] <= m_dims[1] && subdims[2] <= m_dims[2]);

  const auto plane_size_xy = m_dims[0] * m_dims[1];

  for (size_t z = 0; z < subdims[2]; z++)
    for (size_t y = 0; y < subdims[1]; y++) {
      auto beg = m_data_buf.begin() + z * plane_size_xy + y * m_dims[0];
      std::copy(beg, beg + subdims[0], dest);
      dest += subdims[0];
    }
}

// ===== Original (contiguous) lifting kernels =====

void sperr::CDF97::QccWAVCDF97AnalysisSymmetricEvenEven(double* signal, size_t signal_length)
{
  for (size_t i = 1; i < signal_length - 2; i += 2)
    signal[i] += ALPHA * (signal[i - 1] + signal[i + 1]);

  signal[signal_length - 1] += 2.0 * ALPHA * signal[signal_length - 2];

  signal[0] += 2.0 * BETA * signal[1];

  for (size_t i = 2; i < signal_length; i += 2)
    signal[i] += BETA * (signal[i + 1] + signal[i - 1]);

  for (size_t i = 1; i < signal_length - 2; i += 2)
    signal[i] += GAMMA * (signal[i - 1] + signal[i + 1]);

  signal[signal_length - 1] += 2.0 * GAMMA * signal[signal_length - 2];

  for (size_t i = 2; i < signal_length; i += 2)
    signal[i] += DELTA * (signal[i + 1] + signal[i - 1]);

  for (size_t i = 0; i < signal_length; i += 2)
    signal[i] *= EPSILON;

  for (size_t i = 1; i < signal_length; i += 2)
    signal[i] *= -INV_EPSILON;
}

void sperr::CDF97::QccWAVCDF97SynthesisSymmetricEvenEven(double* signal, size_t signal_length)
{
  for (size_t i = 1; i < signal_length; i += 2)
    signal[i] *= -EPSILON;

  for (size_t i = 0; i < signal_length; i += 2)
    signal[i] *= INV_EPSILON;

  for (size_t i = 2; i < signal_length; i += 2)
    signal[i] -= DELTA * (signal[i + 1] + signal[i - 1]);

  for (size_t i = 1; i < signal_length - 2; i += 2)
    signal[i] -= GAMMA * (signal[i - 1] + signal[i + 1]);

  signal[signal_length - 1] -= 2.0 * GAMMA * signal[signal_length - 2];

  signal[0] -= 2.0 * BETA * signal[1];

  for (size_t i = 2; i < signal_length; i += 2)
    signal[i] -= BETA * (signal[i + 1] + signal[i - 1]);

  for (size_t i = 1; i < signal_length - 2; i += 2)
    signal[i] -= ALPHA * (signal[i - 1] + signal[i + 1]);

  signal[signal_length - 1] -= 2.0 * ALPHA * signal[signal_length - 2];
}

void sperr::CDF97::QccWAVCDF97SynthesisSymmetricOddEven(double* signal, size_t signal_length)
{
  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] *= -EPSILON;

  signal[0] = signal[0] * INV_EPSILON - 2.0 * DELTA * signal[1];

  for (size_t i = 2; i < signal_length - 2; i += 2)
    signal[i] = signal[i] * INV_EPSILON -
                DELTA * (signal[i + 1] + signal[i - 1]);

  signal[signal_length - 1] =
      signal[signal_length - 1] * INV_EPSILON - 2.0 * DELTA * signal[signal_length - 2];

  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] -= GAMMA * (signal[i - 1] + signal[i + 1]);

  signal[0] -= 2.0 * BETA * signal[1];

  for (size_t i = 2; i < signal_length - 2; i += 2)
    signal[i] -= BETA * (signal[i + 1] + signal[i - 1]);

  signal[signal_length - 1] -= 2.0 * BETA * signal[signal_length - 2];

  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] -= ALPHA * (signal[i - 1] + signal[i + 1]);
}

void sperr::CDF97::QccWAVCDF97AnalysisSymmetricOddEven(double* signal, size_t signal_length)
{
  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] += ALPHA * (signal[i - 1] + signal[i + 1]);

  signal[0] += 2.0 * BETA * signal[1];

  for (size_t i = 2; i < signal_length - 2; i += 2)
    signal[i] += BETA * (signal[i + 1] + signal[i - 1]);

  signal[signal_length - 1] += 2.0 * BETA * signal[signal_length - 2];

  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] += GAMMA * (signal[i - 1] + signal[i + 1]);

  signal[0] = EPSILON * (signal[0] + 2.0 * DELTA * signal[1]);

  for (size_t i = 2; i < signal_length - 2; i += 2)
    signal[i] = EPSILON * (signal[i] + DELTA * (signal[i + 1] + signal[i - 1]));

  signal[signal_length - 1] =
      EPSILON * (signal[signal_length - 1] + 2.0 * DELTA * signal[signal_length - 2]);

  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] *= -INV_EPSILON;
}

// ===== New stride-aware lifting kernels =====

void sperr::CDF97::QccWAVCDF97AnalysisSymmetricEvenEvenStrided(
    itd_type base, size_t n, ptrdiff_t s)
{
  for (size_t i = 1; i < n - 2; i += 2)
    base[i * s] += ALPHA * (base[(i - 1) * s] + base[(i + 1) * s]);

  base[(n - 1) * s] += 2.0 * ALPHA * base[(n - 2) * s];

  base[0] += 2.0 * BETA * base[1 * s];

  for (size_t i = 2; i < n; i += 2)
    base[i * s] += BETA * (base[(i + 1) * s] + base[(i - 1) * s]);

  for (size_t i = 1; i < n - 2; i += 2)
    base[i * s] += GAMMA * (base[(i - 1) * s] + base[(i + 1) * s]);

  base[(n - 1) * s] += 2.0 * GAMMA * base[(n - 2) * s];

  for (size_t i = 2; i < n; i += 2)
    base[i * s] += DELTA * (base[(i + 1) * s] + base[(i - 1) * s]);

  for (size_t i = 0; i < n; i += 2)
    base[i * s] *= EPSILON;

  for (size_t i = 1; i < n; i += 2)
    base[i * s] *= -INV_EPSILON;
}

void sperr::CDF97::QccWAVCDF97AnalysisSymmetricOddEvenStrided(
    itd_type base, size_t n, ptrdiff_t s)
{
  for (size_t i = 1; i < n - 1; i += 2)
    base[i * s] += ALPHA * (base[(i - 1) * s] + base[(i + 1) * s]);

  base[0] += 2.0 * BETA * base[1 * s];

  for (size_t i = 2; i < n - 2; i += 2)
    base[i * s] += BETA * (base[(i + 1) * s] + base[(i - 1) * s]);

  base[(n - 1) * s] += 2.0 * BETA * base[(n - 2) * s];

  for (size_t i = 1; i < n - 1; i += 2)
    base[i * s] += GAMMA * (base[(i - 1) * s] + base[(i + 1) * s]);

  base[0] = EPSILON * (base[0] + 2.0 * DELTA * base[1 * s]);

  for (size_t i = 2; i < n - 2; i += 2)
    base[i * s] =
        EPSILON * (base[i * s] + DELTA * (base[(i + 1) * s] + base[(i - 1) * s]));

  base[(n - 1) * s] =
      EPSILON * (base[(n - 1) * s] + 2.0 * DELTA * base[(n - 2) * s]);

  for (size_t i = 1; i < n - 1; i += 2)
    base[i * s] *= -INV_EPSILON;
}

void sperr::CDF97::QccWAVCDF97SynthesisSymmetricEvenEvenStrided(
    itd_type base, size_t n, ptrdiff_t s)
{
  for (size_t i = 1; i < n; i += 2)
    base[i * s] *= -EPSILON;

  for (size_t i = 0; i < n; i += 2)
    base[i * s] *= INV_EPSILON;

  for (size_t i = 2; i < n; i += 2)
    base[i * s] -= DELTA * (base[(i + 1) * s] + base[(i - 1) * s]);

  for (size_t i = 1; i < n - 2; i += 2)
    base[i * s] -= GAMMA * (base[(i - 1) * s] + base[(i + 1) * s]);

  base[(n - 1) * s] -= 2.0 * GAMMA * base[(n - 2) * s];

  base[0] -= 2.0 * BETA * base[1 * s];

  for (size_t i = 2; i < n; i += 2)
    base[i * s] -= BETA * (base[(i + 1) * s] + base[(i - 1) * s]);

  for (size_t i = 1; i < n - 2; i += 2)
    base[i * s] -= ALPHA * (base[(i - 1) * s] + base[(i + 1) * s]);

  base[(n - 1) * s] -= 2.0 * ALPHA * base[(n - 2) * s];
}

void sperr::CDF97::QccWAVCDF97SynthesisSymmetricOddEvenStrided(
    itd_type base, size_t n, ptrdiff_t s)
{
  for (size_t i = 1; i < n - 1; i += 2)
    base[i * s] *= -EPSILON;

  base[0] = base[0] * INV_EPSILON - 2.0 * DELTA * base[1 * s];

  for (size_t i = 2; i < n - 2; i += 2)
    base[i * s] =
        base[i * s] * INV_EPSILON -
        DELTA * (base[(i + 1) * s] + base[(i - 1) * s]);

  base[(n - 1) * s] =
      base[(n - 1) * s] * INV_EPSILON - 2.0 * DELTA * base[(n - 2) * s];

  for (size_t i = 1; i < n - 1; i += 2)
    base[i * s] -= GAMMA * (base[(i - 1) * s] + base[(i + 1) * s]);

  base[0] -= 2.0 * BETA * base[1 * s];

  for (size_t i = 2; i < n - 2; i += 2)
    base[i * s] -= BETA * (base[(i + 1) * s] + base[(i - 1) * s]);

  base[(n - 1) * s] -= 2.0 * BETA * base[(n - 2) * s];

  for (size_t i = 1; i < n - 1; i += 2)
    base[i * s] -= ALPHA * (base[(i - 1) * s] + base[(i + 1) * s]);
}
