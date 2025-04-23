#include "CDF97.h"
#include "Sym13Coeffs.h"
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

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[0] * dims[2]), dims[1] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(std::max(m_slice_buf.size() * 2, max_slice));

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

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[0] * dims[2]), dims[1] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(std::max(m_slice_buf.size() * 2, max_slice));

  return RTNType::Good;
}

auto sperr::CDF97::view_data() const -> const vecd_type&
{
  return m_data_buf;
}

auto sperr::CDF97::release_data() -> vecd_type&&
{
  return std::move(m_data_buf);
}

auto sperr::CDF97::get_dims() const -> std::array<size_t, 3>
{
  return m_dims;
}

void sperr::CDF97::dwt1d()
{
  auto num_xforms = sperr::num_of_xforms(m_dims[0]);
  m_dwt1d(m_data_buf.begin(), m_data_buf.size(), num_xforms);
}

void sperr::CDF97::idwt1d()
{
  auto num_xforms = sperr::num_of_xforms(m_dims[0]);
  m_idwt1d(m_data_buf.begin(), m_data_buf.size(), num_xforms);
}

void sperr::CDF97::dwt2d()
{
  auto xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  m_dwt2d(m_data_buf.begin(), {m_dims[0], m_dims[1]}, xy);
}

void sperr::CDF97::idwt2d()
{
  auto xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  m_idwt2d(m_data_buf.begin(), {m_dims[0], m_dims[1]}, xy);
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
      m_dwt1d(m_slice_buf.begin() + x * m_dims[2], m_dims[2], num_xforms_z);

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
    m_dwt2d(m_data_buf.begin() + offset, {m_dims[0], m_dims[1]}, num_xforms_xy);
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
    m_idwt2d(m_data_buf.begin() + offset, {m_dims[0], m_dims[1]}, num_xforms_xy);
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
      m_idwt1d(m_slice_buf.begin() + x * m_dims[2], m_dims[2], num_xforms_z);

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
    m_dwt3d_one_level(m_data_buf.begin(), {x, y, z});
  }
}

void sperr::CDF97::m_idwt3d_dyadic(size_t num_xforms)
{
  for (size_t lev = num_xforms; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev - 1);
    auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev - 1);
    auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev - 1);
    m_idwt3d_one_level(m_data_buf.begin(), {x, y, z});
  }
}

//
// Private Methods
//
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
  std::copy(array, array + array_len, m_qcc_buf.begin());
  if (m_wavelet_type == "sym13"){
    this->Sym13Analysis(m_qcc_buf.data(), array_len);
    if (array_len % 2 == 0) {
      m_gather_even(m_qcc_buf.cbegin(), m_qcc_buf.cbegin() + array_len, array);
    }
    else {
      m_gather_odd(m_qcc_buf.cbegin(), m_qcc_buf.cbegin() + array_len, array);
    }
    
  }
  else{
    if (array_len % 2 == 0) {
      this->QccWAVCDF97AnalysisSymmetricEvenEven(m_qcc_buf.data(), array_len);
      m_gather_even(m_qcc_buf.cbegin(), m_qcc_buf.cbegin() + array_len, array);
    }
    else {
      this->QccWAVCDF97AnalysisSymmetricOddEven(m_qcc_buf.data(), array_len);
      m_gather_odd(m_qcc_buf.cbegin(), m_qcc_buf.cbegin() + array_len, array);
    }
  }
}

void sperr::CDF97::m_idwt1d_one_level(itd_type array, size_t array_len)
{
  if (m_wavelet_type == "sym13"){
    if (array_len % 2 == 0) {
      m_scatter_even(array, array + array_len, m_qcc_buf.begin());
    }
    else {
      m_scatter_odd(array, array + array_len, m_qcc_buf.begin());
    }
    this->Sym13Synthesis(m_qcc_buf.data(), array_len);
  }
  else{
    if (array_len % 2 == 0) {
      m_scatter_even(array, array + array_len, m_qcc_buf.begin());
      this->QccWAVCDF97SynthesisSymmetricEvenEven(m_qcc_buf.data(), array_len);
    }
    else {
      m_scatter_odd(array, array + array_len, m_qcc_buf.begin());
      this->QccWAVCDF97SynthesisSymmetricOddEven(m_qcc_buf.data(), array_len);
    }
  }
  std::copy(m_qcc_buf.cbegin(), m_qcc_buf.cbegin() + array_len, array);
}

void sperr::CDF97::m_dwt2d_one_level(itd_type plane, std::array<size_t, 2> len_xy)
{
  // Note: here we call low-level functions (Qcc*()) instead of
  // m_dwt1d_one_level() because we want to have only one even/odd test at the outer loop.

  const auto max_len = std::max(len_xy[0], len_xy[1]);
  const auto beg = m_qcc_buf.begin();
  const auto beg2 = beg + max_len;

  // First, perform DWT along X for every row

   // Second, perform DWT along Y for every column
  // Note, I've tested that up to 1024^2 planes it is actually slightly slower
  // to transpose the plane and then perform the transforms. This was consistent
  // on both a MacBook and a RaspberryPi 3. Note2, I've tested transpose again
  // on an X86 linux machine using gcc, clang, and pgi. Again the difference is
  // either indistinguishable, or the current implementation has a slight edge.


  if (m_wavelet_type == "sym13"){
    if (len_xy[0] % 2 == 0) {
      for (size_t i = 0; i < len_xy[1]; i++) {
        auto pos = plane + i * m_dims[0];
        std::copy(pos, pos + len_xy[0], beg);
        this->Sym13Analysis(m_qcc_buf.data(), len_xy[0]);
        m_gather_even(beg, beg + len_xy[0], pos);
      }
    }
    else  // Odd length
    {
      for (size_t i = 0; i < len_xy[1]; i++) {
        auto pos = plane + i * m_dims[0];
        std::copy(pos, pos + len_xy[0], beg);
        this->Sym13Analysis(m_qcc_buf.data(), len_xy[0]);
        m_gather_odd(beg, beg + len_xy[0], pos);
      }
    }
    if (len_xy[1] % 2 == 0) {
      for (size_t x = 0; x < len_xy[0]; x++) {
        for (size_t y = 0; y < len_xy[1]; y++)
          m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
        this->Sym13Analysis(m_qcc_buf.data(), len_xy[1]);
        m_gather_even(beg, beg + len_xy[1], beg2);
        for (size_t y = 0; y < len_xy[1]; y++)
          *(plane + y * m_dims[0] + x) = *(beg2 + y);
      }
    }
    else  // Odd length
    {
      for (size_t x = 0; x < len_xy[0]; x++) {
        for (size_t y = 0; y < len_xy[1]; y++)
          m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
        this->Sym13Analysis(m_qcc_buf.data(), len_xy[1]);
        m_gather_odd(beg, beg + len_xy[1], beg2);
        for (size_t y = 0; y < len_xy[1]; y++)
          *(plane + y * m_dims[0] + x) = *(beg2 + y);
      }
    }
  }
  else{
    if (len_xy[0] % 2 == 0) {
      for (size_t i = 0; i < len_xy[1]; i++) {
        auto pos = plane + i * m_dims[0];
        std::copy(pos, pos + len_xy[0], beg);
        this->QccWAVCDF97AnalysisSymmetricEvenEven(m_qcc_buf.data(), len_xy[0]);
        m_gather_even(beg, beg + len_xy[0], pos);
      }
    }
    else  // Odd length
    {
      for (size_t i = 0; i < len_xy[1]; i++) {
        auto pos = plane + i * m_dims[0];
        std::copy(pos, pos + len_xy[0], beg);
        this->QccWAVCDF97AnalysisSymmetricOddEven(m_qcc_buf.data(), len_xy[0]);
        m_gather_odd(beg, beg + len_xy[0], pos);
      }
    }
    if (len_xy[1] % 2 == 0) {
      for (size_t x = 0; x < len_xy[0]; x++) {
        for (size_t y = 0; y < len_xy[1]; y++)
          m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
        this->QccWAVCDF97AnalysisSymmetricEvenEven(m_qcc_buf.data(), len_xy[1]);
        m_gather_even(beg, beg + len_xy[1], beg2);
        for (size_t y = 0; y < len_xy[1]; y++)
          *(plane + y * m_dims[0] + x) = *(beg2 + y);
      }
    }
    else  // Odd length
    {
      for (size_t x = 0; x < len_xy[0]; x++) {
        for (size_t y = 0; y < len_xy[1]; y++)
          m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
        this->QccWAVCDF97AnalysisSymmetricOddEven(m_qcc_buf.data(), len_xy[1]);
        m_gather_odd(beg, beg + len_xy[1], beg2);
        for (size_t y = 0; y < len_xy[1]; y++)
          *(plane + y * m_dims[0] + x) = *(beg2 + y);
      }
    }
  }

 

  
}

void sperr::CDF97::m_idwt2d_one_level(itd_type plane, std::array<size_t, 2> len_xy)
{
  const auto max_len = std::max(len_xy[0], len_xy[1]);
  const auto beg = m_qcc_buf.begin();  // First half of the buffer
  const auto beg2 = beg + max_len;     // Second half of the buffer

  // First, perform IDWT along Y for every column
  if (m_wavelet_type == "sym13"){
    if (len_xy[1] % 2 == 0) {
      for (size_t x = 0; x < len_xy[0]; x++) {
        for (size_t y = 0; y < len_xy[1]; y++)
          m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
        m_scatter_even(beg, beg + len_xy[1], beg2);
        this->Sym13Synthesis(m_qcc_buf.data() + max_len, len_xy[1]);
        for (size_t y = 0; y < len_xy[1]; y++)
          *(plane + y * m_dims[0] + x) = *(beg2 + y);
      }
    }
    else  // Odd length
    {
      for (size_t x = 0; x < len_xy[0]; x++) {
        for (size_t y = 0; y < len_xy[1]; y++)
          m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
        m_scatter_odd(beg, beg + len_xy[1], beg2);
        this->Sym13Synthesis(m_qcc_buf.data() + max_len, len_xy[1]);
        for (size_t y = 0; y < len_xy[1]; y++)
          *(plane + y * m_dims[0] + x) = *(beg2 + y);
      }
    }

    // Second, perform IDWT along X for every row
    if (len_xy[0] % 2 == 0) {
      for (size_t i = 0; i < len_xy[1]; i++) {
        auto pos = plane + i * m_dims[0];
        m_scatter_even(pos, pos + len_xy[0], beg);
        this->Sym13Synthesis(m_qcc_buf.data(), len_xy[0]);
        std::copy(beg, beg + len_xy[0], pos);
      }
    }
    else  // Odd length
    {
      for (size_t i = 0; i < len_xy[1]; i++) {
        auto pos = plane + i * m_dims[0];
        m_scatter_odd(pos, pos + len_xy[0], beg);
        this->Sym13Synthesis(m_qcc_buf.data(), len_xy[0]);
        std::copy(beg, beg + len_xy[0], pos);
      }
    }
  }
  else{
    if (len_xy[1] % 2 == 0) {
      for (size_t x = 0; x < len_xy[0]; x++) {
        for (size_t y = 0; y < len_xy[1]; y++)
          m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
        m_scatter_even(beg, beg + len_xy[1], beg2);
        this->QccWAVCDF97SynthesisSymmetricEvenEven(m_qcc_buf.data() + max_len, len_xy[1]);
        for (size_t y = 0; y < len_xy[1]; y++)
          *(plane + y * m_dims[0] + x) = *(beg2 + y);
      }
    }
    else  // Odd length
    {
      for (size_t x = 0; x < len_xy[0]; x++) {
        for (size_t y = 0; y < len_xy[1]; y++)
          m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
        m_scatter_odd(beg, beg + len_xy[1], beg2);
        this->QccWAVCDF97SynthesisSymmetricOddEven(m_qcc_buf.data() + max_len, len_xy[1]);
        for (size_t y = 0; y < len_xy[1]; y++)
          *(plane + y * m_dims[0] + x) = *(beg2 + y);
      }
    }

    // Second, perform IDWT along X for every row
    if (len_xy[0] % 2 == 0) {
      for (size_t i = 0; i < len_xy[1]; i++) {
        auto pos = plane + i * m_dims[0];
        m_scatter_even(pos, pos + len_xy[0], beg);
        this->QccWAVCDF97SynthesisSymmetricEvenEven(m_qcc_buf.data(), len_xy[0]);
        std::copy(beg, beg + len_xy[0], pos);
      }
    }
    else  // Odd length
    {
      for (size_t i = 0; i < len_xy[1]; i++) {
        auto pos = plane + i * m_dims[0];
        m_scatter_odd(pos, pos + len_xy[0], beg);
        this->QccWAVCDF97SynthesisSymmetricOddEven(m_qcc_buf.data(), len_xy[0]);
        std::copy(beg, beg + len_xy[0], pos);
      }
    }
  }
}

void sperr::CDF97::m_dwt3d_one_level(itd_type vol, std::array<size_t, 3> len_xyz)
{
  // First, do one level of transform on all XY planes.
  const auto plane_size_xy = m_dims[0] * m_dims[1];
  for (size_t z = 0; z < len_xyz[2]; z++) {
    const size_t offset = plane_size_xy * z;
    m_dwt2d_one_level(vol + offset, {len_xyz[0], len_xyz[1]});
  }

  const auto beg = m_qcc_buf.begin();  // First half of the buffer
  const auto beg2 = beg + len_xyz[2];  // Second half of the buffer

  // Second, do one level of transform on all Z columns.  Strategy:
  // 1) extract a Z column to buffer space `m_qcc_buf`
  // 2) use appropriate even/odd Qcc*** function to transform it
  // 3) gather coefficients from `m_qcc_buf` to the second half of `m_qcc_buf`
  // 4) put the Z column back to their locations as a Z column.
  if(m_wavelet_type == "sym13"){
    if (len_xyz[2] % 2 == 0) {  // Even length
      for (size_t y = 0; y < len_xyz[1]; y++) {
        for (size_t x = 0; x < len_xyz[0]; x++) {
          const size_t xy_offset = y * m_dims[0] + x;
          // Step 1
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
          // Step 2
          this->Sym13Analysis(m_qcc_buf.data(), len_xyz[2]);
          // Step 3
          m_gather_even(beg, beg2, beg2);
          // Step 4
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
        }
      }
    }
    else {  // Odd length
      for (size_t y = 0; y < len_xyz[1]; y++) {
        for (size_t x = 0; x < len_xyz[0]; x++) {
          const size_t xy_offset = y * m_dims[0] + x;
          // Step 1
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
          // Step 2
          this->Sym13Analysis(m_qcc_buf.data(), len_xyz[2]);
          // Step 3
          m_gather_odd(beg, beg2, beg2);
          // Step 4
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
        }
      }
    }
  }
  else{
    if (len_xyz[2] % 2 == 0) {  // Even length
      for (size_t y = 0; y < len_xyz[1]; y++) {
        for (size_t x = 0; x < len_xyz[0]; x++) {
          const size_t xy_offset = y * m_dims[0] + x;
          // Step 1
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
          // Step 2
          this->QccWAVCDF97AnalysisSymmetricEvenEven(m_qcc_buf.data(), len_xyz[2]);
          // Step 3
          m_gather_even(beg, beg2, beg2);
          // Step 4
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
        }
      }
    }
    else {  // Odd length
      for (size_t y = 0; y < len_xyz[1]; y++) {
        for (size_t x = 0; x < len_xyz[0]; x++) {
          const size_t xy_offset = y * m_dims[0] + x;
          // Step 1
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
          // Step 2
          this->QccWAVCDF97AnalysisSymmetricOddEven(m_qcc_buf.data(), len_xyz[2]);
          // Step 3
          m_gather_odd(beg, beg2, beg2);
          // Step 4
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
        }
      }
    }
  }
}

void sperr::CDF97::m_idwt3d_one_level(itd_type vol, std::array<size_t, 3> len_xyz)
{
  const auto plane_size_xy = m_dims[0] * m_dims[1];
  const auto beg = m_qcc_buf.begin();  // First half of the buffer
  const auto beg2 = beg + len_xyz[2];  // Second half of the buffer

  // First, do one level of inverse transform on all Z columns.  Strategy:
  // 1) extract a Z column to buffer space `m_qcc_buf`
  // 2) scatter coefficients from `m_qcc_buf` to the second half of `m_qcc_buf`
  // 3) use appropriate even/odd Qcc*** function to transform it
  // 4) put the Z column back to their locations as a Z column.
  if (m_wavelet_type == "sym13"){
    if (len_xyz[2] % 2 == 0) {
      for (size_t y = 0; y < len_xyz[1]; y++) {
        for (size_t x = 0; x < len_xyz[0]; x++) {
          const size_t xy_offset = y * m_dims[0] + x;
          // Step 1
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
          // Step 2
          m_scatter_even(beg, beg2, beg2);
          // Step 3
          this->Sym13Synthesis(m_qcc_buf.data() + len_xyz[2], len_xyz[2]);
          // Step 4
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
        }
      }
    }
    else {
      for (size_t y = 0; y < len_xyz[1]; y++) {
        for (size_t x = 0; x < len_xyz[0]; x++) {
          const size_t xy_offset = y * m_dims[0] + x;
          // Step 1
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
          // Step 2
          m_scatter_odd(beg, beg2, beg2);
          // Step 3
          this->Sym13Synthesis(m_qcc_buf.data() + len_xyz[2], len_xyz[2]);
          // Step 4
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
        }
      }
    }
  }
  else{
    if (len_xyz[2] % 2 == 0) {
      for (size_t y = 0; y < len_xyz[1]; y++) {
        for (size_t x = 0; x < len_xyz[0]; x++) {
          const size_t xy_offset = y * m_dims[0] + x;
          // Step 1
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
          // Step 2
          m_scatter_even(beg, beg2, beg2);
          // Step 3
          this->QccWAVCDF97SynthesisSymmetricEvenEven(m_qcc_buf.data() + len_xyz[2], len_xyz[2]);
          // Step 4
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
        }
      }
    }
    else {
      for (size_t y = 0; y < len_xyz[1]; y++) {
        for (size_t x = 0; x < len_xyz[0]; x++) {
          const size_t xy_offset = y * m_dims[0] + x;
          // Step 1
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
          // Step 2
          m_scatter_odd(beg, beg2, beg2);
          // Step 3
          this->QccWAVCDF97SynthesisSymmetricOddEven(m_qcc_buf.data() + len_xyz[2], len_xyz[2]);
          // Step 4
          for (size_t z = 0; z < len_xyz[2]; z++)
            m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
        }
      }
    }

  }
  
  // Second, do one level of inverse transform on all XY planes.
  for (size_t z = 0; z < len_xyz[2]; z++) {
    const size_t offset = plane_size_xy * z;
    m_idwt2d_one_level(vol + offset, {len_xyz[0], len_xyz[1]});
  }
}

void sperr::CDF97::m_gather_even(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 0);  // This function specifically for even length input
  size_t low_count = len / 2, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *dest = *(begin + i * 2);
    ++dest;
  }
  for (size_t i = 0; i < high_count; i++) {
    *dest = *(begin + i * 2 + 1);
    ++dest;
  }
}

void sperr::CDF97::m_gather_odd(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 1);  // This function specifically for odd length input
  size_t low_count = len / 2 + 1, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *dest = *(begin + i * 2);
    ++dest;
  }
  for (size_t i = 0; i < high_count; i++) {
    *dest = *(begin + i * 2 + 1);
    ++dest;
  }
}

void sperr::CDF97::m_scatter_even(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 0);  // This function specifically for even length input
  size_t low_count = len / 2, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *(dest + i * 2) = *begin;
    ++begin;
  }
  for (size_t i = 0; i < high_count; i++) {
    *(dest + i * 2 + 1) = *begin;
    ++begin;
  }
}

void sperr::CDF97::m_scatter_odd(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 1);  // This function specifically for odd length input
  size_t low_count = len / 2 + 1, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *(dest + i * 2) = *begin;
    ++begin;
  }
  for (size_t i = 0; i < high_count; i++) {
    *(dest + i * 2 + 1) = *begin;
    ++begin;
  }
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

void sperr::CDF97::m_sub_volume(dims_type subdims, itd_type dst) const
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

  signal[0] = EPSILON * (signal[0] + 2.0 * DELTA * signal[1]);

  for (size_t i = 2; i < signal_length; i += 2)
    signal[i] = EPSILON * (signal[i] + DELTA * (signal[i + 1] + signal[i - 1]));

  for (size_t i = 1; i < signal_length; i += 2)
    signal[i] *= -INV_EPSILON;
}

void sperr::CDF97::QccWAVCDF97SynthesisSymmetricEvenEven(double* signal, size_t signal_length)
{
  for (size_t i = 1; i < signal_length; i += 2)
    signal[i] *= (-EPSILON);

  signal[0] = signal[0] * INV_EPSILON - 2.0 * DELTA * signal[1];

  for (size_t i = 2; i < signal_length; i += 2)
    signal[i] = signal[i] * INV_EPSILON - DELTA * (signal[i + 1] + signal[i - 1]);

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
    signal[i] *= (-EPSILON);

  signal[0] = signal[0] * INV_EPSILON - 2.0 * DELTA * signal[1];

  for (size_t i = 2; i < signal_length - 2; i += 2)
    signal[i] = signal[i] * INV_EPSILON - DELTA * (signal[i + 1] + signal[i - 1]);

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
    signal[i] *= (-INV_EPSILON);
}

void CDF97::Sym13Analysis(double* signal, size_t n) {
  //using namespace sym13;
  const int L   = int(sym13::kernel_length);   // 26
  const int pad = L - 1;         // 25

  // 1) 延拓到 ext 长度 n + 2*pad
  std::vector<double> ext(n + 2*pad);
  if (m_padding_mode == "periodic") {
    // 周期延拓
    for (int i = 0; i < pad; ++i) {
      int idx = (i - pad) % int(n);
      if (idx < 0) idx += int(n);
      ext[i] = signal[idx];
    }
    std::copy(signal, signal + n, ext.begin() + pad);
    for (int i = 0; i < pad; ++i) {
      ext[n + pad + i] = signal[i % n];
    }
  } else {
    // 对称延拓
    for (int i = 0; i < pad; ++i) {
      ext[i] = signal[pad - 1 - i];
    }
    std::copy(signal, signal + n, ext.begin() + pad);
    for (int i = 0; i < pad; ++i) {
      ext[n + pad + i] = signal[n - 1 - i];
    }
  }

  // 2) 计算低频/高频输出长度
  size_t nA = (n + 1) / 2;  // ceil(n/2)
  size_t nD = n / 2;        // floor(n/2)

  // 3) 卷积 + 下采样
  std::vector<double> out(n);
  // 低频部分
  for (size_t i = 0; i < nA; ++i) {
    double a = 0;
    size_t start = 2 * i;
    for (int k = 0; k < L; ++k) {
      a += sym13::dec_lo[k] * ext[start + k];
    }
    out[i] = a;
  }
  // 高频部分
  for (size_t i = 0; i < nD; ++i) {
    double d = 0;
    size_t start = 2 * i;
    for (int k = 0; k < L; ++k) {
      d += sym13::dec_hi[k] * ext[start + k];
    }
    out[nA + i] = d;
  }

  // 4) 写回
  std::copy(out.begin(), out.end(), signal);
}

void CDF97::Sym13Synthesis(double* signal, size_t n, bool periodic) {
  //using namespace sym13;
  const int L   = int(sym13::kernel_length);
  const int pad = L - 1;

  // 1) 从 signal[0..nA-1]、[nA..n-1] 拆出低频 cA 和 高频 cD
  size_t nA = (n + 1) / 2;
  // 注意 nD = n - nA
  // 2) 上采样
  std::vector<double> upA(2 * n, 0.0), upD(2 * n, 0.0);
  for (size_t i = 0; i < nA; ++i) {
    upA[2*i] = signal[i];
  }
  for (size_t i = 0; i < n - nA; ++i) {
    upD[2*i] = signal[nA + i];
  }

  // 3) 延拓到 buf 长度 2*n + 2*pad
  std::vector<double> extA(2*n + 2*pad), extD(2*n + 2*pad);
  if (m_padding_mode == "periodic"c) {
    // 周期延拓
    for (int i = 0; i < pad; ++i) {
      int idx = (i - pad) % int(2*n);
      if (idx < 0) idx += int(2*n);
      extA[i] = upA[idx];
      extD[i] = upD[idx];
    }
    std::copy(upA.begin(), upA.end(), extA.begin() + pad);
    std::copy(upD.begin(), upD.end(), extD.begin() + pad);
    for (int i = 0; i < pad; ++i) {
      extA[2*n + pad + i] = upA[i % (2*n)];
      extD[2*n + pad + i] = upD[i % (2*n)];
    }
  } else {
    // 对称延拓
    for (int i = 0; i < pad; ++i) {
      extA[i] = upA[pad - 1 - i];
      extD[i] = upD[pad - 1 - i];
    }
    std::copy(upA.begin(), upA.end(), extA.begin() + pad);
    std::copy(upD.begin(), upD.end(), extD.begin() + pad);
    for (int i = 0; i < pad; ++i) {
      extA[2*n + pad + i] = upA[2*n - 1 - i];
      extD[2*n + pad + i] = upD[2*n - 1 - i];
    }
  }

  // 4) 卷积 + 合并
  std::vector<double> tmp(2 * n);
  for (int i = 0; i < int(2 * n); ++i) {
    double a = 0, d = 0;
    for (int k = 0; k < L; ++k) {
      a += sym13::rec_lo[k] * extA[i + k];
      d += sym13::rec_hi[k] * extD[i + k];
    }
    tmp[i] = a + d;
  }

  // 5) 写回（注意调用者必须保证 signal[] 空间至少能写入 2*n 个元素）
  std::copy(tmp.begin(), tmp.end(), signal);
}
