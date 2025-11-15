#include "CDF97.h"

#include <algorithm>
#include <cassert>
#include <numeric>  // std::accumulate()
#include <type_traits>

#ifdef __AVX2__
#include <immintrin.h>
#endif

// Destructor
sperr::CDF97::~CDF97()
{
  if (m_aligned_buf)
    std::free(m_aligned_buf);
   if (m_aligned_buf_f)
    std::free(m_aligned_buf_f);
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
  if (max_col * sizeof(double) > m_aligned_buf_bytes) {
    if (m_aligned_buf)
      std::free(m_aligned_buf);
    size_t alignment = 32;  // 256 bits
    size_t alloc_chunks = (max_col * 8 + 31) / alignment;
    m_aligned_buf_bytes = alignment * alloc_chunks;
    m_aligned_buf = static_cast<double*>(std::aligned_alloc(alignment, m_aligned_buf_bytes));

    size_t alignment_f = 64;  // 256 bits
    size_t alloc_chunks_f = (max_col * 4 + 63) / alignment_f;
    m_aligned_buf_f_bytes = alignment_f * alloc_chunks_f;
    m_aligned_buf_f = static_cast<float*>(std::aligned_alloc(alignment_f, m_aligned_buf_f_bytes));
  }

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[0] * dims[2]), dims[1] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(max_slice);
  if (max_slice > m_slice_buf_f.size())
    m_slice_buf_f.resize(max_slice);

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
  if (max_col * sizeof(double) > m_aligned_buf_bytes) {
    if (m_aligned_buf)
      std::free(m_aligned_buf);
    size_t alignment = 32;  // 256 bits
    size_t alloc_chunks = (max_col * 8 + 31) / alignment;
    m_aligned_buf_bytes = alignment * alloc_chunks;
    m_aligned_buf = static_cast<double*>(std::aligned_alloc(alignment, m_aligned_buf_bytes));
  }

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[0] * dims[2]), dims[1] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(max_slice);

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
      m_idwt2d_one_level(m_data_buf.data(), {x + xd, y + yd});
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
      m_sub_volume({x, y, z}, buf.data());
      m_idwt3d_one_level({x + xd, y + yd, z + zd});
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
    if(lev ==0)
      m_dwt3d_one_level_f({x, y, z});
    else
      m_dwt3d_one_level({x, y, z});
  }
}

void sperr::CDF97::m_idwt3d_dyadic(size_t num_xforms)
{
  for (size_t lev = num_xforms; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev - 1);
    auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev - 1);
    auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev - 1);
     //if(lev ==0)
    //  m_idwt3d_one_level_f({x, y, z});
    //else
      m_idwt3d_one_level({x, y, z});
  }
}

//
// Private Methods
//
void sperr::CDF97::m_dwt1d(double* array, size_t array_len, size_t num_of_lev)
{
  for (size_t lev = 0; lev < num_of_lev; lev++) {
    m_gather(array, array_len, m_aligned_buf);
    this->QccWAVCDF97AnalysisSymmetric(m_aligned_buf, array_len);
    std::copy(m_aligned_buf, m_aligned_buf + array_len, array);
    array_len -= array_len / 2;
  }
}

void sperr::CDF97::m_idwt1d(double* array, size_t array_len, size_t num_of_lev)
{
  for (size_t lev = num_of_lev; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(array_len, lev - 1);
    this->QccWAVCDF97SynthesisSymmetric(array, x);
    m_scatter(array, x, m_aligned_buf);
    std::copy(m_aligned_buf, m_aligned_buf + x, array);
  }
}

void sperr::CDF97::m_dwt2d(double* plane, std::array<size_t, 2> len_xy, size_t num_of_lev)
{
  for (size_t lev = 0; lev < num_of_lev; lev++) {
    auto [x, xd] = sperr::calc_approx_detail_len(len_xy[0], lev);
    auto [y, yd] = sperr::calc_approx_detail_len(len_xy[1], lev);
    m_dwt2d_one_level(plane, {x, y});
  }
}

void sperr::CDF97::m_idwt2d(double* plane, std::array<size_t, 2> len_xy, size_t num_of_lev)
{
  for (size_t lev = num_of_lev; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(len_xy[0], lev - 1);
    auto [y, yd] = sperr::calc_approx_detail_len(len_xy[1], lev - 1);
    m_idwt2d_one_level(plane, {x, y});
  }
}

void sperr::CDF97::m_dwt2d_one_level(double* plane, std::array<size_t, 2> len_xy)
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

void sperr::CDF97::m_dwt2d_one_level_f(double* plane, std::array<size_t, 2> len_xy)
{
  // First, perform DWT along X for every row
  for (size_t i = 0; i < len_xy[1]; i++) {
    auto* pos = plane + i * m_dims[0];
    m_gather(pos, len_xy[0], m_aligned_buf);
    this->QccWAVCDF97AnalysisSymmetric(m_aligned_buf, len_xy[0]);//this is double, currently
    for(size_t i = 0; i < len_xy[0]; i ++)
      pos[i] = m_aligned_buf_f[i];
  }

  // Second, perform DWT along Y for every column
  for (size_t x = 0; x < len_xy[0]; x++) {
    for (size_t y = 0; y < len_xy[1]; y++)
      m_slice_buf_f[y] = plane[y * m_dims[0] + x];
    m_gather_f(m_slice_buf_f.data(), len_xy[1], m_aligned_buf_f);
    this->QccWAVCDF97AnalysisSymmetric_f(m_aligned_buf_f, len_xy[1]);
    for (size_t y = 0; y < len_xy[1]; y++)
      plane[y * m_dims[0] + x] = m_aligned_buf_f[y];
  }
}



void sperr::CDF97::m_idwt2d_one_level(double* plane, std::array<size_t, 2> len_xy)
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

void sperr::CDF97::m_idwt2d_one_level_f(double* plane, std::array<size_t, 2> len_xy)
{
  // First, perform IDWT along Y for every column
  for (size_t x = 0; x < len_xy[0]; x++) {
    for (size_t y = 0; y < len_xy[1]; y++)
      m_slice_buf_f[y] = plane[y * m_dims[0] + x];
    this->QccWAVCDF97SynthesisSymmetric_f(m_slice_buf_f.data(), len_xy[1]);
    m_scatter_f(m_slice_buf_f.data(), len_xy[1], m_aligned_buf_f);
    for (size_t y = 0; y < len_xy[1]; y++)
      plane[y * m_dims[0] + x] = m_aligned_buf_f[y];
  }

  // Second, perform IDWT along X for every row
  for (size_t i = 0; i < len_xy[1]; i++) {
    auto* pos = plane + i * m_dims[0];
    this->QccWAVCDF97SynthesisSymmetric(pos, len_xy[0]);//this is double, currently
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
  // Note: the reason to process eight columns at a time is that a cache line
  // is usually 64 bytes, or 8 doubles. That means when you pay the cost to retrieve
  // one value from the Z column, its neighboring 7 values are available for free!

   for (size_t y = 0; y < len_xyz[1]; y++) {
    for (size_t x = 0; x < len_xyz[0]; x += 8) {
      const size_t xy_offset = y * m_dims[0] + x;
      const auto stride = std::min(8ul, len_xyz[0] - x);

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


void sperr::CDF97::m_dwt3d_one_level_f(std::array<size_t, 3> len_xyz)
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
  // Note: the reason to process eight columns at a time is that a cache line
  // is usually 64 bytes, or 16 floats. That means when you pay the cost to retrieve
  // one value from the Z column, its neighboring 15 values are available for free!

  for (size_t y = 0; y < len_xyz[1]; y++) {
    for (size_t x = 0; x < len_xyz[0]; x += 8) {
      const size_t xy_offset = y * m_dims[0] + x;
      const auto stride = std::min(16ul, len_xyz[0] - x);

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_slice_buf_f[z + i * col_len] = m_data_buf[z * plane_size_xy + xy_offset + i];
      }

      for (size_t i = 0; i < stride; i++) {
        auto* itr = m_slice_buf_f.data() + i * col_len;
        m_gather_f(itr, col_len, m_aligned_buf_f);
        this->QccWAVCDF97AnalysisSymmetric_f(m_aligned_buf_f, col_len);
        std::copy(m_aligned_buf_f, m_aligned_buf_f + col_len, itr);
      }

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_data_buf[z * plane_size_xy + xy_offset + i] = m_slice_buf_f[z + i * col_len];
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
  // Note: the reason to process eight columns at a time is that a cache line
  // is usually 64 bytes, or 8 doubles. That means when you pay the cost to retrieve
  // one value from the Z column, its neighboring 7 values are available for free!

  for (size_t y = 0; y < len_xyz[1]; y++) {
    for (size_t x = 0; x < len_xyz[0]; x += 8) {
      const size_t xy_offset = y * m_dims[0] + x;
      const auto stride = std::min(8ul, len_xyz[0] - x);

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

void sperr::CDF97::m_idwt3d_one_level_f(std::array<size_t, 3> len_xyz)
{
  const auto plane_size_xy = m_dims[0] * m_dims[1];
  const auto col_len = len_xyz[2];

  // First, do one level of inverse transform on all Z columns.  Strategy:
  // 1) extract eight Z columns to buffer space `m_slice_buf`
  // 2) transform these eight columns
  // 3) put the Z columns back to their appropriate locations in the volume.
  //
  // Note: the reason to process eight columns at a time is that a cache line
  // is usually 64 bytes, or 16 floats. That means when you pay the cost to retrieve
  // one value from the Z column, its neighboring 15 values are available for free!

  for (size_t y = 0; y < len_xyz[1]; y++) {
    for (size_t x = 0; x < len_xyz[0]; x += 8) {
      const size_t xy_offset = y * m_dims[0] + x;
      const auto stride = std::min(16ul, len_xyz[0] - x);

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_slice_buf_f[z + i * col_len] = m_data_buf[z * plane_size_xy + xy_offset + i];
      }

      for (size_t i = 0; i < stride; i++) {
        auto* itr = m_slice_buf_f.data() + i * col_len;
        this->QccWAVCDF97SynthesisSymmetric_f(itr, col_len);
        m_scatter_f(itr, col_len, m_aligned_buf_f);
        std::copy(m_aligned_buf_f, m_aligned_buf_f + col_len, itr);
      }

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_data_buf[z * plane_size_xy + xy_offset + i] = m_slice_buf_f[z + i * col_len];
      }
    }
  }

  // Second, do one level of inverse transform on all XY planes.
  for (size_t z = 0; z < len_xyz[2]; z++) {
    const size_t offset = plane_size_xy * z;
    m_idwt2d_one_level_f(m_data_buf.data() + offset, {len_xyz[0], len_xyz[1]});
  }
}



void sperr::CDF97::m_gather(const double* src, size_t len, double* dst) const
{
#ifdef __AVX2__
  const double* src_end = src + len;
  double* dst_evens = dst;
  double* dst_odds = dst + len - len / 2;

  // Process 8 elements at a time
  for (; src + 8 <= src_end; src += 8) {
    __m256d v0 = _mm256_loadu_pd(src);      // 0, 1, 2, 3
    __m256d v1 = _mm256_loadu_pd(src + 4);  // 4, 5, 6, 7

    __m256d evens = _mm256_unpacklo_pd(v0, v1);  // 0, 4, 2, 6
    __m256d odds = _mm256_unpackhi_pd(v0, v1);   // 1, 5, 3, 7

    __m256d result1 = _mm256_permute4x64_pd(evens, 0b11011000);  // 0, 2, 4, 6
    __m256d result2 = _mm256_permute4x64_pd(odds, 0b11011000);   // 1, 3, 5, 7

    _mm256_store_pd(dst_evens, result1);
    _mm256_storeu_pd(dst_odds, result2);

    dst_evens += 4;
    dst_odds += 4;
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

void sperr::CDF97::m_gather_f(const float* src, size_t len, float* dst) const
{
#ifdef __AVX2__
  const float* src_end = src + len;

  // 前半段存偶下标元素，后半段存奇下标
  const size_t even_count = len - len / 2;  // = ceil(len/2)
  float* dst_evens = dst;
  float* dst_odds  = dst + even_count;

  // idx_even: 取 0,2,4,6；idx_odd: 取 1,3,5,7
  const __m256i idx_even = _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);
  const __m256i idx_odd  = _mm256_setr_epi32(1, 3, 5, 7, 0, 0, 0, 0);

  // 每次处理 16 个 float：src[0..15]
  for (; src + 16 <= src_end; src += 16) {
    __m256 v0 = _mm256_loadu_ps(src);      // 0,1,2,3,4,5,6,7
    __m256 v1 = _mm256_loadu_ps(src + 8);  // 8,9,10,11,12,13,14,15

    // 从每个 8 元向量中取出 0,2,4,6（偶数索引）和 1,3,5,7（奇数索引）
    __m256 ev0_full = _mm256_permutevar8x32_ps(v0, idx_even);  // 0,2,4,6,*,*,*,*
    __m256 ev1_full = _mm256_permutevar8x32_ps(v1, idx_even);  // 8,10,12,14,*,*,*,*
    __m256 od0_full = _mm256_permutevar8x32_ps(v0, idx_odd);   // 1,3,5,7,*,*,*,*
    __m256 od1_full = _mm256_permutevar8x32_ps(v1, idx_odd);   // 9,11,13,15,*,*,*,*

    // 只需要每个向量的低 4 个元素
    __m128 ev0 = _mm256_castps256_ps128(ev0_full);  // 0,2,4,6
    __m128 ev1 = _mm256_castps256_ps128(ev1_full);  // 8,10,12,14
    __m128 od0 = _mm256_castps256_ps128(od0_full);  // 1,3,5,7
    __m128 od1 = _mm256_castps256_ps128(od1_full);  // 9,11,13,15

    // 依次写入偶数下标 / 奇数下标
    _mm_storeu_ps(dst_evens,      ev0);
    _mm_storeu_ps(dst_evens + 4,  ev1);
    _mm_storeu_ps(dst_odds,       od0);
    _mm_storeu_ps(dst_odds  + 4,  od1);

    dst_evens += 8;  // 本轮处理了 8 个偶数下标
    dst_odds  += 8;  // 本轮处理了 8 个奇数下标
  }

  // 尾部处理（最多剩下 0~15 个元素，标量按偶/奇拷贝）
  for (; src + 2 <= src_end; src += 2) {
    *(dst_evens++) = src[0];  // 偶数下标
    *(dst_odds++)  = src[1];  // 奇数下标
  }
  if (src < src_end) {
    // 剩下最后一个偶数下标元素
    *(dst_evens++) = *src;
  }

#else
  // 标量 fallback：和 double 版逻辑一致
  size_t low_count  = len - len / 2;  // 偶数下标个数
  size_t high_count = len / 2;        // 奇数下标个数
  const float* p = src;
  float* q = dst;

  for (size_t i = 0; i < low_count; ++i) {
    *q++ = p[i * 2];       // 0,2,4,...
  }
  for (size_t i = 0; i < high_count; ++i) {
    *q++ = p[i * 2 + 1];   // 1,3,5,...
  }
#endif
}



void sperr::CDF97::m_scatter(const double* begin, size_t len, double* dst) const
{
#ifdef __AVX2__
  const double* even_end = begin + len - len / 2;
  const double* odd_beg = even_end;
  const double* dst_end = dst + len;

  // Process 8 elements at a time
  for (; begin + 4 < even_end; begin += 4) {
    __m256d v0 = _mm256_loadu_pd(begin);    // 0, 1, 2, 3
    __m256d v1 = _mm256_loadu_pd(odd_beg);  // 4, 5, 6, 7

    __m256d evens = _mm256_unpacklo_pd(v0, v1);  // 0, 4, 2, 6
    __m256d odds = _mm256_unpackhi_pd(v0, v1);   // 1, 5, 3, 7

    __m256d result1 = _mm256_permute2f128_pd(evens, odds, 0x20);  // 0, 4, 1, 5
    __m256d result2 = _mm256_permute2f128_pd(evens, odds, 0x31);  // 2, 6, 3, 7

    _mm256_store_pd(dst, result1);
    _mm256_store_pd(dst + 4, result2);

    dst += 8;
    odd_beg += 4;
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

void sperr::CDF97::m_scatter_f(const float* begin, size_t len, float* dst) const
{
#ifdef __AVX2__
  const size_t even_count = len - len / 2;  // ceil(len/2)
  const size_t odd_count  = len / 2;        // floor(len/2)

  const float* even_ptr = begin;
  const float* odd_ptr  = begin + even_count;
  float* dst_end        = dst + len;

  // 每次处理 8 个偶 + 8 个奇 → 16 个输出
  const size_t vec_iters = std::min(even_count / 8, odd_count / 8);

  for (size_t it = 0; it < vec_iters; ++it) {
    // 直接按 128-bit 加载
    __m128 e_lo = _mm_loadu_ps(even_ptr);      // e0 e1 e2 e3
    __m128 e_hi = _mm_loadu_ps(even_ptr + 4);  // e4 e5 e6 e7
    __m128 o_lo = _mm_loadu_ps(odd_ptr);       // o0 o1 o2 o3
    __m128 o_hi = _mm_loadu_ps(odd_ptr + 4);   // o4 o5 o6 o7

    // 交错： [e0 o0 e1 o1], [e2 o2 e3 o3], [e4 o4 e5 o5], [e6 o6 e7 o7]
    __m128 d0 = _mm_unpacklo_ps(e_lo, o_lo); // e0 o0 e1 o1
    __m128 d1 = _mm_unpackhi_ps(e_lo, o_lo); // e2 o2 e3 o3
    __m128 d2 = _mm_unpacklo_ps(e_hi, o_hi); // e4 o4 e5 o5
    __m128 d3 = _mm_unpackhi_ps(e_hi, o_hi); // e6 o6 e7 o7

    _mm_storeu_ps(dst,       d0);
    _mm_storeu_ps(dst + 4,   d1);
    _mm_storeu_ps(dst + 8,   d2);
    _mm_storeu_ps(dst + 12,  d3);

    even_ptr += 8;
    odd_ptr  += 8;
    dst      += 16;
  }

  // scalar 收尾
  while (dst + 1 < dst_end && even_ptr < begin + even_count && odd_ptr < begin + len) {
    *dst       = *even_ptr++;  // 偶
    *(dst + 1) = *odd_ptr++;   // 奇
    dst += 2;
  }

  if (dst < dst_end && even_ptr < begin + even_count) {
    *dst = *even_ptr;
  }

#else
  size_t low_count  = len - len / 2;
  size_t high_count = len / 2;
  const float* p = begin;
  float* q = dst;

  for (size_t i = 0; i < low_count; ++i) {
    q[i * 2] = p[i];
  }
  for (size_t i = 0; i < high_count; ++i) {
    q[i * 2 + 1] = p[low_count + i];
  }
#endif
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

void sperr::CDF97::m_sub_volume(dims_type subdims, double* dst) const
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
void sperr::CDF97::QccWAVCDF97AnalysisSymmetric(double* signal, size_t len)
{

  size_t even_len = len - len / 2;
  size_t odd_len = len / 2;
  double* even = signal;
  double* odd = signal + even_len;
#ifdef __AVX2__
  if(len >= 8){
    const __m256d vALPHA        = _mm256_set1_pd(ALPHA);
    const __m256d vBETA         = _mm256_set1_pd(BETA);
    const __m256d v2BETA        = _mm256_set1_pd(2.0 * BETA);
    const __m256d vGAMMA        = _mm256_set1_pd(GAMMA);
    const __m256d vDELTA        = _mm256_set1_pd(DELTA);
    const __m256d v2DELTA       = _mm256_set1_pd(2.0 * DELTA);
    const __m256d vEPSILON      = _mm256_set1_pd(EPSILON);
    const __m256d vNEG_INV_EPS  = _mm256_set1_pd(-INV_EPSILON);

    // -------------------------
    // 1) Process all the odd elements: ALPHA step
    // odd[i] += ALPHA * (even[i] + even[i+1])
    // -------------------------
    {
      size_t i = 0;
      size_t vec_end = (odd_len > 1) ? ((odd_len - 1) & ~size_t(3)) : 0; // up to odd_len-1 (excluded)

      for (; i < vec_end; i += 4) {
        __m256d vOdd  = _mm256_loadu_pd(odd + i);
        __m256d vE0   = _mm256_loadu_pd(even + i);
        __m256d vE1   = _mm256_loadu_pd(even + i + 1);
        __m256d vSum  = _mm256_add_pd(vE0, vE1);
        __m256d vUpd  = _mm256_mul_pd(vALPHA, vSum);
        vOdd          = _mm256_add_pd(vOdd, vUpd);
        _mm256_storeu_pd(odd + i, vOdd);
      }
      // 尾部（直到 odd_len-2）
      for (; i < odd_len - 1; ++i) {
        odd[i] += ALPHA * (even[i] + even[i + 1]);
      }
      // 最后一个点的对称边界
      odd[odd_len - 1] += ALPHA * (even[odd_len - 1] + even[even_len - 1]);
    }

    // -------------------------
    // 2) Process all the even elements: BETA step
    // even[0] += 2*BETA*odd[0];
    // even[i] += BETA*(odd[i-1] + odd[i]);  i=1..even_len-2
    // even[even_len-1] += BETA*(odd[even_len-2] + odd[odd_len-1]);
    // -------------------------
    {
      even[0] += 2.0 * BETA * odd[0];

      if (even_len > 2) {
        size_t i_start = 1;
        // 内部区间的最大 i：even_len-2
        // 为了 vector load 安全，我们要保证 i+3 <= odd_len-1
        size_t i_max   = (even_len >= 2 ? even_len - 2 : 0);
        size_t vec_end = i_start;
        if (odd_len > 0 && i_max >= i_start) {
          size_t max_i_by_odd = (odd_len >= 4) ? (odd_len - 4) : i_start - 1;
          size_t limit        = (i_max < max_i_by_odd) ? i_max : max_i_by_odd;
          if (limit >= i_start)
            vec_end = ((limit - i_start + 1) & ~size_t(3)) + i_start;
        }

        size_t i = i_start;
        for (; i < vec_end; i += 4) {
          __m256d vEven = _mm256_loadu_pd(even + i);
          __m256d vO0   = _mm256_loadu_pd(odd + i - 1);
          __m256d vO1   = _mm256_loadu_pd(odd + i);
          __m256d vSum  = _mm256_add_pd(vO0, vO1);
          __m256d vUpd  = _mm256_mul_pd(vBETA, vSum);
          vEven         = _mm256_add_pd(vEven, vUpd);
          _mm256_storeu_pd(even + i, vEven);
        }
        // 剩余的 i=vec_end..even_len-2 标量处理
        for (; i < even_len - 1; ++i) {
          even[i] += BETA * (odd[i - 1] + odd[i]);
        }
      }

      if (even_len >= 2) {
        even[even_len - 1] += BETA * (odd[even_len - 2] + odd[odd_len - 1]);
      }
    }

    // -------------------------
    // 3) Process all the odd elements: GAMMA step
    // odd[i] += GAMMA*(even[i] + even[i+1])
    // -------------------------
    {
      size_t i = 0;
      size_t vec_end = (odd_len > 1) ? ((odd_len - 1) & ~size_t(3)) : 0;

      for (; i < vec_end; i += 4) {
        __m256d vOdd  = _mm256_loadu_pd(odd + i);
        __m256d vE0   = _mm256_loadu_pd(even + i);
        __m256d vE1   = _mm256_loadu_pd(even + i + 1);
        __m256d vSum  = _mm256_add_pd(vE0, vE1);
        __m256d vUpd  = _mm256_mul_pd(vGAMMA, vSum);
        vOdd          = _mm256_add_pd(vOdd, vUpd);
        _mm256_storeu_pd(odd + i, vOdd);
      }
      for (; i < odd_len - 1; ++i) {
        odd[i] += GAMMA * (even[i] + even[i + 1]);
      }
      odd[odd_len - 1] += GAMMA * (even[odd_len - 1] + even[even_len - 1]);
    }

    // -------------------------
    // 4) Process all the even elements: DELTA + scaling EPSILON
    // even[0] = EPSILON * (even[0] + 2*DELTA*odd[0]);
    // even[i] = EPSILON * (even[i] + DELTA*(odd[i-1] + odd[i]));
    // even[even_len-1] = EPSILON * (even[...] + DELTA*(odd[...] + odd[...]));
    // -------------------------
    {
      even[0] = EPSILON * (even[0] + 2.0 * DELTA * odd[0]);

      if (even_len > 2) {
        size_t i_start = 1;
        size_t i_max   = (even_len >= 2 ? even_len - 2 : 0);
        size_t vec_end = i_start;
        if (odd_len > 0 && i_max >= i_start) {
          size_t max_i_by_odd = (odd_len >= 4) ? (odd_len - 4) : i_start - 1;
          size_t limit        = (i_max < max_i_by_odd) ? i_max : max_i_by_odd;
          if (limit >= i_start)
            vec_end = ((limit - i_start + 1) & ~size_t(3)) + i_start;
        }

        size_t i = i_start;
        for (; i < vec_end; i += 4) {
          __m256d vEven = _mm256_loadu_pd(even + i);
          __m256d vO0   = _mm256_loadu_pd(odd + i - 1);
          __m256d vO1   = _mm256_loadu_pd(odd + i);
          __m256d vSum  = _mm256_add_pd(vO0, vO1);
          __m256d vUpd  = _mm256_mul_pd(vDELTA, vSum);
          vEven         = _mm256_add_pd(vEven, vUpd);
          vEven         = _mm256_mul_pd(vEven, vEPSILON);
          _mm256_storeu_pd(even + i, vEven);
        }
        for (; i < even_len - 1; ++i) {
          even[i] = EPSILON * (even[i] + DELTA * (odd[i - 1] + odd[i]));
        }
      }

      if (even_len >= 2) {
        even[even_len - 1] =
            EPSILON * (even[even_len - 1] + DELTA * (odd[even_len - 2] + odd[odd_len - 1]));
      }
    }

    // -------------------------
    // 5) Process odd elements: scaling by -INV_EPSILON
    // odd[i] *= -INV_EPSILON;
    // -------------------------
    {
      size_t i = 0;
      size_t vec_end = (odd_len & ~size_t(3)); // multiple of 4
      for (; i < vec_end; i += 4) {
        __m256d vOdd = _mm256_loadu_pd(odd + i);
        vOdd         = _mm256_mul_pd(vOdd, vNEG_INV_EPS);
        _mm256_storeu_pd(odd + i, vOdd);
      }
      for (; i < odd_len; ++i) {
        odd[i] *= -INV_EPSILON;
      }
    }
    return;
  } 

#endif
  

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


void sperr::CDF97::QccWAVCDF97AnalysisSymmetric_f(float* signal, size_t len)
{
  size_t even_len = len - len / 2;
  size_t odd_len  = len / 2;
  float* even     = signal;
  float* odd      = signal + even_len;

#ifdef __AVX2__
  if (len >= 16)
  {
    const __m256 vALPHA       = _mm256_set1_ps(ALPHA_f);
    const __m256 vBETA        = _mm256_set1_ps(BETA_f);
    const __m256 v2BETA       = _mm256_set1_ps(2.0f * BETA_f);
    const __m256 vGAMMA       = _mm256_set1_ps(GAMMA_f);
    const __m256 vDELTA       = _mm256_set1_ps(DELTA_f);
    const __m256 v2DELTA      = _mm256_set1_ps(2.0f * DELTA_f);
    const __m256 vEPSILON     = _mm256_set1_ps(EPSILON_f);
    const __m256 vNEG_INV_EPS = _mm256_set1_ps(-INV_EPSILON_f);

    // -------------------------
    // 1) odd: ALPHA
    // odd[i] += ALPHA_f * (even[i] + even[i+1]);
    // -------------------------
    {
      size_t i = 0;
      size_t vec_end = (odd_len > 1) ? ((odd_len - 1) & ~size_t(7)) : 0;

      for (; i < vec_end; i += 8) {
        __m256 vOdd = _mm256_loadu_ps(odd + i);
        __m256 vE0  = _mm256_loadu_ps(even + i);
        __m256 vE1  = _mm256_loadu_ps(even + i + 1);
        vOdd        = _mm256_add_ps(vOdd, _mm256_mul_ps(vALPHA, _mm256_add_ps(vE0, vE1)));
        _mm256_storeu_ps(odd + i, vOdd);
      }
      for (; i < odd_len - 1; i++) {
        odd[i] += ALPHA_f * (even[i] + even[i + 1]);
      }
      if (odd_len > 0)
        odd[odd_len - 1] += ALPHA_f * (even[odd_len - 1] + even[even_len - 1]);
    }

    // -------------------------
    // 2) even: BETA
    // -------------------------
    {
      if (even_len > 0)
        even[0] += 2.0f * BETA_f * odd[0];

      if (even_len > 2) {
        size_t i_start = 1;
        size_t i_max = even_len - 2;

        size_t vec_end = i_start;
        if (odd_len >= 8) {
          size_t max_i_by_odd = odd_len - 8;
          size_t limit = (i_max < max_i_by_odd) ? i_max : max_i_by_odd;
          if (limit >= i_start)
            vec_end = ((limit - i_start + 1) & ~size_t(7)) + i_start;
        }

        size_t i = i_start;
        for (; i < vec_end; i += 8) {
          __m256 vEven = _mm256_loadu_ps(even + i);
          __m256 vO0   = _mm256_loadu_ps(odd + i - 1);
          __m256 vO1   = _mm256_loadu_ps(odd + i);
          __m256 vSum  = _mm256_add_ps(vO0, vO1);
          vEven        = _mm256_add_ps(vEven, _mm256_mul_ps(vBETA, vSum));
          _mm256_storeu_ps(even + i, vEven);
        }

        for (; i < even_len - 1; ++i)
          even[i] += BETA_f * (odd[i - 1] + odd[i]);
      }

      if (even_len >= 2)
        even[even_len - 1] += BETA_f * (odd[even_len - 2] + odd[odd_len - 1]);
    }

    // -------------------------
    // 3) odd: GAMMA
    // -------------------------
    {
      size_t i = 0;
      size_t vec_end = (odd_len > 1) ? ((odd_len - 1) & ~size_t(7)) : 0;

      for (; i < vec_end; i += 8) {
        __m256 vOdd = _mm256_loadu_ps(odd + i);
        __m256 vE0  = _mm256_loadu_ps(even + i);
        __m256 vE1  = _mm256_loadu_ps(even + i + 1);
        vOdd        = _mm256_add_ps(vOdd, _mm256_mul_ps(vGAMMA, _mm256_add_ps(vE0, vE1)));
        _mm256_storeu_ps(odd + i, vOdd);
      }
      for (; i < odd_len - 1; i++)
        odd[i] += GAMMA_f * (even[i] + even[i + 1]);
      if (odd_len > 0)
        odd[odd_len - 1] += GAMMA_f * (even[odd_len - 1] + even[even_len - 1]);
    }

    // -------------------------
    // 4) even: DELTA + EPSILON
    // -------------------------
    {
      if (even_len > 0)
        even[0] = EPSILON_f * (even[0] + 2.0f * DELTA_f * odd[0]);

      if (even_len > 2) {
        size_t i_start = 1;
        size_t i_max = even_len - 2;

        size_t vec_end = i_start;
        if (odd_len >= 8) {
          size_t max_i_by_odd = odd_len - 8;
          size_t limit = (i_max < max_i_by_odd) ? i_max : max_i_by_odd;
          if (limit >= i_start)
            vec_end = ((limit - i_start + 1) & ~size_t(7)) + i_start;
        }

        size_t i = i_start;
        for (; i < vec_end; i += 8) {
          __m256 vEven = _mm256_loadu_ps(even + i);
          __m256 vO0   = _mm256_loadu_ps(odd + i - 1);
          __m256 vO1   = _mm256_loadu_ps(odd + i);
          __m256 vSum  = _mm256_add_ps(vO0, vO1);
          vEven        = _mm256_add_ps(vEven, _mm256_mul_ps(vDELTA, vSum));
          vEven        = _mm256_mul_ps(vEven, vEPSILON);
          _mm256_storeu_ps(even + i, vEven);
        }

        for (; i < even_len - 1; i++)
          even[i] = EPSILON_f * (even[i] + DELTA_f * (odd[i - 1] + odd[i]));
      }

      if (even_len >= 2)
        even[even_len - 1] =
          EPSILON_f * (even[even_len - 1] +
          DELTA_f * (odd[even_len - 2] + odd[odd_len - 1]));
    }

    // -------------------------
    // 5) odd: scale by -INV_EPSILON
    // -------------------------
    {
      size_t i = 0;
      size_t vec_end = odd_len & ~size_t(7);
      for (; i < vec_end; i += 8) {
        __m256 vOdd = _mm256_loadu_ps(odd + i);
        vOdd        = _mm256_mul_ps(vOdd, vNEG_INV_EPS);
        _mm256_storeu_ps(odd + i, vOdd);
      }
      for (; i < odd_len; i++)
        odd[i] *= -INV_EPSILON_f;
    }

    return;
  }
#endif

  // -------------------------  
  // 标量 fallback（全部使用 *_f 常量）
  // -------------------------

  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] += ALPHA_f * (even[i] + even[i + 1]);
  if (odd_len > 0)
    odd[odd_len - 1] += ALPHA_f * (even[odd_len - 1] + even[even_len - 1]);

  if (even_len > 0)
    even[0] += 2.0f * BETA_f * odd[0];
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] += BETA_f * (odd[i - 1] + odd[i]);
  if (even_len >= 2)
    even[even_len - 1] += BETA_f * (odd[even_len - 2] + odd[odd_len - 1]);

  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] += GAMMA_f * (even[i] + even[i + 1]);
  if (odd_len > 0)
    odd[odd_len - 1] += GAMMA_f * (even[odd_len - 1] + even[even_len - 1]);

  if (even_len > 0)
    even[0] = EPSILON_f * (even[0] + 2.0f * DELTA_f * odd[0]);
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] = EPSILON_f * (even[i] + DELTA_f * (odd[i - 1] + odd[i]));
  if (even_len >= 2)
    even[even_len - 1] =
      EPSILON_f * (even[even_len - 1] + DELTA_f * (odd[even_len - 2] + odd[odd_len - 1]));

  for (size_t i = 0; i < odd_len; i++)
    odd[i] *= -INV_EPSILON_f;
}






void sperr::CDF97::QccWAVCDF97SynthesisSymmetric(double* signal, size_t len)
{
  size_t even_len = len - len / 2;
  size_t odd_len = len / 2;
  double* even = signal;
  double* odd = signal + even_len;
  #ifdef __AVX2__
  if(len >= 8){

    const __m256d vNEG_EPS   = _mm256_set1_pd(-EPSILON);
    const __m256d vINV_EPS   = _mm256_set1_pd(INV_EPSILON);
    const __m256d vDELTA     = _mm256_set1_pd(DELTA);
    const __m256d v2DELTA    = _mm256_set1_pd(2.0 * DELTA);
    const __m256d vGAMMA     = _mm256_set1_pd(GAMMA);
    const __m256d vBETA      = _mm256_set1_pd(BETA);
    const __m256d v2BETA     = _mm256_set1_pd(2.0 * BETA);
    const __m256d vALPHA     = _mm256_set1_pd(ALPHA);
      // ----------------------------------------------------
    // 1) odd[i] *= (-EPSILON)
    // ----------------------------------------------------
    {
      size_t i = 0;
      size_t vec_end = odd_len & ~size_t(3); // 向下取 4 的倍数
      for (; i < vec_end; i += 4) {
        __m256d vOdd = _mm256_loadu_pd(odd + i);
        vOdd         = _mm256_mul_pd(vOdd, vNEG_EPS);
        _mm256_storeu_pd(odd + i, vOdd);
      }
      for (; i < odd_len; ++i) {
        odd[i] *= (-EPSILON);
      }
    }

    // ----------------------------------------------------
    // 2) even：
    // even[0] = even[0] * INV_EPS - 2*DELTA*odd[0];
    // even[i] = even[i] * INV_EPS - DELTA*(odd[i-1] + odd[i]), i=1..even_len-2
    // even[last] = ...
    // 中间那段尽量 AVX 化
    // ----------------------------------------------------
    {
      // 边界 i = 0
      even[0] = even[0] * INV_EPSILON - 2.0 * DELTA * odd[0];

      if (even_len > 2) {
        size_t i_start = 1;
        size_t i_max   = even_len - 2;          // 内部 i 最大值
        // odd 访问 odd[i-1], odd[i]，向量化时 load odd+(i-1) 和 odd+i
        // 要保证 i+3 <= odd_len-1，并且 i-1 >= 0
        size_t max_i_by_odd = (odd_len >= 4) ? (odd_len - 1 - 3) : 0;
        if (max_i_by_odd < i_start) max_i_by_odd = i_start - 1;
        size_t limit = (i_max < max_i_by_odd) ? i_max : max_i_by_odd;

        size_t vec_end = i_start;
        if (limit >= i_start) {
          vec_end = ((limit - i_start + 1) & ~size_t(3)) + i_start; // 4 对齐
        }

        size_t i = i_start;
        for (; i < vec_end; i += 4) {
          __m256d vEven = _mm256_loadu_pd(even + i);
          __m256d vO0   = _mm256_loadu_pd(odd + i - 1);
          __m256d vO1   = _mm256_loadu_pd(odd + i);
          __m256d vSumO = _mm256_add_pd(vO0, vO1);

          vEven = _mm256_mul_pd(vEven, vINV_EPS);      // even * INV_EPS
          __m256d vUpd = _mm256_mul_pd(vDELTA, vSumO); // DELTA*(odd[i-1]+odd[i])
          vEven = _mm256_sub_pd(vEven, vUpd);

          _mm256_storeu_pd(even + i, vEven);
        }

        // 剩余的一点用标量
        for (; i < even_len - 1; ++i) {
          even[i] = even[i] * INV_EPSILON - DELTA * (odd[i - 1] + odd[i]);
        }
      }

      // 右端边界
      if (even_len >= 2) {
        even[even_len - 1] =
            even[even_len - 1] * INV_EPSILON -
            DELTA * (odd[even_len - 2] + odd[odd_len - 1]);
      }
    }

    // ----------------------------------------------------
    // 3) odd:
    // odd[i] -= GAMMA*(even[i] + even[i+1]), i=0..odd_len-2
    // odd[last] -= GAMMA*(even[odd_len-1] + even[even_len-1])
    // 尽量 AVX
    // ----------------------------------------------------
    {
      if (odd_len > 1 && even_len > 1) {
        size_t max_i  = (odd_len >= 2 ? odd_len - 2 : 0);
        size_t max_i2 = (even_len >= 2 ? even_len - 2 : 0);
        if (max_i2 < max_i) max_i = max_i2;         // i+1 <= even_len-1

        size_t vec_end = ((max_i + 1) & ~size_t(3)); // 4 对齐

        size_t i = 0;
        for (; i < vec_end; i += 4) {
          __m256d vOdd  = _mm256_loadu_pd(odd  + i);
          __m256d vE0   = _mm256_loadu_pd(even + i);
          __m256d vE1   = _mm256_loadu_pd(even + i + 1);
          __m256d vSumE = _mm256_add_pd(vE0, vE1);
          __m256d vUpd  = _mm256_mul_pd(vGAMMA, vSumE);
          vOdd = _mm256_sub_pd(vOdd, vUpd);
          _mm256_storeu_pd(odd + i, vOdd);
        }
        // 尾部直到 odd_len-2
        for (; i < odd_len - 1; ++i) {
          odd[i] -= GAMMA * (even[i] + even[i + 1]);
        }
      } else if (odd_len > 1) {
        // 极小长度直接标量
        for (size_t i = 0; i < odd_len - 1; ++i) {
          odd[i] -= GAMMA * (even[i] + even[i + 1]);
        }
      }

      // 最后一个 odd 的边界
      if (odd_len > 0) {
        odd[odd_len - 1] -= GAMMA * (even[odd_len - 1] + even[even_len - 1]);
      }
    }

    // ----------------------------------------------------
    // 4) even:
    // even[0] -= 2*BETA*odd[0];
    // even[i] -= BETA*(odd[i-1] + odd[i]), i=1..even_len-2
    // even[last] -= BETA*(odd[even_len-2] + odd[odd_len-1])
    // 中间段 AVX
    // ----------------------------------------------------
    {
      // 左边界
      even[0] -= 2.0 * BETA * odd[0];

      if (even_len > 2) {
        size_t i_start = 1;
        size_t i_max   = even_len - 2;
        size_t max_i_by_odd = (odd_len >= 2 ? odd_len - 1 : 0); // i <= odd_len-1
        if (max_i_by_odd < i_start) max_i_by_odd = i_start - 1;
        size_t limit = (i_max < max_i_by_odd) ? i_max : max_i_by_odd;

        size_t vec_end = i_start;
        if (limit >= i_start) {
          vec_end = ((limit - i_start + 1) & ~size_t(3)) + i_start; // 4 对齐
        }

        size_t i = i_start;
        for (; i < vec_end; i += 4) {
          __m256d vEven = _mm256_loadu_pd(even + i);
          __m256d vO0   = _mm256_loadu_pd(odd + i - 1);
          __m256d vO1   = _mm256_loadu_pd(odd + i);
          __m256d vSumO = _mm256_add_pd(vO0, vO1);
          __m256d vUpd  = _mm256_mul_pd(vBETA, vSumO);
          vEven = _mm256_sub_pd(vEven, vUpd);
          _mm256_storeu_pd(even + i, vEven);
        }

        for (; i < even_len - 1; ++i) {
          even[i] -= BETA * (odd[i - 1] + odd[i]);
        }
      }

      // 右边界
      if (even_len >= 2) {
        even[even_len - 1] -=
            BETA * (odd[even_len - 2] + odd[odd_len - 1]);
      }
    }

    // ----------------------------------------------------
    // 5) odd:
    // odd[i] -= ALPHA*(even[i] + even[i+1]), i=0..odd_len-2
    // odd[last] -= ALPHA*(even[odd_len-1] + even[even_len-1])
    // AVX
    // ----------------------------------------------------
    {
      if (odd_len > 1 && even_len > 1) {
        size_t max_i  = (odd_len >= 2 ? odd_len - 2 : 0);
        size_t max_i2 = (even_len >= 2 ? even_len - 2 : 0);
        if (max_i2 < max_i) max_i = max_i2;

        size_t vec_end = ((max_i + 1) & ~size_t(3)); // 4 对齐

        size_t i = 0;
        for (; i < vec_end; i += 4) {
          __m256d vOdd  = _mm256_loadu_pd(odd + i);
          __m256d vE0   = _mm256_loadu_pd(even + i);
          __m256d vE1   = _mm256_loadu_pd(even + i + 1);
          __m256d vSumE = _mm256_add_pd(vE0, vE1);
          __m256d vUpd  = _mm256_mul_pd(vALPHA, vSumE);
          vOdd = _mm256_sub_pd(vOdd, vUpd);
          _mm256_storeu_pd(odd + i, vOdd);
        }
        for (; i < odd_len - 1; ++i) {
          odd[i] -= ALPHA * (even[i] + even[i + 1]);
        }
      } else if (odd_len > 1) {
        for (size_t i = 0; i < odd_len - 1; ++i) {
          odd[i] -= ALPHA * (even[i] + even[i + 1]);
        }
      }

      if (odd_len > 0) {
        odd[odd_len - 1] -= ALPHA * (even[odd_len - 1] + even[even_len - 1]);
      }
    }
    return;
  }

  #endif
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

void sperr::CDF97::QccWAVCDF97SynthesisSymmetric_f(float* signal, size_t len)
{
  size_t even_len = len - len / 2;
  size_t odd_len  = len / 2;
  float* even     = signal;
  float* odd      = signal + even_len;

#ifdef __AVX2__
  if (len >= 8) {

    const __m256 vNEG_EPS = _mm256_set1_ps(-EPSILON_f);
    const __m256 vINV_EPS = _mm256_set1_ps(INV_EPSILON_f);
    const __m256 vDELTA   = _mm256_set1_ps(DELTA_f);
    const __m256 v2DELTA  = _mm256_set1_ps(2.0f * DELTA_f);
    const __m256 vGAMMA   = _mm256_set1_ps(GAMMA_f);
    const __m256 vBETA    = _mm256_set1_ps(BETA_f);
    const __m256 v2BETA   = _mm256_set1_ps(2.0f * BETA_f);
    const __m256 vALPHA   = _mm256_set1_ps(ALPHA_f);

    // ----------------------------------------------------
    // 1) odd[i] *= (-EPSILON_f)
    // ----------------------------------------------------
    {
      size_t i = 0;
      size_t vec_end = odd_len & ~size_t(7); // 向下取 8 的倍数
      for (; i < vec_end; i += 8) {
        __m256 vOdd = _mm256_loadu_ps(odd + i);
        vOdd        = _mm256_mul_ps(vOdd, vNEG_EPS);
        _mm256_storeu_ps(odd + i, vOdd);
      }
      for (; i < odd_len; ++i) {
        odd[i] *= -EPSILON_f;
      }
    }

    // ----------------------------------------------------
    // 2) even:
    // even[0] = even[0] * INV_EPSILON_f - 2*DELTA_f*odd[0];
    // even[i] = even[i] * INV_EPSILON_f - DELTA_f*(odd[i-1] + odd[i]), i=1..even_len-2
    // even[last] 同理
    // ----------------------------------------------------
    {
      if (even_len > 0) {
        even[0] = even[0] * INV_EPSILON_f - 2.0f * DELTA_f * odd[0];
      }

      if (even_len > 2) {
        size_t i_start = 1;
        size_t i_max   = even_len - 2;          // 内部 i 最大值

        // odd 访问 odd[i-1], odd[i]，向量化时 load odd+(i-1) 和 odd+i
        // 要保证 i+7 <= odd_len-1（因为我们一次处理 8 个）
        size_t max_i_by_odd = 0;
        if (odd_len >= 8) {
          max_i_by_odd = odd_len - 1 - 7;
        }
        if (max_i_by_odd < i_start) {
          max_i_by_odd = i_start - 1;
        }

        size_t limit   = (i_max < max_i_by_odd) ? i_max : max_i_by_odd;
        size_t vec_end = i_start;
        if (limit >= i_start) {
          vec_end = ((limit - i_start + 1) & ~size_t(7)) + i_start; // 8 对齐
        }

        size_t i = i_start;
        for (; i < vec_end; i += 8) {
          __m256 vEven = _mm256_loadu_ps(even + i);
          __m256 vO0   = _mm256_loadu_ps(odd + i - 1);
          __m256 vO1   = _mm256_loadu_ps(odd + i);
          __m256 vSumO = _mm256_add_ps(vO0, vO1);

          vEven = _mm256_mul_ps(vEven, vINV_EPS);    // even * INV_EPS
          __m256 vUpd = _mm256_mul_ps(vDELTA, vSumO); // DELTA_f*(odd[i-1]+odd[i])
          vEven = _mm256_sub_ps(vEven, vUpd);

          _mm256_storeu_ps(even + i, vEven);
        }

        // 剩余的一点用标量
        for (; i < even_len - 1; ++i) {
          even[i] = even[i] * INV_EPSILON_f - DELTA_f * (odd[i - 1] + odd[i]);
        }
      }

      // 右端边界
      if (even_len >= 2) {
        even[even_len - 1] =
            even[even_len - 1] * INV_EPSILON_f -
            DELTA_f * (odd[even_len - 2] + odd[odd_len - 1]);
      }
    }

    // ----------------------------------------------------
    // 3) odd:
    // odd[i] -= GAMMA_f*(even[i] + even[i+1]), i=0..odd_len-2
    // odd[last] -= GAMMA_f*(even[odd_len-1] + even[even_len-1])
    // ----------------------------------------------------
    {
      if (odd_len > 1 && even_len > 1) {
        size_t max_i  = (odd_len >= 2 ? odd_len - 2 : 0);
        size_t max_i2 = (even_len >= 2 ? even_len - 2 : 0);
        if (max_i2 < max_i) max_i = max_i2;         // i+1 <= even_len-1

        size_t vec_end = ((max_i + 1) & ~size_t(7)); // 8 对齐

        size_t i = 0;
        for (; i < vec_end; i += 8) {
          __m256 vOdd  = _mm256_loadu_ps(odd  + i);
          __m256 vE0   = _mm256_loadu_ps(even + i);
          __m256 vE1   = _mm256_loadu_ps(even + i + 1);
          __m256 vSumE = _mm256_add_ps(vE0, vE1);
          __m256 vUpd  = _mm256_mul_ps(vGAMMA, vSumE);
          vOdd         = _mm256_sub_ps(vOdd, vUpd);
          _mm256_storeu_ps(odd + i, vOdd);
        }
        // 尾部直到 odd_len-2
        for (; i < odd_len - 1; ++i) {
          odd[i] -= GAMMA_f * (even[i] + even[i + 1]);
        }
      } else if (odd_len > 1) {
        // 极小长度直接标量
        for (size_t i = 0; i < odd_len - 1; ++i) {
          odd[i] -= GAMMA_f * (even[i] + even[i + 1]);
        }
      }

      // 最后一个 odd 的边界
      if (odd_len > 0) {
        odd[odd_len - 1] -= GAMMA_f * (even[odd_len - 1] + even[even_len - 1]);
      }
    }

    // ----------------------------------------------------
    // 4) even:
    // even[0] -= 2*BETA_f*odd[0];
    // even[i] -= BETA_f*(odd[i-1] + odd[i]), i=1..even_len-2
    // even[last] -= BETA_f*(odd[even_len-2] + odd[odd_len-1])
    // ----------------------------------------------------
    {
      if (even_len > 0)
        even[0] -= 2.0f * BETA_f * odd[0];

      if (even_len > 2) {
        size_t i_start = 1;
        size_t i_max   = even_len - 2;
        size_t max_i_by_odd = (odd_len >= 2 ? odd_len - 1 : 0); // i <= odd_len-1
        if (max_i_by_odd < i_start) max_i_by_odd = i_start - 1;
        size_t limit = (i_max < max_i_by_odd) ? i_max : max_i_by_odd;

        size_t vec_end = i_start;
        if (limit >= i_start) {
          vec_end = ((limit - i_start + 1) & ~size_t(7)) + i_start; // 8 对齐
        }

        size_t i = i_start;
        for (; i < vec_end; i += 8) {
          __m256 vEven = _mm256_loadu_ps(even + i);
          __m256 vO0   = _mm256_loadu_ps(odd + i - 1);
          __m256 vO1   = _mm256_loadu_ps(odd + i);
          __m256 vSumO = _mm256_add_ps(vO0, vO1);
          __m256 vUpd  = _mm256_mul_ps(vBETA, vSumO);
          vEven        = _mm256_sub_ps(vEven, vUpd);
          _mm256_storeu_ps(even + i, vEven);
        }

        for (; i < even_len - 1; ++i) {
          even[i] -= BETA_f * (odd[i - 1] + odd[i]);
        }
      }

      // 右边界
      if (even_len >= 2) {
        even[even_len - 1] -=
            BETA_f * (odd[even_len - 2] + odd[odd_len - 1]);
      }
    }

    // ----------------------------------------------------
    // 5) odd:
    // odd[i] -= ALPHA_f*(even[i] + even[i+1]), i=0..odd_len-2
    // odd[last] -= ALPHA_f*(even[odd_len-1] + even[even_len-1])
    // ----------------------------------------------------
    {
      if (odd_len > 1 && even_len > 1) {
        size_t max_i  = (odd_len >= 2 ? odd_len - 2 : 0);
        size_t max_i2 = (even_len >= 2 ? even_len - 2 : 0);
        if (max_i2 < max_i) max_i = max_i2;

        size_t vec_end = ((max_i + 1) & ~size_t(7)); // 8 对齐

        size_t i = 0;
        for (; i < vec_end; i += 8) {
          __m256 vOdd  = _mm256_loadu_ps(odd + i);
          __m256 vE0   = _mm256_loadu_ps(even + i);
          __m256 vE1   = _mm256_loadu_ps(even + i + 1);
          __m256 vSumE = _mm256_add_ps(vE0, vE1);
          __m256 vUpd  = _mm256_mul_ps(vALPHA, vSumE);
          vOdd         = _mm256_sub_ps(vOdd, vUpd);
          _mm256_storeu_ps(odd + i, vOdd);
        }
        for (; i < odd_len - 1; ++i) {
          odd[i] -= ALPHA_f * (even[i] + even[i + 1]);
        }
      } else if (odd_len > 1) {
        for (size_t i = 0; i < odd_len - 1; ++i) {
          odd[i] -= ALPHA_f * (even[i] + even[i + 1]);
        }
      }

      if (odd_len > 0) {
        odd[odd_len - 1] -= ALPHA_f * (even[odd_len - 1] + even[even_len - 1]);
      }
    }

    return;
  }
#endif

  // ===== 标量 fallback（用 *_f 常量） =====

  for (size_t i = 0; i < odd_len; i++)
    odd[i] *= -EPSILON_f;

  if (even_len > 0)
    even[0] = even[0] * INV_EPSILON_f - 2.0f * DELTA_f * odd[0];
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] = even[i] * INV_EPSILON_f - DELTA_f * (odd[i - 1] + odd[i]);
  if (even_len >= 2)
    even[even_len - 1] =
        even[even_len - 1] * INV_EPSILON_f
        - DELTA_f * (odd[even_len - 2] + odd[odd_len - 1]);

  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] -= GAMMA_f * (even[i] + even[i + 1]);
  if (odd_len > 0)
    odd[odd_len - 1] -= GAMMA_f * (even[odd_len - 1] + even[even_len - 1]);

  if (even_len > 0)
    even[0] -= 2.0f * BETA_f * odd[0];
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] -= BETA_f * (odd[i - 1] + odd[i]);
  if (even_len >= 2)
    even[even_len - 1] -= BETA_f * (odd[even_len - 2] + odd[odd_len - 1]);

  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] -= ALPHA_f * (even[i] + even[i + 1]);
  if (odd_len > 0)
    odd[odd_len - 1] -= ALPHA_f * (even[odd_len - 1] + even[even_len - 1]);
}

