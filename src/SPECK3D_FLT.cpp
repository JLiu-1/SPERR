#include "SPECK3D_FLT.h"
#include "HuffZstd_INT.h"
#include "LC_INT.h"
#include "SPECK3D_INT_DEC.h"
#include "SPECK3D_INT_ENC.h"

namespace {

template <typename T>
auto make_enc(sperr::IntBackend b) -> std::unique_ptr<sperr::SPECK_INT<T>>
{
  switch (b) {
    case sperr::IntBackend::HuffZstd: return std::make_unique<sperr::HuffZstd_INT<T>>();
    case sperr::IntBackend::LC:       return std::make_unique<sperr::LC_INT<T>>();
    default:                          return std::make_unique<sperr::SPECK3D_INT_ENC<T>>();
  }
}

template <typename T>
auto make_dec(sperr::IntBackend b) -> std::unique_ptr<sperr::SPECK_INT<T>>
{
  switch (b) {
    case sperr::IntBackend::HuffZstd: return std::make_unique<sperr::HuffZstd_INT<T>>();
    case sperr::IntBackend::LC:       return std::make_unique<sperr::LC_INT<T>>();
    default:                          return std::make_unique<sperr::SPECK3D_INT_DEC<T>>();
  }
}

}  // namespace

void sperr::SPECK3D_FLT::m_instantiate_encoder()
{
  switch (m_uint_flag) {
    case UINTType::UINT8:
      m_encoder = make_enc<uint8_t>(m_int_backend);
      break;
    case UINTType::UINT16:
      m_encoder = make_enc<uint16_t>(m_int_backend);
      break;
    case UINTType::UINT32:
      m_encoder = make_enc<uint32_t>(m_int_backend);
      break;
    default:
      m_encoder = make_enc<uint64_t>(m_int_backend);
  }
}

void sperr::SPECK3D_FLT::m_instantiate_decoder()
{
  switch (m_uint_flag) {
    case UINTType::UINT8:
      m_decoder = make_dec<uint8_t>(m_int_backend);
      break;
    case UINTType::UINT16:
      m_decoder = make_dec<uint16_t>(m_int_backend);
      break;
    case UINTType::UINT32:
      m_decoder = make_dec<uint32_t>(m_int_backend);
      break;
    default:
      m_decoder = make_dec<uint64_t>(m_int_backend);
  }
}

void sperr::SPECK3D_FLT::m_wavelet_xform()
{
  m_cdf.dwt3d();
}

void sperr::SPECK3D_FLT::m_inverse_wavelet_xform(bool multi_res)
{
  if (!multi_res)
    m_cdf.idwt3d();
  else
    m_cdf.idwt3d_multi_res(m_hierarchy);
}
