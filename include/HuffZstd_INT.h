#ifndef HUFFZSTD_INT_H
#define HUFFZSTD_INT_H

//
// Alternative backend that replaces SPECK's bitplane integer coding with
// Huffman + ZSTD (the lossless pipeline borrowed from SZ3).  See Phase 1/2
// notes in the project log; the algorithm is the SZ3 default with no
// parameter tuning.
//
// Layout of `m_payload` (also the exact bytes appended by
// `append_encoded_bitstream`):
//   byte  0      num_bitplanes_eq  (= ceil(log2(max_mag+1)); 0 if all-zero)
//                                  Identical role to SPECK's `num_bitplanes`
//                                  so SPECK_FLT's integer-type selection on
//                                  the decode side works unchanged.
//   bytes 1..8   uint64_t  payload_bytes * 8        (let SPECK_INT's
//                                                    `get_stream_full_len`
//                                                    formula yield the
//                                                    correct full length)
//   bytes 9..    ZSTD frame containing:
//                   [Huffman tree serialization]
//                   [Huffman bitstream of signed bins]
//

#include "SPECK_INT.h"

namespace sperr {

template <typename T>
class HuffZstd_INT : public SPECK_INT<T> {
  using uint_type = T;
  using vecui_type = std::vector<uint_type>;

 public:
  ~HuffZstd_INT() override = default;

  void encode() override;
  void decode() override;
  void use_bitstream(const void* p, size_t len) override;

  auto encoded_bitstream_len() const -> size_t override;
  void append_encoded_bitstream(vec8_type& buf) const override;

 protected:
  // The base class declares these as pure virtual.  HuffZstd never visits
  // them because we override encode() / decode() entirely.
  void m_clean_LIS() override {}
  void m_sorting_pass() override {}
  void m_initialize_lists() override {}

 private:
  // Materialised on encode(); appended verbatim by append_encoded_bitstream().
  // On decode, use_bitstream() copies the incoming bytes here before decode()
  // unpacks the Huffman+ZSTD stream into m_coeff_buf / m_sign_array.
  vec8_type m_payload;
};

};  // namespace sperr

#endif
