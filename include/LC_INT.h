#ifndef LC_INT_H
#define LC_INT_H

// Alternative backend that replaces SPECK's bitplane integer coding with
// the LC framework's `BIT_4 → RZE_4 → RZE_1` lossless pipeline + ZSTD outer
// wrap. See Sprint 8 worklog notes; pipeline was chosen by Phase 2 analysis
// + per-file search across the full Sprint 6 grid (only -0.7% gain from
// per-file optimization over this fixed choice; integration overhead of
// per-file pipeline selection not warranted).
//
// Payload layout mirrors HuffZstd_INT for consistency with SPECK_FLT's
// uint-width selection:
//   byte  0      num_bitplanes_eq    (= ceil(log2(max_mag+1)); 0 if all-zero)
//   bytes 1..8   uint64_t payload_bytes * 8     (lets get_stream_full_len work)
//   byte  9      format_version (= 1 for LC)
//   bytes 10..N  ZSTD frame wrapping the LC-pipeline output bytes

#include "SPECK_INT.h"

namespace sperr {

template <typename T>
class LC_INT : public SPECK_INT<T> {
  using uint_type = T;
  using vecui_type = std::vector<uint_type>;

 public:
  ~LC_INT() override = default;

  void encode() override;
  void decode() override;
  void use_bitstream(const void* p, size_t len) override;

  auto encoded_bitstream_len() const -> size_t override;
  void append_encoded_bitstream(vec8_type& buf) const override;

 protected:
  void m_clean_LIS() override {}
  void m_sorting_pass() override {}
  void m_initialize_lists() override {}

 private:
  vec8_type m_payload;
};

};  // namespace sperr

#endif
