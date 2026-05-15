// In-memory wrapper for the LC pipeline `BIT_4 → RZE_4 → RZE_1`.
// Derived from LC-framework's compressor-standalone.cpp / decompressor-standalone.cpp
// (BSD-3, Texas State University). Differences vs the original standalone:
//   - no file I/O — operates on caller-owned in/out byte buffers
//   - no main(), no perf-analysis prints
//   - functions return the encoded/decoded byte length
// CS (chunk size) is fixed at 16 KiB to match the LC standalone defaults; this
// is what Phase 0–2 of Sprint 8 measured against.

#ifndef SPERR_LC_PIPELINE_HPP
#define SPERR_LC_PIPELINE_HPP

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>

namespace sperr_lc {

using byte = unsigned char;
static constexpr int CS = 1024 * 16;  // 16 KiB chunk; must be multiple of 8

}  // namespace sperr_lc

using sperr_lc::byte;
using sperr_lc::CS;

#include "components/h_BIT_4.h"
#include "components/h_RZE_4.h"
#include "components/h_RZE_1.h"

namespace sperr_lc {

// Output layout (mirrors LC-standalone):
//   [8 bytes : long long insize]
//   [2 bytes × chunks : per-chunk compressed size (uint16, with chunks = ceil(insize/CS))]
//   [chunk data, concatenated]
// Returns total encoded byte length.
static inline long long encode(const byte* input,
                               long long insize,
                               byte*       output,
                               long long   output_cap)
{
  const long long chunks = (insize + CS - 1) / CS;
  const long long header_bytes = sizeof(long long) + chunks * sizeof(unsigned short);
  // Worst case: each chunk emits CS raw bytes. Caller must size output >= header + insize.
  if (output_cap < header_bytes + insize) return -1;

  auto* head_out = reinterpret_cast<long long*>(output);
  auto* size_out = reinterpret_cast<unsigned short*>(&head_out[1]);
  byte* data_out = reinterpret_cast<byte*>(&size_out[chunks]);

  // Encode each chunk independently. Carry holds the byte offset within data_out
  // for each chunk's *end*; serial pass after the parallel kernel can be cleaner
  // but we follow the standalone's pattern for determinism.
  // For simplicity (and because per-chunk dependencies are only on prior chunk's
  // offset), use sequential carry computation. (OMP version possible later.)
  long long offs = 0;
  for (long long chunkID = 0; chunkID < chunks; ++chunkID) {
    long long chunk1_buf[CS / sizeof(long long)];
    long long chunk2_buf[CS / sizeof(long long)];
    byte* in_buf  = reinterpret_cast<byte*>(chunk1_buf);
    byte* out_buf = reinterpret_cast<byte*>(chunk2_buf);

    const long long base = chunkID * CS;
    const int osize = static_cast<int>(std::min<long long>(CS, insize - base));
    std::memcpy(out_buf, &input[base], static_cast<size_t>(osize));

    int csize = osize;
    bool good = true;
    if (good) { std::swap(in_buf, out_buf); good = h_BIT_4(csize, in_buf, out_buf); }
    if (good) { std::swap(in_buf, out_buf); good = h_RZE_4(csize, in_buf, out_buf); }
    if (good) { std::swap(in_buf, out_buf); good = h_RZE_1(csize, in_buf, out_buf); }

    if (good && csize < osize) {
      size_out[chunkID] = static_cast<unsigned short>(csize);
      std::memcpy(&data_out[offs], out_buf, static_cast<size_t>(csize));
      offs += csize;
    } else {
      size_out[chunkID] = static_cast<unsigned short>(osize);
      std::memcpy(&data_out[offs], &input[base], static_cast<size_t>(osize));
      offs += osize;
    }
  }

  head_out[0] = insize;
  return header_bytes + offs;
}

// Decode the LC stream. Returns the decoded byte length (== original insize).
// Caller must size output >= insize (read from the stream header).
static inline long long decode(const byte* input,
                               long long   /*insize_compressed*/,
                               byte*       output,
                               long long   output_cap)
{
  const auto* head_in = reinterpret_cast<const long long*>(input);
  const long long outsize = head_in[0];
  if (outsize > output_cap) return -1;

  const long long chunks = (outsize + CS - 1) / CS;
  const auto* size_in = reinterpret_cast<const unsigned short*>(&head_in[1]);
  const byte* data_in = reinterpret_cast<const byte*>(&size_in[chunks]);

  long long pfs = 0;
  for (long long chunkID = 0; chunkID < chunks; ++chunkID) {
    long long chunk1_buf[CS / sizeof(long long)];
    long long chunk2_buf[CS / sizeof(long long)];
    byte* in_buf  = reinterpret_cast<byte*>(chunk1_buf);
    byte* out_buf = reinterpret_cast<byte*>(chunk2_buf);

    const long long base = chunkID * CS;
    const int osize = static_cast<int>(std::min<long long>(CS, outsize - base));
    int csize = size_in[chunkID];

    if (csize == osize) {
      std::memcpy(&output[base], &data_in[pfs], static_cast<size_t>(osize));
    } else {
      std::memcpy(out_buf, &data_in[pfs], static_cast<size_t>(csize));
      std::swap(in_buf, out_buf); h_iRZE_1(csize, in_buf, out_buf);
      std::swap(in_buf, out_buf); h_iRZE_4(csize, in_buf, out_buf);
      std::swap(in_buf, out_buf); h_iBIT_4(csize, in_buf, out_buf);
      if (csize != osize) return -1;
      std::memcpy(&output[base], out_buf, static_cast<size_t>(osize));
    }
    pfs += size_in[chunkID];
  }
  return outsize;
}

// Worst-case output capacity for encode(): chunks header + raw insize bytes.
static inline long long max_encoded_size(long long insize) {
  const long long chunks = (insize + CS - 1) / CS;
  return static_cast<long long>(sizeof(long long))
       + chunks * static_cast<long long>(sizeof(unsigned short))
       + insize;
}

}  // namespace sperr_lc

#endif  // SPERR_LC_PIPELINE_HPP
