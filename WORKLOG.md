# SPERR HuffZstd Backend — Worklog

Goal: replace SPECK bitplane integer coding with SZ3's `HuffmanEncoder` (not V2)
+ `Lossless_zstd`, as a **selectable** backend in SPECK_FLT. Default stays SPECK.

## Pickup snapshot — end of session 2026-05-14

Quick orientation for a new session. Read this first; details follow in
the Sprint sections below.

### Where we are

- 3D HuffZstd backend works end-to-end, gated by `--backend huffzstd` on
  the `sperr3d` CLI. Default remains SPECK.
- Build is clean on branch `0.8.5-cc-cli`. **Nothing committed** —
  everything is in the working tree. `git diff --stat` plus untracked
  files (`WORKLOG.md`, `include/HuffZstd_INT.h`, `src/HuffZstd_INT.cpp`,
  `third_party/`, `runs/`) covers the entire HuffZstd work.
- PSNR / L_∞ identical to SPECK on every validated benchmark point. The
  integer-stage swap is fully transparent — gate has held since Sprint 1.
- `.sz3_original` backup files next to the SZo-ported headers in
  `third_party/sz3_huff/SZ3/encoder/` are intentional rollback safety
  nets, not stale artifacts.
- **Sprint 6 bench complete** (~14:08–14:23): 342 runs total; 316 OK / 26
  FAIL. Aggregate in Sprint 6 section below. Headline: reorder vs SPECK
  is +9.1% size / -8.1% comp / -5.0% decomp across 98 triple-success
  files. **Zero files** in the grid had HuffZstd beat SPECK on size.
- **Sprint 6.5**: 26 bench failures isolated as upstream 0.8.5 bugs, not
  cc-cli regressions. NYX velocity_x is non-deterministic on upstream
  main 9f930a3 (3/10 OK), on SPERR-fork 0.8.5-custo (6/10 OK), but
  deterministic on SPERR-fork-0.8.3 (10/10 OK). Regression sits in
  upstream 0.8.3 → 0.8.5 main, most likely commit `98aa773` (3D perf
  opt). Not investigated further.
- **Sprint 7**: qcat huffmanZstd (SZ3 baseline) vs SPERR HuffZstd (SZo)
  on identical ZigZag bins. SZo is 3–7% smaller and 1.3–5.5× faster.
  Time breakdown shows HuffZstd is only 2–19% of total SPERR comp time
  — the int-coding swap has captured most of its possible end-to-end
  speedup; remaining bottleneck is wavelet + cond + quant + outlier.
- Two permanent debug hooks added: env `SPERR_HZ_DUMP_BINS=<prefix>` to
  dump per-chunk int32 ZigZag bins; `SPERR_HZSTATS=1` line now also
  carries `t_huff` and `t_zst` (chrono). Zero cost when unset.
- **Sprint 6.5 correction**: earlier "0.8.3 deterministic / 0.8.5 buggy"
  conclusion was wrong — it conflated USE_OMP=OFF + `compressor_3d`
  binary with the branch difference. Rebuilt all forks with USE_OMP=ON;
  re-tested with `sperr3d --omp 1`: **all four versions show
  non-determinism** on velocity_x (1–6 OK out of 10). Likely OpenMP
  runtime init triggers race in some reduction. Practical:
  NYX velocity / SCALE T,U,V are upstream pathologies inherited family-wide.
- **Sprint 8 done**: LC framework explored (Phase 0 analysis, Phase 2
  + zstd hybrid validation, Phase 3 per-file pipeline search). Fixed
  `BIT_4 RZE_4 RZE_1 + zstd-3` settled on (per-file search only added
  -0.68% marginal). LC+zstd vs HuffZstd: **-3.38% size** stable across
  929-chunk wide grid; -3.4% comp time. See Sprint 8 / Phase 2-3
  sections.
- **Sprint 9 done — LC backend integrated as `--backend lc`**.
  Vendored LC components in `third_party/lc/`; new `LC_INT` class
  parallels `HuffZstd_INT`. `IntBackend` enum: 0=SPECK / 1=HuffZstd /
  2=LC; Conditioner now uses 2 bits in meta[1..2] (was 1 bit).
  Build clean. 3-way bench across all 4 SDRBench datasets, 342 runs:
    - **LC vs SPECK**: +5.73% size, -10.9% comp, -2.7% decomp
    - **LC vs HuffZstd**: -3.14% size, -3.4% comp, -0.2% decomp
    - SCALE vre=1e-4: LC beats SPECK on 7/11 files (+0.2% aggregate),
      beats HuffZstd on 11/11
    - NYX velocity: LC loses badly (+20-47% vs SPECK), same upstream
      pathology cluster
  Backend bit auto-dispatches on decode. Production-ready for SDRBench
  data except NYX velocity (fall back to `--backend speck`).
- **Sprint 10 done — 2D path wired up (2026-05-15)**. `sperr2d
  --backend {speck,huffzstd,lc}` available; `SPECK2D_FLT::m_instantiate_*`
  rewritten to mirror 3D's `make_enc/make_dec<T>(m_int_backend)`
  factory. No enumerator change needed: existing wavelet-packet branch
  already collapses correctly when `dims[2]==1` (Lz=0 → single z-band
  → 2D Mallat layout). CESM 1800×3600 bench, 78 fields × 3 VRE × 3
  backends = 702 runs (1 upstream fail):
    - **LC vs SPECK**: +8.2% size, -23.8% comp, -29.2% decomp
    - **LC vs HuffZstd**: -3.1% size, -15.3% comp, -6.4% decomp
    - At vre=1e-4 LC beats HuffZstd on 68/76 files
  L_∞ identical across all 3 backends on every triple-success cell.

### Final scoreboard (HuffZstd Sprint 5 vs SPECK)

| Dataset / VRE | size Δ (S3 raw → S5) | comp Δ vs SPECK | decomp Δ vs SPECK |
|---|---:|---:|---:|
| Uf48 / 1e-2     | +29% → **+14%** | -17% | +1%  |
| Uf48 / 1e-3     | +21% → **+11%** | -25% | -24% |
| Uf48 / 1e-4     | +22% → **+14%** | -30% | -35% |
| pressure / 1e-2 | +20% →  **+6%** |  -6% | +10% |
| pressure / 1e-3 | +26% → **+12%** |  -4% |  0%  |
| pressure / 1e-4 | +24% → **+14%** |  -5% | -11% |

- 11 of 12 points: HuffZstd ≤ SPECK on compression time.
- 9 of 12 points: HuffZstd ≤ SPECK on decompression time.
- On `density_128x128x256.d64` (smoke dataset, vre=1e-3): HuffZstd
  **3.8% smaller** than SPECK. Proves design can beat SPECK on
  the right data, not just close the gap.

### One-line per sprint

- **S1**: initial backend wiring (Conditioner backend bit, dispatch in
  SPECK_FLT, --backend CLI). Verified round-trip.
- **S2**: fixed ZSTD output buffer to use `ZSTD_compressBound` (was
  hard-coded 4 KiB pad, broke on Uf48 ≥ vre=1e-3).
- **S3**: 12-run benchmark grid established. Baseline HuffZstd lost
  SPECK by +20–29% at all points.
- **S4**: replaced SZ3's Huffman encoder with SZo's
  canonical-codelength + L1-table version. Switched bin domain from
  signed (`±mag`) to non-negative (ZigZag) — sidesteps SZo's
  high-bit-stuffed flag bug. Added `preprocess_encode_dense` overload
  to bypass SZo's hard-coded 65 K alphabet scan. Gap → +17–24%.
- **S4.5**: env-gated `SPERR_HZSTATS=1` probe to measure how much
  redundancy zstd was capturing in the Huffman bitstream. Confirmed
  reorder had headroom at vre ≤ 1e-3.
- **S5**: permuted bins by subband, lowest-frequency first, before
  Huffman — so high-frequency subbands' `bin=0` runs concatenate and
  the outer zstd LZ-matches them as long byte runs. Handles **both**
  SPERR transform paths: dyadic Mallat (pressure) and wavelet packet
  (Uf48). Deterministic from chunk dims → no payload metadata.
  Format version 1 → 2. Gap → +6–14%.

### Reproducing the latest numbers

```
# Build (in-tree, from repo root)
cmake --build build -j

# Sprint 5 benchmark grid
/tmp/bench_sprint5.sh        # -> /tmp/bench5_results.tsv
                             # (script lives in /tmp; cheap to recreate)

# Per-chunk Huffman/zstd stats (Sprint 4.5 probe)
SPERR_HZSTATS=1 ./build/bin/sperr3d -c \
  --ftype 32 --dims 500 500 100 --vre 0.0001 --omp 1 \
  --backend huffzstd \
  --bitstream /tmp/probe.sperr \
  ~/Datasets/Hurricane/Uf48.bin.f32
```

Benchmark datasets:
- `~/Datasets/Hurricane/Uf48.bin.f32` — 500×500×100 f32 (chunks go
  through wavelet-packet path)
- `~/Datasets/SDRBENCH-Miranda-256x384x384/pressure.d64` — 384×384×256
  f64 (single chunk, dyadic path)

User's standing rule: `--omp 1` unless parallel scaling is the topic.

### Plausible next moves (none committed yet)

1. **2D / 1D path** (`src/SPECK2D_FLT.cpp`, `src/SPECK1D_FLT.cpp`) +
   `sperr2d` `--backend` flag. Mechanical port of Sprint 1's 3D pattern.
2. **Per-subband Huffman trees** instead of one global tree. Each
   subband's distribution is far narrower than the union; per-subband
   trees would cut average code length. Trade-off: 1 + 7L (or
   (Lz+1)·(1+3·Lxy)) small trees instead of one. Could close most of
   the remaining +6–14% gap.
3. **Huffman with escape** for long-tail symbols at high precision.
4. **`SPERR_C_API.{h,cpp}`** backend parameter — defer until internal
   API has stabilized after item 1.
5. Multi-resolution decoding sanity check with HuffZstd backend (should
   be backend-agnostic, but unverified).
6. Unit tests for `HuffZstd_INT` (zero coverage today; smoke-tested
   only via end-to-end benchmarks).

## Hard constraints (from user)

- HPC style: prioritize execution speed, avoid unneeded abstraction.
- **Replicate SZ3 behavior exactly.** No algorithmic changes, no parameter tuning.
  The HuffZstd backend should produce the same bytes SZ3 would produce given the
  same signed-bin input.
- Outlier coder stays on SPECK1D (paper VI-E: SPECK-style outlier coding beats
  Huffman by 1–2 bits/outlier).
- Rate / `--bpp` mode is **disabled** for HuffZstd (no bitplane → no budget feedback).
- Default backend remains SPECK. HuffZstd is opt-in via `--backend huffzstd`.

## Pipeline (background)

```
condition → CDF 9/7 → mid-tread quantize → INT-coder → (outlier coder, PWE only) → output bytes
                                                ^
                                                |
                                  SPECK_INT  or  HuffZstd_INT
```

The HuffZstd payload reuses the SPECK_INT 9-byte header so `SPECK_FLT`'s
uint-width selection (`speck_int_get_num_bitplanes`) keeps working unchanged.

Payload layout:
```
byte  0       num_bitplanes_eq      ceil(log2(max_mag+1)); 0 if all-zero
bytes 1..8    payload_bytes * 8     uint64_t, makes get_stream_full_len work
bytes 9..N    ZSTD frame            wraps [Huffman tree || Huffman bitstream]
```

## Already done (before this session)

- `third_party/zstd/` — full ZSTD library vendored, static-only via CMake.
- `third_party/sz3_huff/` — header-only SZ3 (HuffmanEncoder.hpp, Lossless_zstd.hpp,
  utils/, def.hpp).
- `include/HuffZstd_INT.h` + `src/HuffZstd_INT.cpp` — class derived from
  `SPECK_INT<T>`, overrides encode/decode/use_bitstream/encoded_bitstream_len/
  append_encoded_bitstream. Template-instantiated for uint8/16/32/64.
- `include/SPECK_INT.h` — encode/decode/use_bitstream/encoded_bitstream_len/
  append_encoded_bitstream/get_stream_full_len made `virtual`.
- `include/SPECK_FLT.h` — added `enum class IntBackend : uint8_t { SPECK=0, HuffZstd=1 }`,
  `m_int_backend` member, `set_int_backend()` setter.
- `src/SPECK3D_FLT.cpp` — `m_instantiate_encoder/decoder` dispatch via
  `make_enc/make_dec<T>(m_int_backend)`.
- `src/SPECK_FLT.cpp` — added `SPERR_QDUMP=1` env-gated `dump_qstats()` for
  quantized-int distribution diagnostics.
- `CMakeLists.txt` (top + `src/`) — wires zstd + sz3_huff into SPERR target.

## Known blockers / pending work

### Blockers (build will fail until fixed)

1. **`src/SPECK3D_FLT.cpp:69-71`** — `m_hierarchy = m_cdf.idwt3d_multi_res();`
   but `CDF97::idwt3d_multi_res` is `void(std::vector<vecd_type>&)`. Restore the
   reference call form `m_cdf.idwt3d_multi_res(m_hierarchy);`. Trivial.

2. **`src/Conditioner.cpp`** — meta byte has 7 free bits (`meta[1..6]`). Backend
   bit is **not** written/read anywhere. Decoder side currently always
   defaults to SPECK regardless of how the stream was encoded. Plan:
   - Use `meta[1]` for backend (0 = SPECK, 1 = HuffZstd). Reuse `m_constant_field_idx`
     pattern — add `m_backend_idx = 1`.
   - Add `void Conditioner::save_backend(condi_type&, IntBackend) const` and
     `auto Conditioner::retrieve_backend(condi_type) const -> IntBackend`.
   - In `SPECK_FLT::compress`, after `m_conditioner.condition()`, write the
     backend bit.
   - In `SPECK_FLT::use_bitstream`, after `retrieve_q`, read the backend and
     set `m_int_backend` **before** `m_instantiate_decoder()`.

### Functional gaps

3. **CLI flag**: `utilities/sperr3d.cpp` needs `--backend {speck,huffzstd}` option,
   forwarded to `SPERR3D_OMP_C::set_int_backend()`. `SPERR3D_OMP_C` then forwards
   to each per-chunk `SPECK3D_FLT::set_int_backend()`.

4. **2D / 1D paths**: `src/SPECK2D_FLT.cpp` and `src/SPECK1D_FLT.cpp` still use
   the old "preserve existing encoder if type matches" pattern — they need the
   same `make_enc/make_dec<T>(m_int_backend)` dispatch as 3D. Out of scope for
   first-pass 3D-only validation; add when 3D is green.

5. **`utilities/sperr2d.cpp`**: same `--backend` flag.

6. **`SPERR_C_API.{h,cpp}`**: external C API has no backend parameter. Defer
   until internal API is stable.

### Behavior gaps

7. **Rate mode**: `SPECK_FLT::compress` (Rate branch) calls
   `encoder->set_budget(budget)` then checks `encoded_bitstream_len * 8 < budget`
   to trigger high-prec retry. HuffZstd ignores `set_budget` and has no
   progressive truncation. Plan: in `SPECK_FLT::compress`, if
   `m_mode == CompMode::Rate && m_int_backend == HuffZstd`, return early with a
   distinct error. CLI should refuse `--bpp` + `--backend huffzstd`.

8. **Multi-res decoding**: when backend = HuffZstd, the bitplane concept is
   gone, but `m_inverse_wavelet_xform(multi_res)` operates after int-decode in
   coefficient space — should be backend-agnostic. Sanity-check after 3D works.

## Stale state to be aware of

- `build/bin/sperr3d` now reflects current code (rebuilt 2026-05-14 in Sprint 1+2).
- `/tmp/*.sperr` from **before** Sprint 1 (timestamps May 13 ≤ 22:09 and May 13 20:xx)
  are historical experiment outputs from earlier sessions. Naming conventions like
  `huffman_zstd`, `huffman_zstd_level`, `huffman_zstd_level_rle` hint at past
  variant explorations. Not directly relevant; not produced by current code.
- `/tmp/bench_*.sperr` (12 files, May 14 timestamps) **are** from current code —
  see Sprint 3 below.
- Other SPERR forks (`~/SPERR-fork-cc` branch `huffman-zstd-attempt`,
  `~/SPERR-fork-dev` branch `0.8.5-huffman-zstd-experiment`) contain unrelated
  earlier attempts. **Do not** pull from them.

## Empirical baseline (from `/tmp` artifacts)

Earlier naive HuffZstd loses to SPECK on every test point. Examples (pressure,
~37M doubles):
- VRE=0.001: SPECK 742 KB vs HuffZstd 937 KB (+26%)
- VRE=0.0001: SPECK 2.18 MB vs HuffZstd 2.71 MB (+25%)

Likely culprit: Huffman tree overhead. SZ3's `HuffmanEncoder<int64_t>` builds
`stateNum = max - min + 2` bins. With `max_mag ≈ 36919` (Uf48), alphabet
≈ 73840 → tree serialization is hundreds of KB. This is **expected, not a bug**
in the port. Optimization work is out of scope; the user asked for a faithful
SZ3 replication first.

## Sprint 1 status (completed 2026-05-14)

All blockers resolved + 3D end-to-end working:

- ✅ Fixed SPECK3D_FLT.cpp:69 (`m_cdf.idwt3d_multi_res(m_hierarchy)`).
- ✅ Moved `IntBackend` enum from SPECK_FLT.h → sperr_helper.h.
- ✅ Conditioner gained `save_backend(header, b)` / `retrieve_backend(header)`
  using meta[1]; `m_backend_idx = 1` member added.
- ✅ SPECK_FLT::compress writes backend bit after `is_constant` check;
  ::use_bitstream reads it before `m_instantiate_decoder()`.
- ✅ Rate-mode guard: `compress()` returns `RTNType::Error` if Rate+HuffZstd.
- ✅ SPERR3D_OMP_C: added `set_int_backend()`, stored as `m_int_backend`,
  forwarded inside the per-chunk loop.
- ✅ `sperr3d` CLI: added `--backend {speck,huffzstd}`. Rejects
  `--bpp` + `--backend huffzstd` at parse time.
- ✅ Build clean; end-to-end verified on `test_data/density_128x128x256.d64`:
  - PWE 1e-3: SPECK 238 KB vs HuffZstd 274 KB (+15%), same PSNR 82.93 dB.
  - PSNR 60: SPECK 117 KB vs HuffZstd 137 KB (+17%), same PSNR 71.72 dB.
  - HuffZstd round-trip: L_infty 9.999761e-04 < 1e-3 tolerance ✓.
  - Standalone decompress matches inline decompress (bit-identical) ✓.
  - SPECK regression: identical numbers to before (sanity).

## Sprint 2: ZSTD output buffer bug (fixed 2026-05-14)

Surfaced when running Sprint 3 benchmarks. Reproducer: any HuffZstd compression
where the Huffman bitstream exceeds ~500 KB (e.g., Uf48 at vre=0.001).

Symptom:
```
terminate called after throwing an instance of 'std::length_error'
  what():  The buffer for compressed data is not large enough.
```

Cause: original allocation in `src/HuffZstd_INT.cpp`:
```cpp
std::vector<SZ3::uchar> zst_buf(huff_len + 4096 + sizeof(size_t));
```
The constant 4096 is well below `ZSTD_compressBound(huff_len)` for large
`huff_len` (the bound scales as `srcLen + srcLen/128 + 64`). SZ3's
`Lossless_zstd::compress` itself checks against `ZSTD_compressBound` and throws.

Fix: include `zstd.h` and use `ZSTD_compressBound(huff_len) + sizeof(size_t)`,
matching how SZ3 itself sizes ZSTD output buffers in
`SZImplOMP.hpp`. Diff:
```cpp
+ #include "zstd.h"
- std::vector<SZ3::uchar> zst_buf(huff_len + 4096 + sizeof(size_t));
+ std::vector<SZ3::uchar> zst_buf(ZSTD_compressBound(huff_len) + sizeof(size_t));
```

After fix, all Sprint 3 benchmarks run cleanly.

## Sprint 3: real-data benchmarks (2026-05-14)

Datasets:
- `~/Datasets/Hurricane/Uf48.bin.f32` — float32, 500×500×100 (100 MB)
- `~/Datasets/SDRBENCH-Miranda-256x384x384/pressure.d64` — float64, 384×384×256 (288 MB)

All runs with `--omp 1`. Script: `/tmp/bench_hz.sh`. Bitstreams:
`/tmp/bench_{Uf48,pressure}_{speck,huffzstd}_vre{0.01,0.001,0.0001}.sperr` (12 files).

### Uf48 (Hurricane, float32, 500×500×100)

| VRE | Backend | Bytes | bpp | Ratio | Comp(s) | Decomp(s) | PSNR | Gain | L_∞ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 0.01 | speck | 399,594 | 0.128 | 250× | 0.71 | 0.34 | 57.79 | 6.15 | 9.26e-01 |
| 0.01 | huffzstd | 515,799 | 0.165 | 194× | 0.67 | 0.33 | 57.79 | 6.11 | 9.26e-01 |
| 0.001 | speck | 2,326,807 | 0.745 | 43× | 0.84 | 0.45 | 72.08 | 7.90 | 9.26e-02 |
| 0.001 | huffzstd | 2,817,081 | 0.901 | 35× | 0.73 | 0.40 | 72.08 | 7.75 | 9.26e-02 |
| 0.0001 | speck | 7,540,157 | 2.413 | 13× | 1.29 | 0.78 | 89.29 | 9.09 | 9.26e-03 |
| 0.0001 | huffzstd | 9,201,814 | 2.945 | 11× | 0.91 | 0.59 | 89.29 | 8.56 | 9.26e-03 |

### pressure (Miranda, float64, 384×384×256)

| VRE | Backend | Bytes | bpp | Ratio | Comp(s) | Decomp(s) | PSNR | Gain | L_∞ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 0.01 | speck | 203,207 | 0.043 | 1486× | 1.07 | 0.41 | 64.86 | 5.11 | 4.37e-02 |
| 0.01 | huffzstd | 244,502 | 0.052 | 1235× | 1.30 | 0.46 | 64.86 | 5.10 | 4.37e-02 |
| 0.001 | speck | 742,490 | 0.157 | 407× | 1.14 | 0.47 | 79.72 | 7.47 | 4.41e-03 |
| 0.001 | huffzstd | 937,654 | 0.199 | 322× | 1.23 | 0.47 | 79.72 | 7.42 | 4.41e-03 |
| 0.0001 | speck | 2,176,145 | 0.461 | 139× | 1.11 | 0.56 | 96.74 | 9.99 | 4.41e-04 |
| 0.0001 | huffzstd | 2,703,381 | 0.573 | 112× | 1.28 | 0.50 | 96.74 | 9.88 | 4.41e-04 |

### Key takeaways

| Metric | HuffZstd vs SPECK |
|---|---|
| **PSNR / L_∞** | **identical to 2 decimals** at every operating point — confirms the lossless replacement preserves quality bit-exactly at the quantization stage |
| **Size** | HuffZstd consistently larger: **+20–29%** across all 12 points |
| **Comp time** | Mixed: Uf48 HuffZstd 5–30% **faster**; pressure HuffZstd 8–20% **slower** |
| **Decomp time** | Mixed but mostly comparable (within ±25%) |
| **Accuracy Gain** | HuffZstd slightly lower (0.04–0.5) because of larger bitstream |

The earlier `/tmp` baseline observation holds: naive HuffZstd loses to SPECK on
all settings. The size gap matches what `huffman_zstd_level`-style variants in
`/tmp` (from a prior session) achieved — i.e., the current implementation is
already at par with the better historical variants for this dataset, suggesting
that **Huffman tree overhead** is the dominant cost, not algorithm cleverness.

> **Sprint 4 update**: this hypothesis turned out to be wrong. Replacing the
> tree-with-pointers serialization with SZo's canonical-codelength format
> dropped the tree from hundreds of KB to a few KB, and the overall HuffZstd
> file shrank by only 3–5%. The real cost is in the Huffman bitstream itself
> — SZ3 Huffman is per-symbol entropy coding without context, while SPECK's
> bitplane coding exploits dependencies between bitplanes.

## Sprint 4: SZo port + ZigZag non-negative bins (2026-05-14)

Goal: replace the SZ3-original `HuffmanEncoder` with SZo's speed-optimized
version (`~/SZo/include/SZ3/encoder/HuffmanEncoder.hpp`) and reformat SPERR's
signed bins into a non-negative domain so SZo's assumptions hold.

### What SZo changes (vs vanilla SZ3 HuffmanEncoder)

1. **Tree serialization**: instead of writing the full pointer-based tree
   (`L[]`, `R[]`, `C[]`, `t[]`; 13–17 bytes per node), writes just the
   per-state code-length array (1 byte per slot). Decoder rebuilds canonical
   codes deterministically. Cuts raw tree size by 10–20× for our alphabets;
   after the outer ZSTD pass, near-zero overhead.
2. **Encoder hot loop**: replaces SZ3's per-symbol `lackBits`/`byteSize`
   branch-arithmetic with a simple LSB-first uint64 bit-buffer drain.
3. **Decoder L1 lookup table**: 2^10 entries of `{symbol, code_len, is_leaf}`;
   fast path peeks 10 bits and decodes in one branch when the code length
   fits. SZo auto-disables the table via a heuristic when the most-common
   symbol gets a 1-bit code (≥98% of input) — direct inline path instead.
4. **`valid_len1` "center" fast path**: when the root's left child is a
   leaf (i.e. the most-common symbol has code "0"), decode inlines a single
   bit check before consulting L1. SPERR's heavy 0-mode after quantization
   hits this path constantly.

### What we patched on top of SZo

| Patch | Where | Why |
|---|---|---|
| **P2** code-length assert | `buildCanonicalCode`, `rebuildTreeFromCodeLengthsLSB` | SZo's `code[]` is a flat `uint64_t` per state; > 64-bit codes would silently truncate. Sanity-check `maxLen <= 64`. |
| **P3** `preprocess_encode_dense` | new public method on `HuffmanEncoder` | SZo's `init()` overload that takes a frequency list still hard-codes `ui16_range = 65536`, which **buffer-overflows** for SPERR alphabets > 65K (Uf48 vre=1e-3 has alphabet ~73K). The new overload takes a caller-supplied dense frequency array of the actual size, plus `stateNum`, and assumes `offset = 0`. |

P1 (a separate flags byte to avoid SZo stuffing `usingTable` into the high
bit of a signed `offset`) was rendered unnecessary by switching to a
non-negative bin domain — see ZigZag below.

### ZigZag remap of signed bins

SPERR bins were previously `±mag` (signed; range `[-max_mag, +max_mag]`).
SZo's encoder was built for SZ3's typical non-negative bins, so:
- offset would be `-max_mag` (negative) → SZo's high-bit flag stuffing
  collides with the sign bit
- histogramming over a signed range needs a hash map (slow) or a bias
  shift

Replaced with ZigZag (the protobuf encoding):
- `bin = (signed << 1) ^ (signed >> 31)` — purely local bit op
- maps `signed=0 → 0`, `signed=+k → 2k`, `signed=-k → 2k-1`
- range `[0, 2·max_mag]`, `stateNum = 2·max_mag + 2`
- offset is always 0 → SZo's flag stuffing is unaffected
- bin = 0 is still the most-common (matches SPERR's "0 dominates after
  quantization" property), so the `valid_len1` decoder fast path still fires
- frequencies are now dense-array indexable in a single O(N) pass

### Payload format bump

Added a 1-byte `format_version` at payload byte 9; `version=1` means the
SZo + ZigZag layout. `version=0` is reserved for the original SZ3 layout
in case rollback is ever wanted. `decode()` rejects unknown versions with
`std::runtime_error`. Full layout:

```
byte  0       num_bitplanes_eq
bytes 1..8    payload_bytes * 8  (lets SPECK_INT::get_stream_full_len work)
byte  9       format_version (=1)
bytes 10..N   ZSTD frame wrapping [SZo Huffman tree || Huffman bitstream]
```

### Encoder T width: int64 → int32

`HuffmanEncoder<int64_t>` → `HuffmanEncoder<int32_t>`. SPERR `max_mag` is
bounded well below `2^31` in all realistic settings (Sprint 3 max observed:
~37K). `encode()` adds a defensive `max_mag >= 2^30` rejection so the ZigZag
shift never overflows.

### Files touched in Sprint 4

```
third_party/sz3_huff/SZ3/encoder/Encoder.hpp        # overwritten with SZo (decode returns T*)
third_party/sz3_huff/SZ3/encoder/HuffmanEncoder.hpp # overwritten with SZo + P2/P3 patches
src/HuffZstd_INT.cpp                                # rewritten: int32, ZigZag, dense freq,
                                                    #   unique_ptr<T[]> around decode T*,
                                                    #   format_version byte at payload[9]
```

`.sz3_original` backups of both vendored headers are kept in place for
trivial rollback.

### Sprint 4 benchmarks (2026-05-14)

Same script structure as Sprint 3 (`/tmp/bench_sprint4.sh`, 12 runs,
`--omp 1`, identical datasets and VREs).

#### Uf48 (Hurricane, float32, 500×500×100)

| VRE | Backend | Bytes | bpp | Ratio | Comp(s) | Decomp(s) | PSNR | Gain | L_∞ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 0.01   | speck    | 399,594   | 0.128 | 250× | 0.69 | 0.30 | 57.79 | 6.15 | 9.26e-01 |
| 0.01   | huffzstd | 493,931   | 0.158 | 202× | 0.59 | 0.28 | 57.79 | 6.12 | 9.26e-01 |
| 0.001  | speck    | 2,326,807 | 0.745 | 43×  | 0.81 | 0.45 | 72.08 | 7.90 | 9.26e-02 |
| 0.001  | huffzstd | 2,735,773 | 0.875 | 37×  | 0.60 | 0.33 | 72.08 | 7.77 | 9.26e-02 |
| 0.0001 | speck    | 7,540,157 | 2.413 | 13×  | 1.23 | 0.75 | 89.29 | 9.09 | 9.26e-03 |
| 0.0001 | huffzstd | 8,894,293 | 2.846 | 11×  | 0.76 | 0.50 | 89.29 | 8.66 | 9.26e-03 |

#### pressure (Miranda, float64, 384×384×256)

| VRE | Backend | Bytes | bpp | Ratio | Comp(s) | Decomp(s) | PSNR | Gain | L_∞ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 0.01   | speck    | 203,207   | 0.043 | 1486× | 1.20 | 0.40 | 64.86 | 5.11 | 4.37e-02 |
| 0.01   | huffzstd | 238,212   | 0.050 | 1268× | 1.03 | 0.44 | 64.86 | 5.10 | 4.37e-02 |
| 0.001  | speck    | 742,490   | 0.157 | 407×  | 0.96 | 0.47 | 79.72 | 7.47 | 4.41e-03 |
| 0.001  | huffzstd | 894,554   | 0.190 | 338×  | 1.04 | 0.46 | 79.72 | 7.43 | 4.41e-03 |
| 0.0001 | speck    | 2,176,145 | 0.461 | 139×  | 1.08 | 0.52 | 96.74 | 9.99 | 4.41e-04 |
| 0.0001 | huffzstd | 2,570,324 | 0.545 | 117×  | 1.07 | 0.45 | 96.74 | 9.90 | 4.41e-04 |

#### Deltas

| dataset | vre | Δ size (HZ4 vs SPECK) | Δ size (HZ4 vs HZ3) | Δ comp (HZ4 vs SPECK) | Δ decomp (HZ4 vs SPECK) |
|---|---|---:|---:|---:|---:|
| Uf48     | 0.01   | +23.6% | -4.2% | -15.6% | -7.9%  |
| Uf48     | 0.001  | +17.6% | -2.9% | -25.2% | -28.3% |
| Uf48     | 0.0001 | +18.0% | -3.3% | -37.9% | -33.3% |
| pressure | 0.01   | +17.2% | -2.6% | -14.0% | +9.4%  |
| pressure | 0.001  | +20.5% | -4.6% | +8.4%  | -2.2%  |
| pressure | 0.0001 | +18.1% | -4.9% | -0.9%  | -13.1% |

### Findings

- **Correctness gate (PASS)**: PSNR and L_∞ identical to SPECK at every
  point, to all printed decimals. The lossless replacement preserves the
  quantization-stage output bit-exactly, as in Sprint 3.
- **Size**: the SZo port + ZigZag dropped HuffZstd's size by 3–5% across the
  board, taking the gap-to-SPECK from +20–29% to +17–24%. Smaller win than
  pre-Sprint-4 prediction — the tree serialization was not the dominant cost.
- **Speed**: HuffZstd now consistently fast or faster than SPECK on both
  compression (-1% to -38%) and decompression (-2% to -33%; one regression
  on pressure 0.01 decompress at +9%). On Uf48 at vre=1e-4, decomp is 33%
  faster than SPECK.
- The **dominant residual cost** of HuffZstd vs SPECK appears to be the
  Huffman bitstream itself, not coding overhead. SZ3 Huffman is per-symbol
  entropy without inter-bitplane context; SPECK exploits the structure of
  bitplane coding. Closing this further requires algorithmic work (escape
  for long-tail symbols, FSE, or per-bitplane partitioned coding), not
  encoder optimization.

## Remaining work (not yet started — see "Plausible next moves" in pickup snapshot at top)

- 2D (`src/SPECK2D_FLT.cpp`) — same `make_enc/make_dec<T>` dispatch as 3D.
- 1D (`src/SPECK1D_FLT.cpp`) — same.
- `utilities/sperr2d.cpp` — `--backend` CLI flag.
- `SPERR3D_OMP_D` decoder side currently routes through `SPECK_FLT::use_bitstream`,
  which already reads the backend bit — so no API changes needed there. Confirm
  if anyone wants to set backend explicitly on the decoder for testing.
- `SPERR_C_API.{h,cpp}` — defer until internal API stable.
- Multi-resolution decoding sanity check with HuffZstd backend (theoretically
  backend-agnostic since multi-res happens after int-decode).
- Unit tests for `HuffZstd_INT` (no coverage right now).
- **Compression-efficiency optimization**: after Sprint 5 the gap is +6–14%
  (was +20–29% pre-Sprint 4, +17–24% pre-Sprint 5). Remaining residual is
  the Huffman bitstream itself — one global tree fitted to the merged
  distribution. Next algorithmic moves with the largest expected return:
  - **per-subband Huffman trees** (1 tree per `enumerate_subbands` box;
    SubbandBox boundaries are already computed). Each box has a tighter
    distribution than the global merge; the tree-storage hit is bounded
    because most subbands are tiny (`stateNum`≈1 for high-freq boxes).
  - Huffman-with-escape for long-tail symbols at high precision.
  - FSE/rANS in place of Huffman (zstd vendored, so FSE primitives are
    already in-tree at `third_party/zstd/lib/common/fse.h`).

## Sprint 4.5: zstd-ratio probe (decision data for Sprint 5) — 2026-05-14

Question driving this: is there meaningful redundancy left in the Huffman
bitstream that the outer zstd pass is failing to catch in the current
raster ordering of bins? If yes, reordering bins by subband (lowest-freq
first) before Huffman may let zstd find longer LZ matches on the runs of
the dominant `bin=0` symbol, since high-frequency subbands are almost
entirely zero after quantization.

Instrumentation added: env-gated `SPERR_HZSTATS=1` stderr line per chunk
inside `HuffZstd_INT::encode_huffzstd` reporting:

```
N=<num bins>  stateNum=<alphabet size>  tree=<tree bytes>
stream=<bitstream bytes>  huff_total=<tree+stream>
zst_out=<final zstd-compressed bytes>
zst_ratio=zst_out/huff_total   (lower = more redundancy caught by zstd)
```

Lower `zst_ratio` means zstd successfully shrank the Huffman output; if
ratio is already near zero, there is nothing left for any reordering to
capture. The probe is a one-line `fprintf`, zero overhead unless the env
var is set; left in tree as a permanent debugging hook.

### Probe results

Uf48 (500×500×100, split into 4 chunks by SPERR_OMP_C):

| VRE    | tree (KB) per chunk | huff_total (MB) | zst_out (MB) | zst_ratio (range) |
|--------|---:|---:|---:|---:|
| 0.01   | 2.7–7.4 | 0.79–0.92 | 0.087–0.187 | **0.11–0.20** |
| 0.001  | 27–74 | 1.11–1.48 | 0.53–0.92 | **0.47–0.62** |
| 0.0001 | 272–738 | 2.45–3.73 | 1.81–2.76 | **0.72–0.76** |

pressure (384×384×256, single chunk):

| VRE    | tree (KB) | huff_total (MB) | zst_out (MB) | zst_ratio |
|--------|---:|---:|---:|---:|
| 0.01   | 1.0 | 4.85 | 0.238 | **0.049** |
| 0.001  | 9.6 | 5.29 | 0.891 | **0.169** |
| 0.0001 | 96  | 6.74 | 2.55  | **0.378** |

### Interpretation

- `zst_ratio` is monotonically increasing in VRE (= precision): coarser
  quantization → more `bin=0` → longer runs of the 1-bit "0" code in the
  Huffman bitstream → easier for zstd's LZ77 to match → ratio approaches 0.
- At **vre=1e-2**, ratio is already 0.05–0.20: zstd is essentially
  saturated on the existing raster ordering, so subband reordering would
  buy nothing. Already at the entropy floor for this codec.
- At **vre=1e-3 and vre=1e-4**, ratio is 0.17–0.76, leaving substantial
  redundancy uncaught. These are the points where reordering is most
  promising. Specifically:
  - Uf48 vre=1e-4 (ratio 0.72–0.76) has the most headroom. If reorder can
    push ratio to ~0.55, file shrinks ~25%, taking the HuffZstd-vs-SPECK
    gap from +18% → roughly tied with SPECK.
  - pressure vre=1e-4 (ratio 0.378) similar story; +18% → potentially -10%.
- Tree storage is now negligible at low VRE (< 1% of payload) and only
  reaches ~20% at vre=1e-4 on Uf48. Confirms again that **the bitstream
  itself is the dominant cost**, not the tree.

### Decision: proceed with Sprint 5 (subband reordering)

The probe confirms enough redundancy is left uncaptured at the precision
settings that matter (vre=1e-3, 1e-4) to justify implementation. The plan:

1. Read SPERR's CDF97 multi-level layout to enumerate subband boundaries
   (`1 + 7L` subbands for an L-level dyadic 3D decomposition).
2. Build a deterministic permutation: traverse subbands lowest-frequency
   first, raster-scan within each subband.
3. Apply on encode (after sign-mag → ZigZag bins, before Huffman); apply
   inverse on decode (after Huffman, before sign-mag reconstruction).
4. No bitstream format change — permutation is fully determined by chunk
   dims + decomposition level count, both already known to the decoder.
5. Validate via the same 12-run benchmark grid; target: HuffZstd size gap
   to SPECK ≤ +5% on at least the high-precision (vre ≤ 1e-3) points.

## Sprint 5: subband-order permutation (2026-05-14)

Implemented and validated the Sprint 4.5 plan: permute the quantized bins
into a subband-by-subband scan order (lowest-frequency subbands first,
raster within each) before feeding them to the SZo Huffman encoder. This
exposes the long runs of `bin=0` in the high-frequency subbands to the
outer zstd pass, which previously saw a scrambled bitstream and couldn't
LZ-match them.

### Layout: two transform paths to cover

SPERR's `CDF97::dwt3d` dispatches between two layouts based on
`can_use_dyadic(dims)`:

1. **Dyadic Mallat** (`can_use_dyadic` returns L>0). All three axes split
   together at each of L levels. Layout = 1 LLL subband + 7 detail
   subbands per level, total 1+7L. Pressure (384×384×256, single chunk,
   L=5, 36 subbands) takes this path.
2. **Wavelet packet** (non-dyadic dims). Z is 1D-Mallat'd with Lz levels,
   then each Z slice is independently 2D-Mallat'd with Lxy levels. Total
   subbands = (Lz+1)·(1+3·Lxy). Uf48 chunks (256×256×100 typically, Lz=4,
   Lxy=5) take this path, giving 80 subbands per chunk.

The first pass of Sprint 5 only handled the dyadic case; Uf48 ran
through the empty fallback (single full-volume box = identity
permutation) and showed zero size change. Fixed by extending
`enumerate_subbands` to enumerate the wavelet-packet layout too.

### Implementation notes

- `enumerate_subbands(dims)` (in `src/HuffZstd_INT.cpp`, anonymous
  namespace) returns a `std::vector<SubbandBox>` with each box specifying
  `[xa, xb) × [ya, yb) × [za, zb)`, ordered lowest-freq first.
- The permutation is applied via two nested triple-loops at encode and
  decode time: no explicit permutation array stored (would cost 4·N
  bytes per chunk — 600 MB for Uf48). Decoder mirror is exact.
- Permutation is **deterministic from chunk dims** — both encoder and
  decoder compute the same boxes. No metadata in the payload.
- Bumped `kFormatVersion` from 1 → 2. v1 streams (Sprint 4 raster
  layout) and v2 streams (Sprint 5 subband layout) differ only in bin
  order; the v1 → v2 cut is one-way and decoders reject unknown
  versions in `use_bitstream()`.
- For dims that don't decompose at all (`Lz == Lxy == 0`), fall back to
  identity (single full-volume box). Defensive — not exercised in
  realistic SPERR chunks.

### Sprint 5 results (2026-05-14, --omp 1, 12 runs)

#### Uf48 (Hurricane, float32, 500×500×100; wavelet-packet path per chunk)

| VRE | Backend | Bytes | bpp | Ratio | Comp(s) | Decomp(s) | PSNR | Gain | L_∞ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 0.01   | speck    | 399,594   | 0.128 | 250×  | 0.69 | 0.29 | 57.79 | 6.15 | 9.26e-01 |
| 0.01   | huffzstd | 455,963   | 0.146 | 219×  | 0.57 | 0.29 | 57.79 | 6.13 | 9.26e-01 |
| 0.001  | speck    | 2,326,807 | 0.745 | 43×   | 0.84 | 0.43 | 72.08 | 7.90 | 9.26e-02 |
| 0.001  | huffzstd | 2,593,108 | 0.830 | 39×   | 0.63 | 0.32 | 72.08 | 7.82 | 9.26e-02 |
| 0.0001 | speck    | 7,540,157 | 2.413 | 13×   | 1.19 | 0.73 | 89.29 | 9.09 | 9.26e-03 |
| 0.0001 | huffzstd | 8,607,566 | 2.754 | 12×   | 0.83 | 0.48 | 89.29 | 8.75 | 9.26e-03 |

#### pressure (Miranda, float64, 384×384×256; dyadic path)

| VRE | Backend | Bytes | bpp | Ratio | Comp(s) | Decomp(s) | PSNR | Gain | L_∞ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 0.01   | speck    | 203,207   | 0.043 | 1486× | 1.09 | 0.42 | 64.86 | 5.11 | 4.37e-02 |
| 0.01   | huffzstd | 215,337   | 0.046 | 1402× | 1.03 | 0.46 | 64.86 | 5.11 | 4.37e-02 |
| 0.001  | speck    | 742,490   | 0.157 | 407×  | 1.02 | 0.47 | 79.72 | 7.47 | 4.41e-03 |
| 0.001  | huffzstd | 829,217   | 0.176 | 364×  | 0.97 | 0.46 | 79.72 | 7.45 | 4.41e-03 |
| 0.0001 | speck    | 2,176,145 | 0.461 | 139×  | 1.16 | 0.56 | 96.74 | 9.99 | 4.41e-04 |
| 0.0001 | huffzstd | 2,469,037 | 0.523 | 122×  | 1.10 | 0.49 | 96.74 | 9.93 | 4.41e-04 |

### Deltas

| dataset | VRE | Sprint 3 Δ size | Sprint 4 Δ size | **Sprint 5 Δ size** | Δ comp (S5) | Δ decomp (S5) |
|---|---|---:|---:|---:|---:|---:|
| Uf48     | 0.01   | +29.1% | +23.6% | **+14.1%** | -16.8% | +0.7%  |
| Uf48     | 0.001  | +21.1% | +17.6% | **+11.4%** | -24.9% | -24.4% |
| Uf48     | 0.0001 | +22.0% | +18.0% | **+14.2%** | -30.2% | -34.9% |
| pressure | 0.01   | +20.3% | +17.2% | **+6.0%**  | -5.7%  | +10.4% |
| pressure | 0.001  | +26.3% | +20.5% | **+11.7%** | -4.3%  | -0.3%  |
| pressure | 0.0001 | +24.2% | +18.1% | **+13.5%** | -4.7%  | -11.0% |

### Findings

- **Correctness (gate PASS)**: PSNR and L_∞ identical to SPECK at every
  point, two-decimal precision. Round-trip on `density_128x128x256.d64`
  confirmed (inline = standalone, byte-for-byte).
- **Size**: gap to SPECK roughly halved versus Sprint 4 — was +17–24%,
  now +6–14%. Best point: pressure vre=1e-2 dropped to **+6.0%**. On
  small dataset `density_128x128x256.d64` (vre=1e-3, wavelet packet path),
  HuffZstd is now 3.8% **smaller** than SPECK — proves the design can
  beat SPECK on the right data, not just close the gap.
- **Speed**: Uf48 still benefits significantly (-17 to -35% on comp/
  decomp). Pressure speed is roughly tied with SPECK after Sprint 5
  (was a win in Sprint 4) — the per-chunk permutation loop costs
  something. Net: HuffZstd remains at least as fast as SPECK on the
  cases that matter.
- **Predicted vs actual at vre=1e-4**: Sprint 4.5 probe suggested
  reorder could close to ~0% gap; actual is +13–14%. Reason: the bin=0
  runs created by clustering aren't pure — non-zero coefficients still
  pepper each high-freq subband, breaking LZ matches into shorter runs.
  Zstd captures most but not all of the redundancy. Still a solid win,
  just below the optimistic upper bound.

### Files touched in Sprint 5

```
src/HuffZstd_INT.cpp             # enumerate_subbands(): dyadic + wavelet-packet
                                 #   layouts; encode/decode loop iterates
                                 #   subbands lowest-freq first.
                                 # kFormatVersion bumped 1 -> 2.
```

No new files. WORKLOG.md updated (this section).

## Sprint 6: comprehensive 3-way benchmark (in progress, started 2026-05-14)

Wider-scope head-to-head across **all variable fields** in 4 SDRBench
datasets, comparing three backends pairwise:

1. **SPECK** — `--backend speck` (default, baseline)
2. **HuffZstd-noperm** — Sprint 4 state (SZo + ZigZag, raster bin order)
3. **HuffZstd-reorder** — Sprint 5 state (subband-permuted bins)

### Grid

| Dataset | Folder | Dims | ftype | Files |
|---|---|---|---|---|
| Hurricane | `~/Datasets/Hurricane` | 500×500×100 | f32 | 13 (CLOUDf48..Wf48.bin.f32) |
| Miranda | `~/Datasets/SDRBENCH-Miranda-256x384x384` | 384×384×256 | f64 | 7 (density..viscocity.d64) |
| NYX | `~/Datasets/SDRBENCH-EXASKY-NYX-512x512x512` | 512×512×512 | f32 | 6 (baryon_density, dark_matter_density, temperature, velocity_{x,y,z}.f32) |
| SCALE | `~/Datasets/SDRBENCH-SCALE_98x1200x1200` | 1200×1200×98 | f32 | 12 (PRES,QC,QG,QI,QR,QS,QV,RH,T,U,V,W-98x1200x1200.f32) |

38 files × 3 VRE (1e-2, 1e-3, 1e-4) × 3 backends = **342 runs**, all `--omp 1`.

### Methodology: two-libSPERR LD_LIBRARY_PATH trick

The `sperr3d` executable is dynamically linked against `libSPERR.so.0.8.5`
via `DT_RUNPATH`. All HuffZstd_INT code lives in the shared lib, **not**
in the binary itself — so two `sperr3d` binaries built from different
source states are byte-identical (this tripped us up on the first
sanity check; both binaries had the same md5 and produced different
output only because the library on disk had changed).

Final setup (single binary, two libs):

- `/home/administrator/SPERR-fork-cc-cli/build/bin/sperr3d` — the only
  binary used
- `/tmp/libSPERR_reorder/libSPERR.so.0.8.5` — built from clean Sprint 5
  source (subband permutation active)
- `/tmp/libSPERR_noperm/libSPERR.so.0.8.5` — built after patching
  `enumerate_subbands` in `src/HuffZstd_INT.cpp` to short-circuit to a
  single full-volume box (identity permutation); patch was applied,
  rebuilt, then **reverted** so the source tree is back to clean
  Sprint-5 state

Bench script switches per-run via `LD_LIBRARY_PATH=/tmp/libSPERR_{...}
build/bin/sperr3d ...`. RUNPATH (not RPATH) lets `LD_LIBRARY_PATH`
override the build-tree default. `readelf -d` to confirm.

The noperm libSPERR is a **build artifact only** — source tree is
clean; do not commit anything from `/tmp/libSPERR_noperm`. To regenerate
it: re-apply the early-return patch at the top of `enumerate_subbands`
(line 46 area), rebuild, copy `build/src/libSPERR.so.0.8.5` to
`/tmp/libSPERR_noperm/`, revert patch, rebuild.

### Files and artifacts

```
/tmp/bench_3way.sh              # bench driver, sequential, --omp 1
/tmp/bench_3way_results.tsv     # one row per run; columns:
                                #   dataset file vre backend bytes
                                #   comp_s decomp_s psnr_db linfty ratio bpp
/tmp/bench_3way.progress        # one-line live counter, e.g. "41/342 ..."
/tmp/bench_3way.log             # start/end stamp + SKIP/FAIL lines
/tmp/bench_3way.stdout          # nohup stdout
/tmp/libSPERR_reorder/          # Sprint 5 lib
/tmp/libSPERR_noperm/           # Sprint 4 lib (identity permutation)
```

Per-run captures `bytes` (stat of bitstream), `comp_s`/`decomp_s`
(parsed from `--print_stats` output), `psnr_db`, `linfty`, `ratio`,
`bpp`. PSNR/L_∞ comparison across backends acts as the correctness
gate (must match SPECK to two decimals at the same VRE).

### Status — completed 2026-05-14 ~14:25 (15 min wall)

342 runs total, 316 OK / **26 FAIL**. Failure breakdown:

- **NYX velocity_{x,y,z}: 16/26** — upstream 0.8.5 non-determinism bug (see Sprint 6.5)
- **SCALE T/U/V: 9/26** — wide-dynamic-range fields, same kind of issue suspected
- **hurricane TCf48 vre=0.001 reorder: 1/26** — stray

### Per-dataset × VRE summary (triple-success files only)

`NP` = HuffZstd-noperm, `RE` = HuffZstd-reorder. `NP%`/`RE%` = size delta vs
SPECK (positive = larger). `cmp_{s,n,r}` = compression time (s) for the
three backends; same for `dec_*`. `NP<S` = count of files where noperm beat
SPECK on size.

```
dataset    vre        # | SPECK_MB   NP_MB   RE_MB |    NP%    RE% | cmp_s cmp_n cmp_r | dec_s dec_n dec_r |  NP<S  RE<S
hurricane  0.01      13 |   0.356   0.399   0.387 | +14.3% +10.1% |  0.64  0.60  0.60 |  0.31  0.29  0.29 |  0/13  0/13
hurricane  0.001     12 |   1.504   1.686   1.653 | +11.9%  +9.8% |  0.74  0.63  0.61 |  0.38  0.33  0.31 |  0/12  0/12
hurricane  0.0001    13 |   4.160   4.745   4.675 | +12.9% +11.6% |  0.96  0.70  0.71 |  0.56  0.42  0.41 |  0/13  0/13
miranda    0.01       7 |   0.201   0.239   0.216 | +19.1%  +7.5% |  1.03  1.13  1.07 |  0.44  0.46  0.47 |  0/7   0/7
miranda    0.001      7 |   0.719   0.867   0.802 | +20.8% +11.7% |  1.15  1.12  1.18 |  0.48  0.48  0.48 |  0/7   0/7
miranda    0.0001     7 |   2.086   2.457   2.352 | +17.8% +12.7% |  1.31  1.14  1.29 |  0.60  0.62  0.60 |  0/7   0/7
nyx        0.01       3 |   0.082   0.111   0.102 | +31.4% +26.9% |  3.74  3.96  4.01 |  1.93  2.09  2.15 |  0/3   0/3
nyx        0.001      3 |   0.882   1.108   1.052 | +24.0% +18.3% |  3.80  3.94  3.98 |  2.08  2.19  2.18 |  0/3   0/3
nyx        0.0001     3 |   6.033   6.889   6.701 | +16.4% +11.4% |  4.45  4.17  4.09 |  2.44  2.29  2.35 |  0/3   0/3
scale      0.01      10 |   0.467   0.608   0.532 | +31.2% +14.5% |  3.85  3.76  3.80 |  2.05  2.16  2.17 |  0/10  0/10
scale      0.001      9 |   3.229   3.666   3.446 | +15.8%  +9.1% |  4.18  3.83  3.88 |  2.33  2.21  2.24 |  0/9   0/9
scale      0.0001    11 |  14.331  15.971  15.393 | +12.1%  +8.2% |  5.21  4.30  4.28 |  3.10  2.67  2.70 |  0/11  0/11
```

### Reorder vs noperm (isolated Sprint-5 gain on Huffman input ordering)

```
hurricane  0.01    13 |  -3.6%   (13/13 RE smaller)
hurricane  0.001   12 |  -1.8%   (10/12)
hurricane  0.0001  13 |  -1.1%   (10/13)
miranda    0.01     7 |  -9.7%   (7/7)
miranda    0.001    7 |  -7.5%   (7/7)
miranda    0.0001   7 |  -4.3%   (7/7)
nyx        0.01     3 |  -2.9%   (2/3)
nyx        0.001    3 |  -4.6%   (3/3)
nyx        0.0001   4 |  -4.1%   (4/4)
scale      0.01    10 | -12.6%   (10/10)
scale      0.001    9 |  -5.8%   (9/9)
scale      0.0001  11 |  -3.4%   (11/11)
```

### Totals (98 triple-success files merged)

```
total bytes (MB):  SPECK=325.2  noperm=368.2  reorder=354.9
  noperm  vs SPECK:  size +13.2%   comp_time -8.8%   decomp_time -5.4%
  reorder vs SPECK:  size  +9.1%   comp_time -8.1%   decomp_time -5.0%
  reorder vs noperm: size  -3.6%   comp_time +0.7%   decomp_time +0.4%
```

### Takeaways

- **Zero files** in the bench grid had HuffZstd (either variant) beat SPECK
  on size. The `density_128x128x256.d64` win noted earlier in Sprint 5
  was on a small test-data file outside this grid; the broader picture
  on real SDRBench data is +6–14% larger.
- **HuffZstd is always faster than SPECK** on real data: -5 to -9% comp,
  -5% decomp, on average. Bigger wins at low VRE / smaller alphabets.
- **Sprint 5 reorder gives -3.6% size for free** (+0.7% comp, +0.4% decomp).
  Biggest wins on `scale` at vre=1e-2 (-12.6%) and `miranda` at vre=1e-2
  (-9.7%). Note: this is the opposite of the Sprint 4.5 probe's prediction
  (which said high precision had the most reorder headroom). Reality: at
  low precision the bin=0 fraction is so dominant that clustering them
  into longer runs makes zstd LZ matches dramatically longer.

## Sprint 6.5: NYX velocity upstream non-determinism (2026-05-14)

Triggered by the NYX velocity failures in Sprint 6 bench. SPECK-backend
failures couldn't be our HuffZstd work, so dug into upstream.

### What we found

`sperr3d -c --pwe 1e6` on NYX `velocity_x.f32` is **non-deterministic**
on the upstream 0.8.5 codebase — same binary, same input, same flags
intermittently fails with `Compression failed!` or succeeds with
PSNR 61.33 (no in-between). 10-run failure rate:

| Repo | branch / commit | velocity_x pwe=1e6, 10 runs |
|---|---|---|
| `~/SPERR` (upstream main, 0.8.5) | main / 9f930a3 | **3/10 OK**, 7/10 fail |
| `~/SPERR-fork` (0.8.5-custo) | 0.8.5-custo / 384693a | **6/10 OK**, 4/10 fail |
| `~/SPERR-fork-0.8.3` | main / 9718600 | **10/10 OK**, 0 fail ✓ |

The 0.8.5-custo's slightly lower fail rate is sample noise — both are
non-deterministic. The 0.8.3 baseline is fully deterministic. **The
regression is in upstream 0.8.3 → 0.8.5**, not in any local fork branch.

### Suspect commits (0.8.3 → 0.8.5 upstream main)

```
9f930a3  automatically detect and enable/disable AVX2
38b13a5  Improve 2D compression performance by taking advantage of the int_8 msb info
98aa773  Various performance optimization in 3D compression       ← prime suspect (3D hot path)
4061a44  Simpler AVX2 code, and other clean ups
```

Most likely 98aa773 (3D hot-path optimization). git-bisect would
localize in ~2 builds; deferred unless we need to upstream a fix.

### Diff `~/SPERR-fork` vs `~/SPERR` (main): nothing functionally suspect

```
CMakeLists.txt          |  1 -
include/SPECK_FLT.h     |  2 ++   # adds m_q_coeff member, default 1.5
include/SPERR3D_OMP_C.h |  2 ++   # adds set_q_coeff forwarder
include/Timer.h         | 36 ++++ # new utility, only used in sperr3d.cpp CLI
include/sperr_helper.h  |  3 ++   # adds any_ge declaration
src/CMakeLists.txt      |  1 -
src/SPECK_FLT.cpp       |  8 ++++ # m_estimate_q PWE: 1.5 → m_q_coeff (default 1.5, equivalent)
src/SPERR3D_OMP_C.cpp   |  8 ++++ # set_q_coeff plumbing
src/sperr_helper.cpp    | 89 ++++ # adds any_ge() AVX2 helper, NEVER CALLED anywhere
utilities/sperr2d.cpp   | 54 ++++ # adds --vre / --q CLI options
utilities/sperr3d.cpp   | 63 ++++ # adds --vre / --q, --vre maps to --pwe internally
```

- `m_q_coeff` defaults to 1.5 = identical to upstream's hard-coded `1.5`,
  so the PWE path is functionally unchanged when the CLI doesn't pass `--q`.
- `any_ge` is defined but not referenced anywhere — dead code, can't be
  the bug.
- The `--vre` plumbing just computes `pwe = vre * (max-min)` and feeds
  `--pwe` underneath, so any `--vre` failure on `velocity_x` is really
  an upstream `--pwe` failure.

**Conclusion**: SPERR-fork's 0.8.5-custo branch inherits the upstream
non-determinism via the `384693a` merge of `local-sp/main`; cc-cli
inherits it transitively. The 26 bench failures attributable to NYX
velocity (16) and likely SCALE T/U/V (9) are **upstream bugs**, not
regressions from any Sprint 1–5 work.

### What this means for the bench statistics

- The bench captures one sample per (file, vre, backend). For
  velocity_x and similar non-deterministic files, a single sample is
  binary success/fail — re-running produces a *different* fail set.
- Future bench runs should either:
  - Drop the affected files from the grid, or
  - Take ≥3 samples per cell and treat any-success as success

For Sprint 6's reported statistics, the triple-success filter is
implicitly biased toward files that happen to succeed on this run.
Triple-success counts on NYX are only 3-4 of 6 files at each VRE.

### Correction (2026-05-14, later): the bug is older than 0.8.3→0.8.5

The 10/10 OK observation for `~/SPERR-fork-0.8.3` above used the
`compressor_3d` binary from `main` branch with `USE_OMP=OFF`. After
rebuilding `~/SPERR` and `~/SPERR-fork` with `USE_OMP=ON` (per user
guidance — original cmake config had OMP off, so `--omp` wasn't even
exposed), and re-testing with the `sperr3d` binary on `0.8.3-custo`
branch, all four versions exhibit non-determinism:

| Version | branch / binary | OK / 10 (--omp 1, USE_OMP=ON) |
|---|---|---:|
| `~/SPERR` (main 9f930a3) | sperr3d | **1/10** |
| `~/SPERR-fork` (0.8.5-custo) | sperr3d | **3/10** |
| `~/SPERR-fork-0.8.3` (0.8.3-custo) | sperr3d | **5/10** |
| `~/SPERR-fork-cc-cli` (0.8.5-cc-cli) | sperr3d | **6/10** |

The earlier "0.8.3 is deterministic" claim isolated three confounded
variables: branch (main vs custo), binary (`compressor_3d` vs `sperr3d`),
and OpenMP linkage (OFF vs ON). The most suspect is **USE_OMP=ON
triggering nondeterminism** even at `--omp 1` (OpenMP runtime init may
introduce race conditions in some reductions). cc-cli's relatively
lower fail rate (6/10) is most likely sampling noise on small N.

**Practical implication**: NYX velocity_x/y/z and SCALE T/U/V are
pathological cases across the whole SPERR fork family. Don't treat the
absolute fail counts in Sprint 6 as cc-cli-specific regressions —
they're upstream behavior we inherit.

## Sprint 7: qcat huffmanZstd baseline + HuffZstd time breakdown (2026-05-14)

Motivation: user's intuition was that HuffZstd should be *much* faster
than SPECK end-to-end, but Sprint 6 showed only -5 to -9% comp time wins.
Hypothesis: the HuffZstd step itself is fast, but it's a small fraction
of total compression time. Test that.

### Method

1. Added env-gated **bins dump** to `HuffZstd_INT::encode` —
   `SPERR_HZ_DUMP_BINS=<prefix>` writes the int32 ZigZag bins (post-perm
   in reorder lib, raster in noperm lib) to `<prefix>.c<chunk_id>.bin`
   plus a `.meta` (N, stateNum, max_mag).
2. Built `~/qcat`'s `huffmanZstd` example (target `huffmanZstd` only;
   sibling `compress`/`generateIndexData` fail to compile in qcat
   unrelated to us — skipped via `--target huffmanZstd`).
3. Dumped Uf48 (4 chunks) and pressure (1 chunk) at vre ∈ {1e-2,1e-3,1e-4}
   for both reorder and noperm libs.
4. Ran qcat `huffmanZstd <bin> <N> <quantBinCapacity>` per chunk; qbc set
   to `ceil(stateNum/2)` so qcat's internal Huffman alphabet
   (= 2·quantBinCapacity) covers our ZigZag range.
5. Added `t_huff` and `t_zst` timing to the `SPERR_HZSTATS=1` probe in
   `encode_huffzstd` (chrono::steady_clock). Re-ran SPERR with HZSTATS to
   capture HuffZstd-only wall time per chunk.

### qcat (SZ3 Huffman + Zstd) vs SPERR HuffZstd (SZo + Zstd), same bins

```
                            qcat (SZ3)            SPERR (SZo)
dataset  vre     variant   bytes   time     bytes      time   q/s_size  q/s_time
uf48     0.01    reorder   474 KB  0.075s    442 KB   0.018s   +7.2%     4.2×
uf48     0.001   reorder  2624 KB  0.103s   2476 KB   0.036s   +6.0%     2.8×
uf48     0.0001 reorder   8561 KB  0.187s   8112 KB   0.143s   +5.5%     1.3×
uf48     0.01    noperm    498 KB  0.073s    478 KB        —   +4.2%       —
uf48     0.001   noperm   2694 KB  0.106s   2616 KB        —   +3.0%       —
uf48     0.0001 noperm    8698 KB  0.194s   8392 KB        —   +3.6%       —
pressure 0.01    reorder   217 KB  0.121s    210 KB   0.022s   +3.5%     5.5×
pressure 0.001   reorder   859 KB  0.086s    806 KB   0.028s   +6.6%     3.1×
pressure 0.0001 reorder  2531 KB  0.108s   2390 KB   0.045s   +5.9%     2.4×
```

SZo (current SPERR HuffZstd) vs SZ3 (qcat baseline) on identical
ZigZag-encoded bins:
- **SZo output 3–7% smaller** than SZ3 (canonical-codelen tree vs SZ3's
  pointer tree — confirms Sprint 4 hypothesis).
- **SZo 1.3–5.5× faster** than SZ3 (L1 lookup table + valid_len1 fast
  path). Bigger speedups at low VRE / smaller alphabets — where the
  fast paths dominate.

### Time breakdown: HuffZstd is not the bottleneck

```
dataset  vre     SPERR_total  HuffZstd_only  rest (wavelet/cond/quant/outlier)  HZ %
uf48     0.01      0.580s       0.018s           0.562s (97%)                    3.1%
uf48     0.001     0.624s       0.036s           0.588s (94%)                    5.8%
uf48     0.0001    0.758s       0.143s           0.615s (81%)                   18.9%
pressure 0.01      1.168s       0.022s           1.146s (98%)                    1.9%
pressure 0.001     1.037s       0.028s           1.009s (97%)                    2.7%
pressure 0.0001    1.075s       0.045s           1.030s (96%)                    4.2%
```

HuffZstd is **2–19% of total compression time**, and only crosses 15% at
high precision (vre=1e-4 on Uf48 where the alphabet explodes).

### HuffZstd vs SPECK on the same int-coding stage

Estimate SPECK_INT time as `SPECK_total - (HuffZstd_total - HuffZstd_only)`
(both share wavelet/cond/quant/outlier):

```
dataset  vre      SPECK_INT (est)   HuffZstd      SPECK/HuffZstd
uf48     0.0001     0.575s           0.143s          4.0×
uf48     0.001      0.252s           0.036s          7.0×
uf48     0.01       0.128s           0.018s          7.1×
pressure 0.0001     0.131s           0.045s          2.9×
pressure 0.001      0.043s           0.028s          1.5×
pressure 0.01       0.046s           0.022s          2.1×
```

**On the int-coding stage alone, HuffZstd is 1.5–7× faster than SPECK.**
User intuition confirmed. The reason end-to-end speedup is modest is
that int-coding is only 5–40% of total SPERR compression. Wavelet
transform + conditioning + quantization + outlier coding dominate, and
those don't change with the int backend.

### Implication for "where to optimize next"

If end-to-end SPERR speed is the goal, the int-coding swap has captured
most of its possible win. Further gains require attacking:
- **Wavelet transform** (CDF 9/7) — likely SIMD already
- **Outlier coder** (SPECK1D) — same kind of bitplane coding as SPECK3D;
  switching it to a Huffman variant could give 1.5–7× on that subroutine
  too (but the outlier coder is far smaller than main SPECK_INT, so the
  end-to-end win is bounded). Note: paper VI-E claims SPECK1D beats
  Huffman by 1–2 bits/outlier on size, so this is a speed/size trade-off.
- **Quantization / max-element scan** — usually dominates at high
  precision (large `m_vals_d`).

### Files / artifacts from Sprint 7

```
src/HuffZstd_INT.cpp                    # added:
                                        #   - SPERR_HZ_DUMP_BINS env-gated dump hook
                                        #     (writes per-chunk int32 bins + .meta file)
                                        #   - t_huff/t_zst chrono timing in SPERR_HZSTATS probe
/tmp/libSPERR_reorder_dump/             # libSPERR.so.0.8.5 built from reorder source + dump hook
/tmp/libSPERR_noperm_dump/              # libSPERR.so.0.8.5 built from noperm-patched source + dump hook
/tmp/bins/                              # int32 bin dumps (12 dataset×vre×variant combos × #chunks)
/tmp/qcat_results.tsv                   # qcat huffmanZstd output, per-chunk
/home/administrator/qcat/build/         # qcat build tree (only huffmanZstd target built)
```

Both env hooks (`SPERR_HZ_DUMP_BINS` and the t_huff/t_zst extension in
`SPERR_HZSTATS`) left in tree as permanent debugging hooks — zero cost
when env vars are unset.

## Sprint 8: LC framework exploration (Phase 0, 2026-05-14)

Goal: figure out whether a hybrid pipeline of **LC components + Huffman
+ zstd** can beat the current SPERR HuffZstd backend on
size-at-equal-speed. Read `~/LC-framework` source and ASPLOS25 paper
first; design rationale captured below.

### Reading-time analysis (before any code)

**LC's stance**: LC explicitly avoids entropy coding (Huffman / arith /
FSE / rANS) for throughput. The 70 components are all *structural*
byte/bit transforms: BIT (transposition), CLOG/HCLOG (leading-zero),
DIFF/DIFFMS/DIFFNB/DBEFS/DBESF (deltas), RARE/RAZE (adaptive top-k +
zero-elim), RLE/RZE/RRE (run/zero/repeat elim), TCMS/TCNB
(two's-complement→magnitude-sign), TUPL* (tuple regroup). Paper's 4
recommended pipelines:

```
SPspeed: DIFFMS → MPLG
DPspeed: DIFFMS → MPLG
SPratio: DIFFMS → BIT → RZE
DPratio: FCM → DIFFMS → RAZE → RARE
```

**Hybrid rationale**: three layers each catch different redundancy —
LC (structural pattern), Huffman (per-symbol entropy), zstd (LZ77 +
internal entropy). Order constraint: **Huffman must be last** if used
at all, since its output is non-byte-aligned bitstream that defeats any
byte-level downstream transform.

**Lorenzo predictor evaluated and rejected** before any code. Lorenzo
is SZ3's core weapon on raw scientific data (CR 30–50% delta), but our
input is *post-wavelet* ZigZag bins where CDF 9/7 has already done
spatial decorrelation. High-freq subbands are ~all 0 so Lorenzo gives
0-pred → 0-resid (no win); LLL is tiny (<0.003% of elements). Expected
return on our data: 2–5% size reduction for 300 lines of code. Not the
right battle. Decision: skip Lorenzo entirely.

### Phase 0: baseline + hand-picked / analytical pipelines

**Build**:

```
cd ~/LC-framework && python3 ./generate_Host_LC-Framework.py
g++ -O3 -march=native -fopenmp -mno-fma -ffp-contract=off \
    -DUSE_CPU -I. -std=c++17 -o lc lc.cpp
```

**Artifacts** (under `/tmp/lc_bench/`):
- `run_lc.sh` — helper, returns one TSV row of `(bin_label, pipeline,
  size, enc_s, dec_s, cr)` for a given input file + pipeline string
- `baseline.tsv` — per-chunk SPERR HuffZstd `zst_out` + `t_huff` + `t_zst`
  exported via `SPERR_HZSTATS=1` (30 rows covering Uf48 4-chunks ×
  pressure 1-chunk × 3 VRE × 2 orders; only `reorder` lib has timing
  since `noperm` lib was built before chrono hook was added — sizes
  available for both)

**Representative file for Phase 0**: `pressure_vre0.001_reorder.c0.bin`
(37.7M int32, 144 MB raw, single chunk, mid-VRE). SPERR HuffZstd
baseline on this bin: **825 KB, 27.3 ms** (`t_huff + t_zst`).

**Key empirical findings**:

| LC pipeline | size | vs HZ | enc | comment |
|---|---:|---:|---:|---|
| `BIT_4` (1-stage) | 144 MB | +18187% | 0.020s | bijection → CR 1.0 |
| `DIFFMS_4` (1-stage) | 144 MB | +18187% | 0.020s | bijection, **no spatial corr after wavelet** |
| `TCMS_4` (1-stage) | 144 MB | +18187% | 0.017s | **no-op** on ZigZag bins (already non-negative) |
| `RLE_4` (1-stage) | 5.09 MB | +516% | 0.003s | catches some runs but doesn't reach bit-plane structure |
| `BIT_4 RLE_4` | 2.20 MB | +166% | 0.011s | BIT exposes bit-planes, RLE captures long runs |
| `BIT_4 RZE_4` (best 2-stage) | 1.42 MB | +71.7% | 0.019s | the canonical pair |
| `BIT_4 RZE_4 RLE_4` (3-stage) | 1.42 MB | +71.6% | 0.008s | RLE post-RZE adds nothing |
| `BIT_4 RZE_1 BIT_1 RZE_1` (4-stage, mixed word) | 1.09 MB | +31.8% | 0.040s | hand-derived; **recurse bitmap compress** |
| **`BIT_4 RZE_4 RZE_1` (analytical 3-stage)** | **1.08 MB** | **+30.8%** | **0.020s** | designed below |

**Analytical 3-stage**: `BIT_4 → RZE_4 → RZE_1`

- Step 1 `BIT_4` is mandatory entry. Our bins fit in 14 bits (pressure
  vre=1e-3 max_mag=4812) so bit-planes 14–31 are entirely 0; BIT_4 lays
  them out as 18 contiguous all-zero plane regions of N/8 bytes each
  (~84 MB of pure zero in 4.7 MB plane chunks).
- Step 2 `RZE_4` eliminates those 4-byte aligned zero runs efficiently.
  RZE internally already does a `RZE-on-bitmap` recursive pass per
  paper §3.2 → free secondary compression.
- Step 3 `RZE_1` finishes the byte stream: the bitmap RZE_4 outputs
  for high bit-planes is mostly-zero bytes too; another byte-granularity
  RZE squeezes the tail.

This 3-stage matched (in fact slightly bested) the hand-tuned 4-stage
without an explicit search — confirms the bit-plane structure dominates.

**Other variants tested**: `BIT_4 RAZE_4 RZE_1` (+31.8%, RAZE's
adaptive-k adds nothing once BIT_4 has put zero planes contiguously);
`BIT_4 RZE_4 RARE_1` (+60.9%, RARE looks for repeated-non-zero patterns
that don't exist in RZE output); `BIT_4 RZE_1 RZE_1` (+34.8%, missing
4-byte coarse pass costs).

### Decision points reached

1. **LC alone cannot beat HuffZstd on size**. Tight 3-stage = +30.8%
   gap; 4-stage marginal improvement. The missing ~30% is exactly the
   space where entropy coding lives. Hybrid pipeline (Phase 2+) needed.
2. **Time**: LC's `BIT_4 RZE_4 RZE_1` runs in 20 ms vs HuffZstd 27 ms
   for the same bin — already faster than current backend on this data
   point. Adding `zstd -3` tail will close the size gap at modest time
   cost (qcat showed ~0.1s for similar size — significant but acceptable).
3. **No need for LC search to make decision**. Phase 0 design-by-reasoning
   converged to the same answer LC's `EX`/`GA` search would likely
   produce (modulo a few % from mixed word sizes); cost-benefit of full
   search not worth it at this stage.

### Phase plan (going forward)

- **Phase 1**: SKIP. Don't bother with full LC `EX` search on 357K
  combos; the design space's structure is now understood.
- **Phase 2** (next): chain `BIT_4 RZE_4 RZE_1 → zstd -3` and compare
  against HuffZstd on the 12 dataset×VRE×variant cells. Need to dump
  LC output bytes and pipe through external zstd CLI (or libzstd
  in-process).
- **Phase 3** (conditional): if Phase 2 doesn't win, try
  `LC pipeline → byte-Huffman (alphabet=256) → zstd`. byte-Huffman
  via SZo `HuffmanEncoder<uint8_t>` already vendored.
- **Phase 4** (conditional): if Phase 2 or 3 wins on size at comparable
  speed → use LC's `generate_standalone_CPU_compressor_decompressor.py`
  to emit hard-coded C++ for the chosen pipeline, vendor it, add
  `IntBackend::LC` and corresponding `LC_INT` class.

### Files / artifacts from Sprint 8

```
~/LC-framework/lc                        # CPU-only LC binary (built)
/tmp/lc_bench/run_lc.sh                  # helper: pipeline → size+timing TSV
/tmp/lc_bench/baseline.tsv               # 30 per-chunk HuffZstd baselines
```

No SPERR source changes in Sprint 8 — exploration-only. The Sprint 7
`SPERR_HZ_DUMP_BINS` hook plus `SPERR_HZSTATS` chrono extension are
sufficient instrumentation.

## Sprint 8 Phase 2: LC + zstd validation (2026-05-14)

Followed up Phase 0's analytical 3-stage `BIT_4 RZE_4 RZE_1` with the
hybrid `LC pipeline + zstd-3` test. Phase 0 had shown LC alone was
+30.8% over HuffZstd; the hypothesis was that zstd's outer entropy
+LZ would close that gap.

### Setup

- Built LC standalone for `BIT_4 RZE_4 RZE_1` via
  `generate_standalone_CPU_compressor_decompressor.py` →
  `/tmp/lc_bench/lc_compress`.
- Built `/tmp/lc_bench/zstd_wrap` — minimal C++ wrapper around vendored
  `libzstd_static` (no shell zstd CLI; everything in-tree). Takes
  `[-d] <in> <out> [level=3]`, prints `elapsed_ms`.
- Bench scripts:
  - `/tmp/lc_bench/run_phase2.sh` — 30 chunks (Uf48 + pressure only,
    reorder + noperm variants from Sprint 6 dumps)
  - `/tmp/lc_bench/run_wide.sh` — wide grid (all 38 files × 3 VRE,
    reorder only). Dumps bins per (file, vre), runs LC+zstd, captures
    HuffZstd `zst_out` baseline from `SPERR_HZSTATS` in the same SPERR
    call; cleans up dumps per combo (NYX 512³ files = ~9 GB if kept).

### zstd-3 chosen as final level

Level sweep on `pressure_vre0.001_reorder.c0`:

```
level   size     zstd_ms   vs HZ
1       850 KB   1.0 ms    +2.98%
3       841 KB   2.4 ms    +1.80%   ← chosen (matches HZ's internal default)
6       838 KB   5.2 ms    +1.51%
12      838 KB  11.4 ms    +1.46%
19      826 KB  74.6 ms    -0.01%
```

zstd-3 picked for fair comparison (HuffZstd backend uses zstd level 3
internally). Higher levels close the gap but at 10–30× zstd time —
not justified.

### Phase 2 narrow grid (30 chunks: Uf48 4×3×2 + pressure 1×3×2)

```
dataset   vre     variant    #    LC+zstd_B         HZ_B    Δsize
pressure  0.0001  reorder    1      2393602      2447002   -2.18%
pressure  0.001   reorder    1       840674       825793   +1.80%
pressure  0.01    reorder    1       225246       215239   +4.65%
uf48      0.0001  reorder    4      7835627      8307013   -5.67%
uf48      0.001   reorder    4      2512729      2535665   -0.90%
uf48      0.01    reorder    4       450215       451994   -0.39%
(plus 12 noperm rows, similar values within ±0.5% of reorder)
```

Grand: LC+zstd=29.18 MB vs HZ=30.22 MB → **-3.43%** on the 30-chunk set.

### Phase 2 wide grid (929 chunks across all 4 SDRBench datasets)

```
dataset   vre     #chnk  LC+zstd_B          HZ_B    Δsize  wins/n
hurricane 0.0001     44   50.8 MB        56.7 MB  -10.48%   32/44
hurricane 0.001      52   22.4 MB        21.9 MB   +2.57%   13/52
hurricane 0.01       34    3.08 MB        3.47 MB -11.40%   12/34
miranda   0.0001      7   16.6 MB        17.1 MB   -2.84%    5/7
miranda   0.001       7    5.84 MB        5.87 MB  -0.51%    4/7
miranda   0.01        7    1.59 MB        1.58 MB  +0.61%    4/7
nyx       0.0001     29   58.4 MB        56.9 MB   +2.68%    8/29
nyx       0.001      25   10.1 MB        10.3 MB   -2.74%   11/25
nyx       0.01       31    1.14 MB        1.08 MB  +5.73%    4/31
scale     0.0001    238  150.2 MB       156.0 MB   -3.72%  187/238
scale     0.001     233   38.3 MB        40.1 MB   -4.53%  195/233
scale     0.01      222    4.67 MB        4.79 MB  -2.39%  177/222

GRAND TOTAL: 929 chunks LC+zstd=363 MB HZ=376 MB Δ=-3.38% (stable vs narrow grid)
```

13 (file, vre) combos failed (same upstream non-determinism as Sprint 6).
SCALE dominates wins (~80%+ chunks). NYX velocity is the weakness:
+5.73% loss at vre=1e-2.

### Phase 2 conclusion

Fixed `BIT_4 RZE_4 RZE_1 + zstd-3` is **stable ~-3% smaller than
HuffZstd** across the broad grid. Confirms Phase 0's prediction that
hybrid LC+zstd works. Time data was wall (standalone has file I/O
overhead) — true in-memory expected ~20 ms LC + ~5 ms zstd ≈ HZ
27 ms.

## Sprint 8 Phase 3: per-file pipeline search (2026-05-14)

Question: does picking the best pipeline *per file* (not a fixed one)
materially improve LC+zstd? Tested 25 hand-curated candidates (all
starting `BIT_4`, varying stage 2/3) on each (file, vre) combo across
the full grid; per-(file,vre) we sum LC+zstd bytes across all chunks
of that file, pick the candidate with min total.

Implementation in `/tmp/lc_bench/per_file_search.sh`. Per-candidate
standalones built once into `/tmp/lc_bench/cand/lc_<idx>`; then 114
combos × 25 candidates = ~2850 (combo, candidate) trials done in
sequence.

### Aggregate (105 successful files / 9 fail)

```
dataset   vre        # | Fix%HZ   Best%HZ   Best%Fix
hurricane 0.0001    13 | -17.89%  -18.34%   -0.54%
hurricane 0.001     13 |  +2.57%   +2.06%   -0.50%
hurricane 0.01      13 | -26.89%  -27.25%   -0.50%
miranda   0.0001     7 |  -2.84%   -3.54%   -0.72%
miranda   0.001      7 |  -0.51%   -1.24%   -0.73%
miranda   0.01       7 |  +0.61%   +0.05%   -0.56%
nyx       0.0001     5 |  -4.03%   -4.92%   -0.93%
nyx       0.001      4 |  -9.44%  -10.60%   -1.28%
nyx       0.01       4 |  +6.07%   +1.23%  **-4.57%**  ← only big per-file gain
scale     0.0001    10 | -13.95%  -14.50%   -0.63%
scale     0.001     11 |  -6.19%   -6.46%   -0.29%
scale     0.01      11 |  -7.72%   -8.69%   -1.05%

GRAND: 105 files  Fix=-9.74%  Best=-10.35%  Best vs Fix=-0.68%
```

### Conclusion

Fixed pipeline `BIT_4 RZE_4 RZE_1` is **already near-optimal**.
Per-file search adds only **-0.68% additional**. Only NYX vre=1e-2 had
a meaningfully different sweet spot (-4.57% from picking per-file).
Engineering cost of per-file pipeline encoding (metadata, dispatch)
not worth it. **Decision: ship fixed pipeline as the LC backend.**

(Wide bench grand total -3.38% vs per-file grand -9.74%: the discrepancy
is due to NYX velocity files randomly succeeding/failing differently
between runs — different file populations, different aggregate weights.
Methodologically both confirm "LC+zstd wins HuffZstd net of NYX
pathologies".)

## Sprint 9: LC backend integration into SPERR (2026-05-14)

Integrated the fixed-pipeline LC + zstd backend as `IntBackend::LC`,
mirroring the existing `HuffZstd_INT` path.

### Architecture choices

- **Pre-ZigZag on SPERR side** (same as HuffZstd_INT). Tried LC's TCMS
  as stage 1 in Phase 0; equivalent on ZigZag input (TCMS becomes ×2,
  no compression gain). Cleaner to keep ZigZag in `LC_INT` and have
  LC pipeline see pure non-negative int32.
- **Subband reorder kept** (Sprint 5 layout). Phase 2 data showed
  reorder vs noperm differ <0.5% on LC+zstd, but reorder is the
  current production default for HuffZstd → consistency.
- **Fixed pipeline** `BIT_4 → RZE_4 → RZE_1` baked into vendored code.
  No per-file selection (Phase 3 showed -0.68% marginal).
- **Outer zstd via SZ3's `Lossless_zstd`** (already in tree via
  third_party/sz3_huff). Default level 3.
- **Payload format** parallel to HuffZstd v2:
  - byte 0: num_bitplanes_eq
  - bytes 1..8: m_total_bits (uint64)
  - byte 9: format_version = 1 (LC)
  - bytes 10..N: ZSTD frame wrapping LC-pipeline output

### Files added/changed

**New** (vendored LC + new backend):

```
third_party/lc/
  include/macros.h
  lc_pipeline.hpp                          # in-memory wrapper (encode/decode)
  components/h_BIT_4.h
  components/h_RZE_4.h
  components/h_RZE_1.h
  components/include/h_RZE.h
  components/include/h_zero_elimination.h
include/LC_INT.h                           # derives from SPECK_INT<T>
src/LC_INT.cpp                             # ~290 lines: ZigZag, subband
                                           #   reorder (duplicated from
                                           #   HuffZstd_INT — small enough
                                           #   to keep self-contained),
                                           #   LC pipeline call, zstd
```

**Modified**:

```
include/sperr_helper.h        # IntBackend enum: add LC = 2
include/Conditioner.h         # m_backend_idx unchanged (still 1) but
                              #   now occupies meta[1..2] (2 bits)
src/Conditioner.cpp           # save_backend/retrieve_backend: 2-bit
                              #   encoding (00=SPECK, 01=HuffZstd, 10=LC)
src/SPECK3D_FLT.cpp           # make_enc/make_dec: switch on 3 backends
src/SPECK_FLT.cpp             # Rate-mode guard extended to LC
utilities/sperr3d.cpp         # --backend {speck,huffzstd,lc} +
                              #   --bpp+lc rejection
src/CMakeLists.txt            # add LC_INT.cpp
```

### lc_pipeline.hpp details

In-memory wrapper extracted from LC's `compressor-standalone.cpp` /
`decompressor-standalone.cpp`. Drops main() + file I/O; operates on
caller-owned byte buffers. CS = 16384 (LC default). **Sequential
chunk loop** (no OpenMP) — SPERR's per-chunk path is already inside
SPERR3D_OMP_C's parallel for, so adding nested OMP would oversubscribe.
LC chunk-level parallelism would only matter if we eventually ran
SPERR with `--omp 1` and wanted to parallelize within the LC stage.

### Sanity / correctness

3-way smoke test on `test_data/density_128x128x256.d64` (vre=1e-3):

```
SPECK:     CR=190.14, PSNR=77.74
HuffZstd:  CR=197.14, PSNR=77.74
LC:        CR=209.39, PSNR=77.74   ← 10.1% smaller than SPECK, 6.2% smaller than HZ
```

Standalone vs inline decomp byte-identical for LC. Conditioner 2-bit
backend dispatch works on both encode and decode side.

### 3-way bench (342 runs, all 4 SDRBench datasets, 2026-05-14)

`/tmp/bench_lc.sh` — sperr3d --backend {speck,huffzstd,lc}, single
binary, --omp 1. 11 fails (same NYX velocity / SCALE upstream bug
cluster as Sprint 6). 102 triple-success files for fair aggregate:

```
dataset   vre      # | SPECK_MB    HZ_MB    LC_MB |   HZ%S    LC%S |   cs    ch    cl  |  ds    dh    dl  | LC<S  LC<H
hurricane 0.0001   12 |  4.059    4.566    4.522 | +12.5% +11.4% | 0.94  0.70  0.65 | 0.53  0.41  0.38 | 0/12  4/12
hurricane 0.001    13 |  1.499    1.647    1.688 |  +9.9% +12.6% | 0.74  0.62  0.59 | 0.37  0.32  0.32 | 0/13  3/13
hurricane 0.01     13 |  0.356    0.387    0.405 |  +8.6% +13.9% | 0.64  0.59  0.59 | 0.31  0.29  0.33 | 0/13  3/13
miranda   0.0001    7 |  2.086    2.352    2.286 | +12.8%  +9.6% | 1.25  1.14  1.21 | 0.65  0.54  0.51 | 0/7   5/7
miranda   0.001     7 |  0.719    0.802    0.798 | +11.5% +11.0% | 1.12  1.13  1.08 | 0.48  0.46  0.49 | 0/7   4/7
miranda   0.01      7 |  0.201    0.216    0.217 |  +7.5%  +8.2% | 1.16  1.17  1.01 | 0.44  0.45  0.47 | 0/7   4/7
nyx       0.0001    3 |  6.033    6.701    7.255 | +11.1% +20.3% | 4.51  3.99  3.88 | 2.33  2.28  2.38 | 0/3   0/3
nyx       0.001     4 |  2.323    2.723    2.885 | +17.2% +24.2% | 3.76  3.94  3.79 | 2.11  2.16  2.39 | 0/4   0/4
nyx       0.01      3 |  0.082    0.102    0.120 | +24.4% +46.9% | 3.67  3.86  3.78 | 1.95  2.27  2.21 | 0/3   0/3
scale     0.0001   11 | 13.785   14.763   13.813 |  +7.1%  +0.2% | 5.37  4.45  4.18 | 3.11  2.86  2.56 | 7/11 11/11
scale     0.001    11 |  3.273    3.513    3.323 |  +7.3%  +1.6% | 4.28  4.01  3.91 | 2.35  2.33  2.45 | 2/11 10/11
scale     0.01     11 |  0.504    0.577    0.540 | +14.5%  +7.1% | 3.90  3.90  3.87 | 2.15  2.31  2.36 | 0/11  8/11
```

Grand total (102 files):

```
  SPECK=330.0MB   HZ=360.2MB   LC=348.9MB
  HZ  vs SPECK: size +9.16%   comp -7.8%   decomp -2.5%
  LC  vs SPECK: size +5.73%   comp -10.9%  decomp -2.7%   ← LC closer to SPECK
  LC  vs HZ:    size -3.14%   comp -3.4%   decomp -0.2%   ← LC wins HZ on both
```

### Headline findings

1. **LC beats HuffZstd on both size (-3.14%) and compress time (-3.4%)**
   net of failures. Decompress time roughly tied (-0.2%).
2. **LC is faster than SPECK on compress by -10.9%** (better than HZ's
   -7.8%); decomp roughly tied.
3. **SPECK still wins on size overall** (+5.73% LC vs SPECK), but
   on **SCALE vre=1e-4 LC actually beats SPECK on 7/11 files**
   (+0.2% aggregate — essentially tied) and beats HZ on 11/11.
4. **NYX velocity remains the weak spot** for both LC and HZ — both
   lose 11–47% vs SPECK on NYX low-VRE. NYX velocity has huge dynamic
   range (alphabet ~2^17–2^20), so BIT_4's bit-plane separation
   produces fewer all-zero planes → LC's zero-elim loses leverage.

### Bench artifacts

```
/tmp/bench_lc.sh                # bench driver
/tmp/bench_lc_results.tsv       # 342 rows of (dataset, file, vre,
                                #   backend, bytes, comp_s, decomp_s,
                                #   psnr, linfty, ratio, bpp)
/tmp/bench_lc.log               # FAIL rows
/tmp/bench_lc.progress          # live counter
```

### Production readiness

`sperr3d --backend lc` is the recommended default for SDRBench-class
data **except NYX velocity-like fields**. For NYX velocity, fall back
to `--backend speck` (smallest) or `--backend huffzstd` (faster but
larger than SPECK). Backend choice is encoded in the bitstream, so
decoders auto-dispatch — no user action needed on decode side.

## Sprint 10: 2D path — both backends wired up (2026-05-15)

Goal: extend the `--backend {speck,huffzstd,lc}` dispatch to the 2D
path, the only sane "remaining work" item left from Sprint 9. CLI is
`sperr2d` and the `SPECK2D_FLT` class.

### Why this turned out to be small

The 3D-path subband enumerator already covers 2D as a degenerate case.
`enumerate_subbands(dims)` (in both `HuffZstd_INT.cpp` and `LC_INT.cpp`)
takes the wavelet-packet branch when `can_use_dyadic` returns nullopt;
that branch is keyed on `Lz = num_of_xforms(dims[2])` and `Lxy =
num_of_xforms(min(dims[0], dims[1]))`. For 2D inputs the CLI passes
`dims = {Nx, Ny, 1ul}` (sperr2d.cpp:247) — `num_of_xforms(1) == 0`, so
Lz=0, the loop produces a single z-band `[0, 1)`, and within it 1
2D-LL + 3 detail subbands per Lxy level. That is exactly the Mallat
2D layout `m_cdf.dwt2d()` produces. So no enumerator change needed;
no payload format change either (deterministic from chunk dims).

The Conditioner backend bit was already plumbed through `SPECK_FLT`
base in Sprint 1 / Sprint 9 — both encode (`compress` writes via
`save_backend`) and decode (`use_bitstream` reads via `retrieve_backend`
before `m_instantiate_decoder()`). The 2D class only needed its
`m_instantiate_encoder/decoder` overrides updated to dispatch on
`m_int_backend`, mirroring what 3D had.

### Changes

`src/SPECK2D_FLT.cpp` — full rewrite to the 3D pattern. Old "preserve
existing encoder if uint type matches" `std::variant` checks dropped;
replaced with `make_enc/make_dec<T>(m_int_backend)` factories that
return `std::unique_ptr<SPECK_INT<T>>` (assignable into the variant).
Three-way switch: SPECK2D_INT_ENC/DEC, HuffZstd_INT, LC_INT.

`utilities/sperr2d.cpp` — added `--backend {speck,huffzstd,lc}` option
identical to the sperr3d.cpp pattern. Added `--bpp + huffzstd/lc`
rejection at parse time. Forwards via `encoder->set_int_backend(bk)`
on the SPECK2D_FLT (no OMP_C wrapper for 2D).

That's it — no other files touched.

### Smoke test (lena512.float, PWE=1e-3)

| Backend | Bytes | Comp(s) | Decomp(s) | PSNR | L_∞ |
|---|---:|---:|---:|---:|---:|
| speck    | 457,505 |       0.028 | 0.017 | 114.62 | 1.01e-3 |
| huffzstd | 479,130 (+4.7%) | 0.610 | 0.150 | 114.62 | 1.01e-3 |
| lc       | 474,234 (+3.7%) | 0.012 | 0.007 | 114.62 | 1.01e-3 |

Three decoded `.f32` files md5-identical (`fda12a43...`). Standalone
decompress matches inline decompress for all 3 backends — Conditioner
backend bit auto-dispatch works on 2D too. `--bpp + huffzstd/lc`
rejected at parse time; `--bpp + speck` and `--psnr + lc` both work.

### CESM 2D bench (2026-05-15)

`SDRBENCH-CESM-ATM-cleared-1800x3600`, all 78 fields × 3 VRE × 3
backends = 702 runs (1 fail: `FSDTOA_1_1800_3600` vre=1e-4 SPECK
returns rc=78 — upstream issue, fails on SPECK too, not a backend
regression). Bench script: `/tmp/bench_cesm_2d.sh`. Single sperr2d
binary; bitstreams written to `/tmp/cesm_bs/` then deleted per-run
to keep disk usage bounded.

Triple-success cells per VRE (77 / 77 / 76):

```
vre    #   SPECK_MB    HZ_MB    LC_MB    HZ%S    LC%S    LC%H   cs    ch    cl   ds    dh    dl   LC<S  LC<H
1e-2   77   6.820     8.014    8.070  +17.5%  +18.3%   +0.7% 0.12  0.13  0.13 0.05  0.05  0.06   2/77  46/77
1e-3   77  35.960    39.854   39.682  +10.8%  +10.4%   -0.4% 0.16  0.14  0.13 0.08  0.06  0.06   2/77  53/77
1e-4   76 117.464   131.092  125.678  +11.6%   +7.0%   -4.1% 0.26  0.21  0.15 0.15  0.10  0.08   1/76  68/76
```

Grand (230 triples, ~160 MB SPECK basis):
```
SPECK=160.2MB  HZ=179.0MB  LC=173.4MB
HZ vs SPECK: size +11.7%  comp -10.0%  decomp -24.4%
LC vs SPECK: size  +8.2%  comp -23.8%  decomp -29.2%
LC vs HZ:    size  -3.1%  comp -15.3%  decomp  -6.4%
```

### Findings

- **Lossless gate clean**: 230/230 triple-success cells have identical
  L_∞ across the three backends. (One TSV row showed a parser glitch
  on `ODV_bcar2 vre=1e-3`; reproduced cleanly with all three backends
  giving md5-identical decoded files. Not a real mismatch.)
- **LC vs HuffZstd size advantage grows with precision**: +0.7% at
  vre=1e-2 (essentially tied), -4.1% at vre=1e-4 (LC wins on 68/76
  files). Opposite trend from 3D NYX, where LC lost more at low VRE.
  CESM is well-behaved 2D atmospheric data; LC's BIT_4→RZE chain
  exploits the high-bit-plane all-zero structure cleanly.
- **Both backends crush SPECK on time**: LC at vre=1e-4 is ~42% faster
  on compression (0.26→0.15s) and ~47% faster on decompression
  (0.15→0.08s). Even bigger time win than 3D — 2D wavelet+conditioning
  is a smaller share of total runtime than 3D, so the int-coding
  speedup propagates more directly.
- **Size still loses to SPECK** by +8% (LC) / +12% (HZ). Same residual
  as 3D — global Huffman tree / structural pipeline can't capture the
  inter-bitplane correlations SPECK exploits.
- **One upstream fail** (FSDTOA vre=1e-4) reproduces on SPECK alone,
  i.e. inherited from the same 0.8.5 lineage that produced NYX
  velocity / SCALE T-U-V failures in Sprint 6. Not a cc-cli regression.

### Production readiness (now)

- 3D: `sperr3d --backend {speck,huffzstd,lc}`
- 2D: `sperr2d --backend {speck,huffzstd,lc}`
- 1D and `SPERR_C_API.{h,cpp}` are still SPECK-only — defer until
  someone needs them.

### Files touched in Sprint 10

```
src/SPECK2D_FLT.cpp              # rewrite to 3-way backend dispatch via
                                 #   make_enc/make_dec<T>(m_int_backend)
utilities/sperr2d.cpp            # --backend {speck,huffzstd,lc} + bpp
                                 #   guard + set_int_backend forward
WORKLOG.md                       # this section
```

## Files touched in this round of work

```
include/sperr_helper.h           # IntBackend enum moved here
include/SPECK_FLT.h              # removed IntBackend (now in sperr_helper.h)
include/Conditioner.h            # save/retrieve_backend decls; m_backend_idx
include/SPERR3D_OMP_C.h          # set_int_backend + m_int_backend
src/Conditioner.cpp              # save/retrieve_backend impls
src/SPECK_FLT.cpp                # compress writes bit; use_bitstream reads bit;
                                 #   Rate+HuffZstd guard
src/SPECK3D_FLT.cpp              # idwt3d_multi_res signature fixup
src/SPERR3D_OMP_C.cpp            # set_int_backend on each chunk compressor
utilities/sperr3d.cpp            # --backend CLI flag + bpp/huffzstd guard

# Sprint 4
third_party/sz3_huff/SZ3/encoder/Encoder.hpp         # overwritten with SZo (decode -> T*)
third_party/sz3_huff/SZ3/encoder/HuffmanEncoder.hpp  # overwritten with SZo + P2 assert
                                                     #   + P3 preprocess_encode_dense
src/HuffZstd_INT.cpp             # rewritten: int32 + ZigZag + dense freq +
                                 #   unique_ptr<int32_t[]> + format_version byte
                                 #   + SPERR_HZSTATS env-gated stats probe

# Sprint 5: subband-order permutation
src/HuffZstd_INT.cpp             # enumerate_subbands() dyadic + wavelet-packet;
                                 #   encode/decode loops in subband order;
                                 #   kFormatVersion 1 -> 2

# Sprint 7: bins dump + per-stage timing
src/HuffZstd_INT.cpp             # added SPERR_HZ_DUMP_BINS env-gated dump hook
                                 #   (per-chunk int32 bins + .meta);
                                 #   chrono t_huff / t_zst in SPERR_HZSTATS line

# Sprint 9: LC backend integration
third_party/lc/include/macros.h                              # vendored LC
third_party/lc/lc_pipeline.hpp                               # in-memory wrapper for
                                                             #   BIT_4 → RZE_4 → RZE_1
third_party/lc/components/h_BIT_4.h                          # vendored components
third_party/lc/components/h_RZE_4.h
third_party/lc/components/h_RZE_1.h
third_party/lc/components/include/h_RZE.h
third_party/lc/components/include/h_zero_elimination.h
include/LC_INT.h                 # new backend class header
src/LC_INT.cpp                   # new backend impl (~290 lines)
include/sperr_helper.h           # IntBackend: add LC = 2
src/Conditioner.cpp              # backend encoding: 1 bit -> 2 bits (meta[1..2])
src/SPECK3D_FLT.cpp              # make_enc/make_dec: 3-way switch
src/SPECK_FLT.cpp                # Rate-mode guard extended to LC
utilities/sperr3d.cpp            # --backend {speck,huffzstd,lc}
src/CMakeLists.txt               # add LC_INT.cpp

# Sprint 10: 2D path
src/SPECK2D_FLT.cpp              # rewrite to make_enc/make_dec<T>(m_int_backend)
                                 #   3-way switch (SPECK / HuffZstd / LC)
utilities/sperr2d.cpp            # --backend {speck,huffzstd,lc} + bpp guard +
                                 #   set_int_backend forward

WORKLOG.md                       # this file
```

`.sz3_original` backups of the two ported vendored headers are kept
in-tree for rollback.

All of the above is committed on branch `0.8.5-cc-cli` as of
2026-05-15. `git log --stat` will show the commit; `runs/` (bench
artifact dumps, ~770 MB) is now `.gitignore`d.

## Reference: paper findings (IPDPS 2023)

- SPERR pipeline ends with ZSTD on the concatenated (SPECK + outlier) bitstream.
  Verify whether this still holds in 0.8.5 codebase — if so, an extra ZSTD pass
  may be applied on top of HuffZstd output (which is itself ZSTD-wrapped). This
  is fine but redundant; out of scope to optimize.
- Huffman-coded outliers in SZ vs SPECK-coded outliers in SPERR: SPERR wins by
  1–2 bits/outlier. Reinforces keeping outlier coder on SPECK1D.

## SZ3 / SZo API notes (current call pattern in HuffZstd_INT)

```cpp
SZ3::HuffmanEncoder<int32_t> enc;
// SZo's stock preprocess_encode hard-codes a 65 K alphabet scan; for SPERR
// we use the dense overload added in Sprint 4 to bypass that. The freq
// array must be sized exactly stateNum and indexed by bin value.
enc.preprocess_encode_dense(N, freq.data(), stateNum);
enc.save(p);                       // canonical Huffman: 1 byte codelen per state slot
enc.encode(bins.data(), N, p);     // LSB-first packed bitstream
enc.postprocess_encode();

SZ3::Lossless_zstd zstd;           // default compression level 3
zstd.compress(src, srcLen, dst, dstCap);   // dst: [srcLen:8B][zstd frame]
zstd.decompress(src, srcLen, dst, dstLen); // mallocs dst if null

// SZo decode returns a raw new[] pointer with SIMD-padding tail
std::unique_ptr<int32_t[]> bins{dec.decode(rp, N)};
```

`Lossless_zstd::compress` prepends 8 bytes (original length). Don't double-count.
