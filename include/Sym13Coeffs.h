#pragma once
#include <array>
namespace sym13 {
  constexpr size_t kernel_length = 26;

// Symlet-13 分解低通滤波器（dec_lo）
static const std::array<double, kernel_length> dec_lo = {
  0.00006820325263075319, -0.00003573862364868901, -0.0011360634389281183,
  -0.0001709428585302221,   0.007526225389968105,   0.005296359738725025,
  -0.02021676813338983,    -0.01721164272629905,    0.013862497435849205,
  -0.05975062771794370,   -0.12436246075153011,     0.19770481877117801,
   0.69573915056149642,     0.64456438390118563,     0.27269020736982133,
  -0.05194583810787573,    -0.07576571478927332,     0.01519718516421139,
   0.009652395002718613,   -0.0015754644723494324,  -0.0002771129173878957,
   0.00001600369053734232,  0.00000116423758507660, -0.00000007042986690694402,
  -0.00000000036905373423196, 0.000000000035401737968842
};

// 分解高通滤波器 dec_hi[i] = (-1)^(i) * dec_lo[length-1-i]
static const std::array<double, kernel_length> dec_hi = []{
  std::array<double, kernel_length> a{};
  for(size_t i=0;i<kernel_length;++i)
    a[i] = ((i&1)? 1.0:-1.0) * dec_lo[kernel_length-1-i];
  return a;
}();

// 重构低通 rec_lo 直接是 dec_lo 逆序
static const std::array<double, kernel_length> rec_lo = []{
  std::array<double, kernel_length> a{};
  for(size_t i=0;i<kernel_length;++i)
    a[i] = dec_lo[kernel_length-1-i];
  return a;
}();

// 重构高通 rec_hi 直接是 dec_hi 逆序
static const std::array<double, kernel_length> rec_hi = []{
  std::array<double, kernel_length> a{};
  for(size_t i=0;i<kernel_length;++i)
    a[i] = dec_hi[kernel_length-1-i];
  return a;
}();

} // namespace sym13
