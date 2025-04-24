#pragma once
#include <array>
namespace sym13 {
  constexpr size_t kernel_length = 26;

// Symlet-13 分解低通滤波器（rec_lo）
static const std::array<double, kernel_length> rec_lo = {
 7.0429866906944016e-005, 3.6905373423196241e-005, -0.0007213643851362283, 0.00041326119884196064, 0.0056748537601224395, -0.0014924472742598532, -0.020749686325515677, 0.017618296880653084, 0.092926030899137119, 0.0088197576704205465, -0.14049009311363403, 0.11023022302137217, 0.64456438390118564, 0.69573915056149638, 0.19770481877117801, -0.12436246075153011, -0.059750627717943698, 0.013862497435849205, -0.017211642726299048, -0.02021676813338983, 0.0052963597387250252, 0.0075262253899680996, -0.00017094285853022211, -0.0011360634389281183, -3.5738623648689009e-005, 6.8203252630753188e-005
};

// 重构低通 dec_lo 直接是 rec_lo 逆序
static const std::array<double, kernel_length> dec_lo = []{
  std::array<double, kernel_length> a{};
  for(size_t i=0;i<kernel_length;++i)
    a[i] = rec_lo[kernel_length-1-i];
  return a;
}();



// 分解高通滤波器 dec_hi[i] = (-1)^(i) * dec_lo[length-1-i]
static const std::array<double, kernel_length> dec_hi = []{
  std::array<double, kernel_length> a{};
  for(size_t i=0;i<kernel_length;++i)
    a[i] = ((i&1)? 1.0:-1.0) * dec_lo[kernel_length-1-i];
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
