#ifndef NUMALG_INVERSE_NORM_ESTIMATE_HPP
#define NUMALG_INVERSE_NORM_ESTIMATE_HPP

#include "matrix.hpp"

namespace numalg {

// 估计 ||B||_1（当 B = A^{-T} 时，即 ||A^{-1}||_∞ 的估计值）
// 基于优化法（Hager's method）
double estimate_inverse_norm_inf(const Matrix& B);

}  // namespace numalg

#endif  // NUMALG_INVERSE_NORM_ESTIMATE_HPP
