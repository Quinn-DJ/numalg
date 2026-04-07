#ifndef NUMALG_INVERSE_NORM_ESTIMATE_HPP
#define NUMALG_INVERSE_NORM_ESTIMATE_HPP

#include "matrix.hpp"

namespace numalg {

// Estimate ||B||_1 (when B = A^{-T}, this gives ||A^{-1}||_\infty)
// Based on Hager's method (optimization approach)
double estimate_inverse_norm_inf(const Matrix& B);

}  // namespace numalg

#endif  // NUMALG_INVERSE_NORM_ESTIMATE_HPP
