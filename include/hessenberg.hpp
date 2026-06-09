#ifndef NUMALG_HESSENBERG_HPP
#define NUMALG_HESSENBERG_HPP

#include "matrix.hpp"
#include <vector>

namespace numalg {

// Householder reflector: compute v and beta such that
// (I - beta * v * v^T) * x = ||x||_2 * e_1
// Returns {v (with v[0] = 1), beta}
struct HouseholderResult {
    std::vector<double> v;
    double beta;
};

HouseholderResult householder(const std::vector<double>& x);

// Reduce A to upper Hessenberg form using Householder reflectors.
// On output:
//   H = Q^T * A * Q   (upper Hessenberg)
//   Q = H_1 * H_2 * ... * H_{n-2}  (orthogonal)
void hessenbergReduction(const Matrix& A, Matrix& H, Matrix& Q);

}  // namespace numalg

#endif  // NUMALG_HESSENBERG_HPP
