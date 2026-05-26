#ifndef NUMALG_POWER_METHOD_HPP
#define NUMALG_POWER_METHOD_HPP

#include "matrix.hpp"
#include <vector>
#include <complex>

namespace numalg {

// Build the companion matrix of a monic polynomial.
// Given coefficients c = [c_0, c_1, ..., c_{n-1}] where
//   p(x) = x^n + c_{n-1} x^{n-1} + ... + c_1 x + c_0,
// returns the n×n Frobenius companion matrix:
//   C = [0 0 ... 0  -c_0]
//       [1 0 ... 0  -c_1]
//       [0 1 ... 0  -c_2]
//       [        ...     ]
//       [0 0 ... 1  -c_{n-1}]
Matrix companionMatrix(const std::vector<double>& coeffs);

// Result of the power method iteration
struct PowerMethodResult {
    double eigenvalue;         // dominant eigenvalue (largest magnitude)
    Matrix eigenvector;        // corresponding eigenvector (unit 2-norm)
    std::size_t iterations;    // number of iterations performed
    bool converged;            // whether tolerance was met
    std::vector<double> history;  // |λ_k - λ_{k-1}| history
};

// Power method: find the dominant eigenvalue of matrix A.
//   v0 - initial vector (optional, random if empty)
//   maxIter - maximum iterations
//   tol - convergence tolerance on |λ_k - λ_{k-1}|
PowerMethodResult powerMethod(const Matrix& A,
                              std::size_t maxIter = 1000,
                              double tol = 1e-10,
                              const Matrix& v0 = Matrix(0));

// Convenience: find the largest-magnitude root of a monic polynomial
// using the companion matrix + power method.
PowerMethodResult dominantRoot(const std::vector<double>& coeffs,
                               std::size_t maxIter = 1000,
                               double tol = 1e-10);

}  // namespace numalg

#endif  // NUMALG_POWER_METHOD_HPP
