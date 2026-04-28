#ifndef NUMALG_QR_HPP
#define NUMALG_QR_HPP

#include "matrix.hpp"
#include <vector>
#include <utility>

namespace numalg {

// Householder transformation: given vector x, compute v and beta such that
// (I - beta * v * v^T) * x = ||x|| * e_1
// Returns {v (with v[0] = 1), beta}
std::pair<std::vector<double>, double> house(const std::vector<double>& x);

// QR decomposition via Householder method
// For A (m x n), computes Q (m x m) and R (m x n) such that A = Q * R
struct QRResult {
    Matrix Q;
    Matrix R;
};

QRResult qrDecomposition(const Matrix& A);

// Solve linear system Ax = b using QR decomposition (A is m x m, square)
Matrix qrSolve(const Matrix& A, const Matrix& b);

// Solve linear least squares min ||Ax - b||_2 using QR decomposition (A is m x n, m >= n)
Matrix qrLeastSquares(const Matrix& A, const Matrix& b);

}  // namespace numalg

#endif  // NUMALG_QR_HPP
