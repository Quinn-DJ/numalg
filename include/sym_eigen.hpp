#ifndef NUMALG_SYM_EIGEN_HPP
#define NUMALG_SYM_EIGEN_HPP

#include "matrix.hpp"
#include <vector>

namespace numalg {

// Givens rotation
struct Givens {
    double c;  // cos(θ)
    double s;  // sin(θ)
};

// Result of symmetric spectral decomposition
struct SymEigenResult {
    Matrix T;           // Tridiagonal → diagonal (eigenvalues on diagonal)
    Matrix Q;           // Orthogonal, T = Q^T * A * Q
    std::vector<double> eigenvalues;
    bool converged;
    std::size_t iterations;
};

// Compute Givens rotation to zero out z in [x, z]^T
Givens givens(double x, double z);

// Algorithm 7.2.1: Reduce symmetric A to tridiagonal form T
// T = Q^T * A * Q
void tridiagonalReduction(const Matrix& A, Matrix& T, Matrix& Q);

// Algorithm 7.2.2: One implicit symmetric QR step with Wilkinson shift
// on the active tridiagonal block T[first..first+blk-1]
void symQRStep(Matrix& T, std::size_t first, std::size_t blk, Matrix& Q);

// Algorithm 7.2.3: Full spectral decomposition of symmetric A
// Returns eigenvalues and eigenvectors
SymEigenResult symmetricEigensolver(const Matrix& A,
                                    double tol = 1e-14,
                                    std::size_t maxIter = 300);

}  // namespace numalg

#endif  // NUMALG_SYM_EIGEN_HPP
