#ifndef NUMALG_QR_ALGORITHM_HPP
#define NUMALG_QR_ALGORITHM_HPP

#include "matrix.hpp"
#include <vector>

namespace numalg {

// Result of the real Schur decomposition
struct SchurResult {
    Matrix T;       // Quasi-triangular (real Schur form)
    Matrix Q;       // Orthogonal matrix, T = Q^T * A * Q
    std::vector<double> eigenvalues;  // Real eigenvalues
    bool converged;
    std::size_t iterations;
};

// Apply one implicit double-shift QR step (Algorithm 6.4.2)
// to the active (n × n) submatrix of H starting at (start, start).
// Q accumulates the orthogonal transformations.
void implicitQRStep(Matrix& H, std::size_t start, std::size_t n, Matrix& Q);

// Compute the full real Schur decomposition of a real matrix A.
// A is first reduced to Hessenberg form, then the implicit QR algorithm
// iterates until convergence.
SchurResult realSchurDecomposition(const Matrix& A,
                                   double tol = 1e-14,
                                   std::size_t maxIter = 500);

}  // namespace numalg

#endif  // NUMALG_QR_ALGORITHM_HPP
