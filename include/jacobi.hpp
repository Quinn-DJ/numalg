#ifndef NUMALG_JACOBI_HPP
#define NUMALG_JACOBI_HPP

#include "matrix.hpp"
#include <vector>

namespace numalg {

// Result of the Jacobi eigenvalue decomposition
struct JacobiResult {
    Matrix eigenvalues;    // diagonal matrix of eigenvalues
    Matrix Q;             // orthogonal matrix, eigenvalues = Q^T * A * Q
    std::vector<double> eigenvalues_list;  // eigenvalue vector
    std::size_t sweeps;   // number of sweeps performed
    bool converged;       // whether tolerance was met
};

/**
 * @brief Classic Jacobi algorithm
 *
 * At each iteration, find the largest off-diagonal element |A(p,q)|,
 * then apply a Givens rotation J(p,q,θ) to zero it out.
 * Iterates until the off-diagonal Frobenius norm falls below tol.
 *
 * @param A symmetric matrix
 * @param tol convergence tolerance
 * @param maxSweeps maximum number of sweeps (a sweep means zeroing all n(n-1)/2 rotations)
 */
JacobiResult classicJacobi(const Matrix& A,
                           double tol = 1e-14,
                           std::size_t maxSweeps = 100);

/**
 * @brief Cyclic Jacobi algorithm
 *
 * Sweeps through the upper triangle in row-wise order (p,q) and
 * zeros each off-diagonal entry with a Givens rotation.
 * One sweep = n(n-1)/2 rotations.
 *
 * @param A symmetric matrix
 * @param tol convergence tolerance
 * @param maxSweeps maximum number of sweeps
 */
JacobiResult cyclicJacobi(const Matrix& A,
                          double tol = 1e-14,
                          std::size_t maxSweeps = 100);

/**
 * @brief Threshold Jacobi algorithm
 *
 * Like cyclic Jacobi, but only zeros off-diagonal entries whose magnitude
 * exceeds a threshold that is reduced each sweep.
 *
 * @param A symmetric matrix
 * @param tol convergence tolerance
 * @param maxSweeps maximum number of sweeps
 */
JacobiResult thresholdJacobi(const Matrix& A,
                             double tol = 1e-14,
                             std::size_t maxSweeps = 100);

}  // namespace numalg

#endif  // NUMALG_JACOBI_HPP
