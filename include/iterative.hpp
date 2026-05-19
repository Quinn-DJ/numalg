#ifndef NUMALG_ITERATIVE_HPP
#define NUMALG_ITERATIVE_HPP

#include "matrix.hpp"
#include <vector>

namespace numalg {

struct IterResult {
    Matrix x;                // solution vector
    std::size_t iterations;  // number of iterations performed
    bool converged;          // whether convergence was achieved
    std::vector<double> residual_history;  // ||x^{(k)} - x^{(k-1)}||_inf at each step
};

// Jacobi iteration: x^{(k+1)} = D^{-1} (b - (L+U) x^{(k)})
// Convergence criterion: ||x^{(k+1)} - x^{(k)}||_inf < tol
// Throws if any diagonal element of A is zero.
IterResult jacobiSolve(const Matrix& A, const Matrix& b,
                       double tol = 1e-10, std::size_t maxIter = 10000);

// Gauss-Seidel iteration: uses updated components immediately
// x_i^{(k+1)} = (b_i - sum_{j<i} a_{ij} x_j^{(k+1)} - sum_{j>i} a_{ij} x_j^{(k)}) / a_{ii}
// Throws if any diagonal element of A is zero.
IterResult gaussSeidelSolve(const Matrix& A, const Matrix& b,
                            double tol = 1e-10, std::size_t maxIter = 10000);

// SOR (Successive Over-Relaxation) iteration:
// x_i^{(k+1)} = (1 - omega) * x_i^{(k)} + omega / a_{ii} * (b_i - sum_{j!=i} a_{ij} x_j)
// omega in (0, 2) for convergence; omega=1 reduces to Gauss-Seidel
// Throws if any diagonal element of A is zero.
IterResult sorSolve(const Matrix& A, const Matrix& b,
                    double omega = 1.5,
                    double tol = 1e-10, std::size_t maxIter = 10000);

// CG (Conjugate Gradient) for symmetric positive definite systems
// Implements Algorithm 5.3.1 (实用共轭梯度法)
// Starting from x = 0, iterates until ||r||_2 <= tol * ||b||_2
// residual_history stores ||r||_2 / ||b||_2 at each iteration
IterResult cgSolve(const Matrix& A, const Matrix& b,
                   double tol = 1e-10, std::size_t maxIter = 10000);

// PCG (Preconditioned Conjugate Gradient) with Jacobi preconditioner
// Implements Algorithm 5.4.1 (预优共轭梯度法)
// M_diag[i] = diagonal entry of the preconditioner M
// Solves Mz = r as z[i] = r[i] / M_diag[i] (Jacobi preconditioner)
// residual_history stores ||r||_2 / ||b||_2 at each iteration
IterResult pcgSolve(const Matrix& A, const Matrix& b,
                    const Matrix& M_diag,
                    double tol = 1e-10, std::size_t maxIter = 10000);

}  // namespace numalg

#endif  // NUMALG_ITERATIVE_HPP
