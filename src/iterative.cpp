#include "iterative.hpp"
#include "matrix.hpp"
#include <cmath>
#include <stdexcept>

namespace numalg {

// Helper: compute ||x - x_prev||_inf
static double infNormDiff(const Matrix& x, const Matrix& x_prev) {
    double d = 0.0;
    for (std::size_t i = 0; i < x.lines(); ++i)
        d = std::max(d, std::abs(x(i) - x_prev(i)));
    return d;
}

// Check basic requirements: square, compatible, non-zero diagonal
static void checkIterInput(const Matrix& A, const Matrix& b) {
    if (A.lines() != A.rows())
        throw std::invalid_argument("matrix A must be square");
    if (A.lines() != b.lines())
        throw std::invalid_argument("A and b must have the same size");
    for (std::size_t i = 0; i < A.lines(); ++i)
        if (std::abs(A(i, i)) < 1e-15)
            throw std::invalid_argument("zero diagonal element at row " + std::to_string(i));
}

IterResult jacobiSolve(const Matrix& A, const Matrix& b,
                       double tol, std::size_t maxIter) {
    checkIterInput(A, b);
    std::size_t n = A.lines();

    Matrix x(n);      // x^{(k+1)}
    Matrix x_prev(n); // x^{(k)}, initialized to zero

    std::vector<double> residual_history;
    residual_history.reserve(maxIter);

    for (std::size_t iter = 0; iter < maxIter; ++iter) {
        for (std::size_t i = 0; i < n; ++i) {
            double sigma = 0.0;
            for (std::size_t j = 0; j < n; ++j) {
                if (j != i) sigma += A(i, j) * x_prev(j);
            }
            x(i) = (b(i) - sigma) / A(i, i);
        }

        double diff = infNormDiff(x, x_prev);
        residual_history.push_back(diff);
        x_prev = x;

        if (diff < tol) {
            return {x, iter + 1, true, std::move(residual_history)};
        }
    }

    return {x, maxIter, false, std::move(residual_history)};
}

IterResult gaussSeidelSolve(const Matrix& A, const Matrix& b,
                            double tol, std::size_t maxIter) {
    checkIterInput(A, b);
    std::size_t n = A.lines();

    Matrix x(n);  // initialized to zero

    std::vector<double> residual_history;
    residual_history.reserve(maxIter);

    for (std::size_t iter = 0; iter < maxIter; ++iter) {
        double diff = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double x_old = x(i);
            double sigma = 0.0;
            for (std::size_t j = 0; j < n; ++j) {
                if (j != i) sigma += A(i, j) * x(j);
            }
            x(i) = (b(i) - sigma) / A(i, i);
            diff = std::max(diff, std::abs(x(i) - x_old));
        }

        residual_history.push_back(diff);

        if (diff < tol) {
            return {x, iter + 1, true, std::move(residual_history)};
        }
    }

    return {x, maxIter, false, std::move(residual_history)};
}

IterResult sorSolve(const Matrix& A, const Matrix& b,
                    double omega, double tol, std::size_t maxIter) {
    if (omega <= 0.0 || omega >= 2.0)
        throw std::invalid_argument("omega must be in (0, 2) for SOR convergence");
    checkIterInput(A, b);
    std::size_t n = A.lines();

    Matrix x(n);  // initialized to zero

    std::vector<double> residual_history;
    residual_history.reserve(maxIter);

    for (std::size_t iter = 0; iter < maxIter; ++iter) {
        double diff = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double x_old = x(i);
            double sigma = 0.0;
            for (std::size_t j = 0; j < n; ++j) {
                if (j != i) sigma += A(i, j) * x(j);
            }
            double gs = (b(i) - sigma) / A(i, i);
            x(i) = (1.0 - omega) * x_old + omega * gs;
            diff = std::max(diff, std::abs(x(i) - x_old));
        }

        residual_history.push_back(diff);

        if (diff < tol) {
            return {x, iter + 1, true, std::move(residual_history)};
        }
    }

    return {x, maxIter, false, std::move(residual_history)};
}

}  // namespace numalg
