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

// ---- CG helper: dot product ----
static double cgDot(const numalg::Matrix& a, const numalg::Matrix& b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.lines(); ++i)
        sum += a(i) * b(i);
    return sum;
}

// ---- CG helper: 2-norm of a vector ----
static double cgNorm2(const numalg::Matrix& x) {
    return std::sqrt(cgDot(x, x));
}

IterResult cgSolve(const Matrix& A, const Matrix& b,
                   double tol, std::size_t maxIter) {
    if (A.lines() != A.rows())
        throw std::invalid_argument("A must be square for CG");
    if (A.lines() != b.lines())
        throw std::invalid_argument("A and b must have compatible sizes");

    std::size_t n = A.lines();
    Matrix x(n);                   // x = 0 (initial guess)
    Matrix r = b;                  // r = b - A*x = b (since x = 0)
    double rho = cgDot(r, r);      // $\rho = r^{\mathrm{T}} r$
    double b_norm = cgNorm2(b);    // $\|b\|_2$

    if (b_norm < 1e-15)
        return {x, 0, true, {0.0}};

    std::vector<double> residual_history;
    residual_history.reserve(maxIter);

    Matrix p(n);
    Matrix w(n);
    double rho_old = 0.0;  // $\widetilde{\rho}$

    for (std::size_t k = 1; k <= maxIter; ++k) {
        // Check convergence: $\sqrt{\rho} \le \varepsilon \|b\|_2$
        double r_norm = std::sqrt(rho);
        residual_history.push_back(r_norm / b_norm);

        if (r_norm <= tol * b_norm)
            return {x, k - 1, true, std::move(residual_history)};

        if (k == 1) {
            // $p = r$
            p = r;
        } else {
            // $\beta = \rho / \widetilde{\rho}$, $p = r + \beta p$
            double beta = rho / rho_old;
            for (std::size_t i = 0; i < n; ++i)
                p(i) = r(i) + beta * p(i);
        }

        // $w = Ap$
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < n; ++j)
                sum += A(i, j) * p(j);
            w(i) = sum;
        }

        double pTw = cgDot(p, w);   // $p^{\mathrm{T}} w$
        // $\alpha = \rho / p^{\mathrm{T}} w$; if p^T A p ≈ 0, CG has exhausted
        // all search directions — treat as converged (residual is at best attainable)
        if (std::abs(pTw) < std::numeric_limits<double>::min() * 1e3)
            return {x, k - 1, true, std::move(residual_history)};

        double alpha = rho / pTw;

        // $x = x + \alpha p$, $r = r - \alpha w$
        for (std::size_t i = 0; i < n; ++i) {
            x(i) += alpha * p(i);
            r(i) -= alpha * w(i);
        }

        // $\widetilde{\rho} = \rho$, $\rho = r^{\mathrm{T}} r$
        rho_old = rho;
        rho = cgDot(r, r);
    }

    double r_norm = std::sqrt(rho);
    residual_history.push_back(r_norm / b_norm);
    return {x, maxIter, false, std::move(residual_history)};
}

IterResult pcgSolve(const Matrix& A, const Matrix& b,
                    const Matrix& M_diag,
                    double tol, std::size_t maxIter) {
    if (A.lines() != A.rows())
        throw std::invalid_argument("A must be square for PCG");
    if (A.lines() != b.lines())
        throw std::invalid_argument("A and b must have compatible sizes");

    std::size_t n = A.lines();
    Matrix x(n);                   // x = 0 (initial guess)
    Matrix r = b;                  // r = b - A*x = b (since x = 0)
    double b_norm = cgNorm2(b);    // $\|b\|_2$

    if (b_norm < 1e-15)
        return {x, 0, true, {0.0}};

    std::vector<double> residual_history;
    residual_history.reserve(maxIter);

    Matrix z(n);                   // $z = M^{-1} r$
    Matrix p(n);
    Matrix w(n);
    double rho = 0.0;              // $\rho = r^{\mathrm{T}} z$
    double rho_old = 0.0;          // $\widetilde{\rho}$

    for (std::size_t k = 1; k <= maxIter; ++k) {
        double r_norm = cgNorm2(r);  // $\sqrt{r^{\mathrm{T}} r}$
        residual_history.push_back(r_norm / b_norm);

        if (r_norm <= tol * b_norm)
            return {x, k - 1, true, std::move(residual_history)};

        // Solve $Mz = r$: Jacobi preconditioner $z_i = r_i / M_{ii}$
        for (std::size_t i = 0; i < n; ++i)
            z(i) = r(i) / M_diag(i);

        if (k == 1) {
            // $p = z$, $\rho = r^{\mathrm{T}} z$
            p = z;
            rho = cgDot(r, z);
        } else {
            // $\widetilde{\rho} = \rho$, $\rho = r^{\mathrm{T}} z$
            // $\beta = \rho / \widetilde{\rho}$, $p = z + \beta p$
            rho_old = rho;
            rho = cgDot(r, z);
            double beta = rho / rho_old;
            for (std::size_t i = 0; i < n; ++i)
                p(i) = z(i) + beta * p(i);
        }

        // $w = Ap$
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < n; ++j)
                sum += A(i, j) * p(j);
            w(i) = sum;
        }

        double pTw = cgDot(p, w);   // $p^{\mathrm{T}} w$
        if (std::abs(pTw) < std::numeric_limits<double>::min() * 1e3)
            return {x, k - 1, true, std::move(residual_history)};

        double alpha = rho / pTw;   // $\alpha = \rho / p^{\mathrm{T}} w$

        // $x = x + \alpha p$, $r = r - \alpha w$
        for (std::size_t i = 0; i < n; ++i) {
            x(i) += alpha * p(i);
            r(i) -= alpha * w(i);
        }
    }

    double r_norm = cgNorm2(r);
    residual_history.push_back(r_norm / b_norm);
    return {x, maxIter, false, std::move(residual_history)};
}

}  // namespace numalg
