#include "hessenberg.hpp"
#include "operations.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace numalg {

// ======================================================================
// Householder reflector
// ======================================================================
// Given vector x (size m), compute v (with v[0] = 1) and beta such that
//   (I - beta * v * v^T) * x = ||x||_2 * e_1
// ======================================================================

HouseholderResult householder(const std::vector<double>& x) {
    std::size_t m = x.size();
    HouseholderResult result;
    result.v.resize(m);

    // eta = ||x||_inf
    double eta = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
        eta = std::max(eta, std::abs(x[i]));
    }

    if (eta == 0.0) {
        // x is zero vector → P = I
        result.v[0] = 1.0;
        result.beta = 0.0;
        return result;
    }

    // Scale x to avoid overflow
    std::vector<double> xp(m);
    for (std::size_t i = 0; i < m; ++i) {
        xp[i] = x[i] / eta;
    }

    // sigma = xp[1:].dot(xp[1:])
    double sigma = 0.0;
    for (std::size_t i = 1; i < m; ++i) {
        sigma += xp[i] * xp[i];
    }

    if (sigma == 0.0) {
        // x is already a multiple of e_1
        result.v[0] = 1.0;
        result.beta = 0.0;
        return result;
    }

    // alpha = ||xp||_2
    double alpha = std::sqrt(xp[0] * xp[0] + sigma);

    if (xp[0] <= 0.0) {
        result.v[0] = xp[0] - alpha;
    } else {
        result.v[0] = -sigma / (xp[0] + alpha);
    }

    for (std::size_t i = 1; i < m; ++i) {
        result.v[i] = xp[i];
    }

    result.beta = 2.0 * result.v[0] * result.v[0] / (sigma + result.v[0] * result.v[0]);

    // Normalize so that v[0] = 1
    double v0 = result.v[0];
    for (std::size_t i = 0; i < m; ++i) {
        result.v[i] /= v0;
    }

    return result;
}

// ======================================================================
// Householder reduction to upper Hessenberg form
// ======================================================================
// Algorithm: for k = 1 : n-2 (1-indexed), i.e. k = 0 : n-3 (0-indexed)
//   [v, beta] = house(A(k+1:n, k))
//   A(k+1:n, k:n) = (I - beta*v*v^T) * A(k+1:n, k:n)
//   A(1:n, k+1:n) = A(1:n, k+1:n) * (I - beta*v*v^T)
// ======================================================================

void hessenbergReduction(const Matrix& A, Matrix& H, Matrix& Q) {
    if (A.lines() != A.rows()) {
        throw std::invalid_argument("hessenbergReduction: A must be square");
    }

    std::size_t n = A.lines();
    H = A;
    Q = Matrix(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        Q(i, i) = 1.0;
    }

    for (std::size_t k = 0; k + 2 < n; ++k) {
        std::size_t m = n - k - 1;  // length of subdiagonal section: rows k+1..n-1

        // Extract x = H(k+1 : n-1, k)
        std::vector<double> x(m);
        for (std::size_t i = 0; i < m; ++i) {
            x[i] = H(k + 1 + i, k);
        }

        auto [v, beta] = householder(x);

        if (beta == 0.0) continue;

        // === Left multiplication: H(k+1:n, k:n) = (I - beta*v*v^T) * H(k+1:n, k:n) ===
        // w^T = beta * v^T * H(k+1:n, k:n)
        // H(k+1:n, k:n) -= v * w^T

        std::size_t rightCols = n - k;  // columns k..n-1
        std::vector<double> w(rightCols, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            double vi = v[i];
            for (std::size_t j = 0; j < rightCols; ++j) {
                w[j] += beta * vi * H(k + 1 + i, k + j);
            }
        }
        for (std::size_t i = 0; i < m; ++i) {
            double vi = v[i];
            for (std::size_t j = 0; j < rightCols; ++j) {
                H(k + 1 + i, k + j) -= vi * w[j];
            }
        }

        // === Right multiplication: H(1:n, k+1:n) = H(1:n, k+1:n) * (I - beta*v*v^T) ===
        // p = beta * H(1:n, k+1:n) * v
        // H(1:n, k+1:n) -= p * v^T

        std::vector<double> p(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < m; ++j) {
                sum += H(i, k + 1 + j) * v[j];
            }
            p[i] = beta * sum;
        }
        for (std::size_t i = 0; i < n; ++i) {
            double pi = p[i];
            for (std::size_t j = 0; j < m; ++j) {
                H(i, k + 1 + j) -= pi * v[j];
            }
        }

        // === Update Q: Q = Q * (I - beta*v*v^T) ===
        // Q affects columns k+1:n similarly
        // q = beta * Q(1:n, k+1:n) * v
        // Q(1:n, k+1:n) -= q * v^T

        std::vector<double> q(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < m; ++j) {
                sum += Q(i, k + 1 + j) * v[j];
            }
            q[i] = beta * sum;
        }
        for (std::size_t i = 0; i < n; ++i) {
            double qi = q[i];
            for (std::size_t j = 0; j < m; ++j) {
                Q(i, k + 1 + j) -= qi * v[j];
            }
        }
    }
}

}  // namespace numalg
