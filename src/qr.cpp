#include "qr.hpp"
#include "operations.hpp"
#include "solve.hpp"
#include <cmath>
#include <stdexcept>

namespace numalg {

// Householder transformation
// Given x, compute v and beta such that (I - beta * v * v^T) * x = ||x|| * e_1
std::pair<std::vector<double>, double> house(const std::vector<double>& x) {
    std::size_t n = x.size();

    // eta = ||x||_inf
    double eta = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        eta = std::max(eta, std::abs(x[i]));
    }

    if (eta == 0.0) {
        // x is zero vector
        std::vector<double> v(n, 0.0);
        return {v, 0.0};
    }

    // x = x / eta
    std::vector<double> xp(n);
    for (std::size_t i = 0; i < n; ++i) {
        xp[i] = x[i] / eta;
    }

    // sigma = x(2:n)^T * x(2:n)
    double sigma = 0.0;
    for (std::size_t i = 1; i < n; ++i) {
        sigma += xp[i] * xp[i];
    }

    std::vector<double> v(n);
    for (std::size_t i = 1; i < n; ++i) {
        v[i] = xp[i];
    }

    double beta;
    if (sigma == 0.0) {
        beta = 0.0;
        v[0] = 1.0;
    } else {
        double alpha = std::sqrt(xp[0] * xp[0] + sigma);
        if (xp[0] <= 0.0) {
            v[0] = xp[0] - alpha;
        } else {
            v[0] = -sigma / (xp[0] + alpha);
        }
        beta = 2.0 * v[0] * v[0] / (sigma + v[0] * v[0]);
        // Normalize so that v[0] = 1
        double v0 = v[0];
        for (std::size_t i = 0; i < n; ++i) {
            v[i] /= v0;
        }
    }

    return {v, beta};
}

// QR decomposition via Householder method
QRResult qrDecomposition(const Matrix& A) {
    std::size_t m = A.lines();
    std::size_t n = A.rows();

    // Work on a copy of A
    Matrix R = A;

    // Store beta values for Q reconstruction
    std::vector<double> beta_vals(n, 0.0);

    std::size_t min_mn = std::min(m, n);
    for (std::size_t j = 0; j < min_mn; ++j) {
        std::size_t len = m - j;
        if (len <= 1) break;

        // Extract column A(j:m, j)
        std::vector<double> col(len);
        for (std::size_t i = 0; i < len; ++i) {
            col[i] = R(j + i, j);
        }

        auto [v, beta] = house(col);
        beta_vals[j] = beta;

        // Apply Householder: A(j:m, j:n) = (I - beta * v * v^T) * A(j:m, j:n)
        // This zeros out sub-diagonal of column j and updates columns j:n
        // w^T = beta * v^T * A(j:m, j:n)
        std::vector<double> w(n - j, 0.0);
        for (std::size_t k = 0; k < len; ++k) {
            for (std::size_t l = 0; l < n - j; ++l) {
                w[l] += beta * v[k] * R(j + k, j + l);
            }
        }

        // A(j:m, j:n) = A(j:m, j:n) - v * w^T
        for (std::size_t k = 0; k < len; ++k) {
            for (std::size_t l = 0; l < n - j; ++l) {
                R(j + k, j + l) -= v[k] * w[l];
            }
        }

        // Store v[1:] below the diagonal for Q reconstruction
        // (R(j,j) keeps the diagonal value, sub-diagonal was zeroed by Householder)
        for (std::size_t i = 1; i < len; ++i) {
            R(j + i, j) = v[i];
        }
    }

    // Extract upper triangular part as R
    Matrix R_tri(m, n);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i <= j) {
                R_tri(i, j) = R(i, j);
            }
        }
    }

    // Reconstruct Q = H_1 * H_2 * ... * H_k by applying Householder vectors in reverse
    // Q starts as I_m
    Matrix Q(m, m);
    for (std::size_t i = 0; i < m; ++i) {
        Q(i, i) = 1.0;
    }

    // Build v vectors from stored data
    for (std::size_t j = min_mn; j-- > 0;) {
        std::size_t len = m - j;
        if (len <= 1) continue;

        // Reconstruct v: v[0] = 1, v[1:] = R(j+1:m, j)
        std::vector<double> v(len);
        v[0] = 1.0;
        for (std::size_t i = 1; i < len; ++i) {
            v[i] = R(j + i, j);
        }
        double beta = beta_vals[j];

        // Apply H_j on the left: Q = H_j * Q
        // H_j affects rows j:m across ALL columns
        for (std::size_t col = 0; col < m; ++col) {
            double p = 0.0;
            for (std::size_t k = 0; k < len; ++k) {
                p += v[k] * Q(j + k, col);
            }
            p *= beta;
            for (std::size_t k = 0; k < len; ++k) {
                Q(j + k, col) -= v[k] * p;
            }
        }
    }

    return {Q, R_tri};
}

// Solve linear system Ax = b using QR decomposition (A is square)
Matrix qrSolve(const Matrix& A, const Matrix& b) {
    if (A.lines() != A.rows()) {
        throw std::invalid_argument("matrix A must be square for qrSolve");
    }
    if (A.lines() != b.lines()) {
        throw std::invalid_argument("A and b must have the same number of lines");
    }

    auto [Q, R] = qrDecomposition(A);
    // Ax = b => QRx = b => Rx = Q^T b
    Matrix Qtb = tr(Q) * b;
    return solve_upper_triangular(R, Qtb);
}

// Solve linear least squares min ||Ax - b||_2 using QR decomposition
// A is m x n (m >= n), b is m x 1
Matrix qrLeastSquares(const Matrix& A, const Matrix& b) {
    std::size_t m = A.lines();
    std::size_t n = A.rows();

    if (m < n) {
        throw std::invalid_argument("A must have at least as many lines as rows for least squares");
    }
    if (A.lines() != b.lines()) {
        throw std::invalid_argument("A and b must have the same number of lines");
    }

    auto [Q, R] = qrDecomposition(A);

    // Q is m x m, R is m x n
    // R1 = R(1:n, 1:n), c = Q^T * b (first n components)
    // Solve R1 * x = c
    Matrix Qtb = tr(Q) * b;

    // Extract upper triangular part of R (n x n)
    Matrix R1(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i <= j) {
                R1(i, j) = R(i, j);
            }
        }
    }

    // Extract first n components of Q^T b
    Matrix c(n, b.rows());
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < b.rows(); ++j) {
            c(i, j) = Qtb(i, j);
        }
    }

    return solve_upper_triangular(R1, c);
}

}  // namespace numalg
