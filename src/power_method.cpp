#include "power_method.hpp"
#include <random>
#include <cmath>

namespace numalg {

// ======================================================================
// Companion matrix construction
// ======================================================================
// For p(x) = x^n + c_{n-1}x^{n-1} + ... + c_1 x + c_0,
// the companion matrix is n×n:
//   C = [0 0 ... 0  -c_0]
//       [1 0 ... 0  -c_1]
//       [0 1 ... 0  -c_2]
//       [        ...     ]
//       [0 0 ... 1  -c_{n-1}]
// ======================================================================

Matrix companionMatrix(const std::vector<double>& coeffs) {
    std::size_t n = coeffs.size();
    if (n < 1) {
        throw std::invalid_argument("companionMatrix: need at least one coefficient");
    }
    // The companion matrix is n×n; the polynomial is degree n (monic, coeff of x^n = 1)
    // coeffs[0] = c_0 (constant term), ..., coeffs[n-1] = c_{n-1}
    // So the matrix is n×n.
    Matrix C(n, n);

    // First row: [0, 0, ..., 0, -c_0]
    C(0, n - 1) = -coeffs[0];

    // Subdiagonal: 1's at (i, i-1) for i = 1, ..., n-1
    for (std::size_t i = 1; i < n; ++i) {
        C(i, i - 1) = 1.0;
        // Last column entries: -c_i for i = 1, ..., n-1
        C(i, n - 1) = -coeffs[i];
    }

    return C;
}

// ======================================================================
// Random vector helper
// ======================================================================

static Matrix randomVector(std::size_t n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    Matrix v(n);
    double norm = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        v(i) = dist(gen);
        norm += v(i) * v(i);
    }
    norm = std::sqrt(norm);
    if (norm > 1e-15) {
        for (std::size_t i = 0; i < n; ++i) v(i) /= norm;
    }
    return v;
}

// ======================================================================
// Power method
// ======================================================================
// Algorithm: Starting from v_0, iterate:
//   w_k = A * v_{k-1}
//   v_k = w_k / ||w_k||_2
//   λ_k = v_k^T * A * v_k
// Stop when |λ_k - λ_{k-1}| < tol.
// ======================================================================

PowerMethodResult powerMethod(const Matrix& A,
                              std::size_t maxIter,
                              double tol,
                              const Matrix& v0) {
    std::size_t n = A.lines();
    if (n != A.rows()) {
        throw std::invalid_argument("powerMethod: matrix must be square");
    }

    Matrix v = (v0.lines() == n) ? v0 : randomVector(n);
    double lambda = 0.0;
    bool converged = false;
    std::size_t k;

    PowerMethodResult result;
    result.history.reserve(maxIter);

    for (k = 0; k < maxIter; ++k) {
        // w = A * v
        Matrix w(n);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < n; ++j) {
                sum += A(i, j) * v(j);
            }
            w(i) = sum;
        }

        // ||w||_2
        double w_norm = 0.0;
        for (std::size_t i = 0; i < n; ++i) w_norm += w(i) * w(i);
        w_norm = std::sqrt(w_norm);

        if (w_norm < 1e-15) {
            // Zero vector — A has a zero eigenvalue;
            // v is already the eigenvector, λ = 0.
            result.eigenvalue = 0.0;
            result.eigenvector = v;
            result.iterations = k;
            result.converged = true;
            result.history.push_back(0.0);
            return result;
        }

        // v = w / ||w||
        for (std::size_t i = 0; i < n; ++i) v(i) = w(i) / w_norm;

        // λ = v^T A v = v^T w (since w = Av and v is now w/||w||,
        // but careful: w = A * v_old, and v = w / ||w||.
        // We need v_new^T A v_new. Since A v = A (w / ||w||) = (A w) / ||w||,
        // we can compute λ = v_new^T (A v_new) = (w^T / ||w||) * (A w / ||w||)
        // = (w^T A w) / ||w||^2.
        // But simpler: just recompute:
        Matrix Av(n);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < n; ++j) sum += A(i, j) * v(j);
            Av(i) = sum;
        }
        double lambda_new = 0.0;
        for (std::size_t i = 0; i < n; ++i) lambda_new += v(i) * Av(i);

        double diff = std::abs(lambda_new - lambda);
        result.history.push_back(diff);

        lambda = lambda_new;

        if (diff < tol) {
            converged = true;
            ++k;  // count this iteration
            break;
        }
    }

    result.eigenvalue = lambda;
    result.eigenvector = v;
    result.iterations = k;
    result.converged = converged;

    return result;
}

// ======================================================================
// Dominant root of a polynomial via companion matrix + power method
// ======================================================================

PowerMethodResult dominantRoot(const std::vector<double>& coeffs,
                               std::size_t maxIter,
                               double tol) {
    Matrix C = companionMatrix(coeffs);
    return powerMethod(C, maxIter, tol);
}

}  // namespace numalg
