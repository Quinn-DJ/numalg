#include "sym_eigen.hpp"
#include "hessenberg.hpp"
#include "operations.hpp"
#include <cmath>
#include <limits>
#include <algorithm>

namespace numalg {

constexpr double EPS_SYM = std::numeric_limits<double>::epsilon();

// ======================================================================
// Givens rotation
// ======================================================================

Givens givens(double x, double z) {
    if (z == 0.0) return {1.0, 0.0};
    double r = std::sqrt(x * x + z * z);
    return {x / r, z / r};
}

// ======================================================================
// Algorithm 7.2.1: Householder tridiagonalization
// ======================================================================
void tridiagonalReduction(const Matrix& A, Matrix& T, Matrix& Q) {
    std::size_t n = A.lines();
    T = A;
    Q = Matrix(n, n);
    for (std::size_t i = 0; i < n; ++i) Q(i, i) = 1.0;

    for (std::size_t k = 0; k + 2 < n; ++k) {
        std::size_t m = n - k - 1;  // size of the subvector

        // Extract x = T(k+1 : n-1, k)
        std::vector<double> x(m);
        for (std::size_t i = 0; i < m; ++i)
            x[i] = T(k + 1 + i, k);

        auto [v, beta] = householder(x);
        if (beta == 0.0) continue;

        // u = β * T(k+1:n, k+1:n) * v  (size m)
        std::vector<double> u(m, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < m; ++j)
                sum += T(k + 1 + i, k + 1 + j) * v[j];
            u[i] = beta * sum;
        }

        // w = u - (β/2 * u^T * v) * v  (size m)
        double utv = 0.0;
        for (std::size_t i = 0; i < m; ++i)
            utv += u[i] * v[i];
        double factor = beta * utv / 2.0;

        std::vector<double> w(m);
        for (std::size_t i = 0; i < m; ++i)
            w[i] = u[i] - factor * v[i];

        // Update: T(k+1:n, k+1:n) -= v*w^T + w*v^T
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t j = 0; j < m; ++j) {
                T(k + 1 + i, k + 1 + j) -= v[i] * w[j] + w[i] * v[j];
            }
        }

        // Set tridiagonal element: ||x||_2
        double nrm = 0.0;
        for (std::size_t i = 0; i < m; ++i)
            nrm += x[i] * x[i];
        nrm = std::sqrt(nrm);

        T(k + 1, k) = nrm;
        T(k, k + 1) = nrm;

        // Zero out elements below subdiagonal in column k
        for (std::size_t i = 2; i < m; ++i) {
            T(k + 1 + i, k) = 0.0;
            T(k, k + 1 + i) = 0.0;
        }

        // Update Q = Q * P  where P = I - β*v*v^T
        // Q(:, k+1:n) = Q(:, k+1:n) - p * v^T  where p = β * Q(:, k+1:n) * v
        std::vector<double> p(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < m; ++j)
                sum += Q(i, k + 1 + j) * v[j];
            p[i] = beta * sum;
        }
        for (std::size_t i = 0; i < n; ++i) {
            double pi = p[i];
            for (std::size_t j = 0; j < m; ++j)
                Q(i, k + 1 + j) -= pi * v[j];
        }
    }

    // Symmetrize: ensure strict tridiagonal by zeroing everything beyond
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i > j + 1 || j > i + 1) {
                T(i, j) = 0.0;
            }
        }
    }
}

// ======================================================================
// Algorithm 7.2.2: Implicit symmetric QR step with Wilkinson shift
// ======================================================================
void symQRStep(Matrix& T, std::size_t first, std::size_t blk, Matrix& Q) {
    if (blk <= 1) return;

    std::size_t n = blk;
    std::size_t a = first;

    // Wilkinson shift
    double d = (T(a + n - 2, a + n - 2) - T(a + n - 1, a + n - 1)) / 2.0;
    double mu;
    if (d == 0.0) {
        mu = T(a + n - 1, a + n - 1) - std::abs(T(a + n - 1, a + n - 2));
    } else {
        double sign = (d > 0.0) ? 1.0 : -1.0;
        mu = T(a + n - 1, a + n - 1)
             - T(a + n - 1, a + n - 2) * T(a + n - 1, a + n - 2)
               / (d + sign * std::sqrt(d * d + T(a + n - 1, a + n - 2) * T(a + n - 1, a + n - 2)));
    }

    // Initial bulge from (T - μI) e_1
    double x = T(a, a) - mu;
    double z = T(a + 1, a);

    // Givens chase
    for (std::size_t k = 0; k + 1 < n; ++k) {
        auto [c, s] = givens(x, z);

        // Apply G_k * T * G_k^T
        // Left: G_k * T (rows k, k+1)
        for (std::size_t j = 0; j < n; ++j) {
            double t1 = T(a + k, a + j);
            double t2 = T(a + k + 1, a + j);
            T(a + k, a + j) = c * t1 + s * t2;
            T(a + k + 1, a + j) = -s * t1 + c * t2;
        }

        // Right: (G_k*T) * G_k^T (columns k, k+1)
        for (std::size_t i = 0; i < n; ++i) {
            double t1 = T(a + i, a + k);
            double t2 = T(a + i, a + k + 1);
            T(a + i, a + k) = c * t1 + s * t2;
            T(a + i, a + k + 1) = -s * t1 + c * t2;
        }

        // Update Q = Q * G_k^T
        std::size_t mQ = Q.lines();
        for (std::size_t i = 0; i < mQ; ++i) {
            double q1 = Q(i, a + k);
            double q2 = Q(i, a + k + 1);
            Q(i, a + k) = c * q1 + s * q2;
            Q(i, a + k + 1) = -s * q1 + c * q2;
        }

        // Get next bulge
        if (k + 2 < n) {
            x = T(a + k + 1, a + k);
            z = T(a + k + 2, a + k);
        }
    }

    // Clean up numerical noise
    for (std::size_t i = 0; i + 2 < n; ++i) {
        T(a + i + 2, a + i) = 0.0;
        T(a + i, a + i + 2) = 0.0;
    }
}

// ======================================================================
// Algorithm 7.2.3: Full spectral decomposition
// ======================================================================
SymEigenResult symmetricEigensolver(const Matrix& A, double /*tol*/, std::size_t maxIter) {
    std::size_t n = A.lines();

    Matrix T, Q;
    tridiagonalReduction(A, T, Q);

    std::size_t active = n;
    std::size_t iter = 0;

    while (active > 1 && iter < maxIter) {
        iter++;

        // Zero out negligible off-diagonal elements
        for (std::size_t i = active - 1; i > 0; --i) {
            double thresh = (std::abs(T(i - 1, i - 1)) + std::abs(T(i, i))) * EPS_SYM * 100;
            if (std::abs(T(i, i - 1)) <= thresh) {
                T(i, i - 1) = 0.0;
                T(i - 1, i) = 0.0;
            }
        }

        // Find top of irreducible tridiagonal block
        std::size_t first = 0;
        for (std::size_t i = active - 1; i > 0; --i) {
            if (T(i, i - 1) == 0.0) {
                first = i;
                break;
            }
        }

        std::size_t blk = active - first;
        if (blk <= 1) {
            active = first;
            continue;
        }

        symQRStep(T, first, blk, Q);
    }

    // Extract eigenvalues (diagonal of T)
    std::vector<double> eigenvalues;
    for (std::size_t i = 0; i < n; ++i)
        eigenvalues.push_back(T(i, i));

    SymEigenResult result;
    result.T = T;
    result.Q = Q;
    result.eigenvalues = eigenvalues;
    result.converged = (active <= 1 || iter < maxIter);
    result.iterations = iter;
    return result;
}

}  // namespace numalg
