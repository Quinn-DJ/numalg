#include "qr_algorithm.hpp"
#include "hessenberg.hpp"
#include "operations.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>

namespace numalg {

constexpr double EPS = std::numeric_limits<double>::epsilon();

// ======================================================================
// Compute Householder reflector for 3-vector [x, y, z]^T
// ======================================================================
static std::pair<std::vector<double>, double> house3(double x, double y, double z) {
    std::vector<double> v(3, 0.0);
    v[0] = 1.0;
    double beta = 0.0;
    double sigma = y * y + z * z;
    if (sigma == 0.0) return {v, 0.0};
    double mu = std::sqrt(x * x + sigma);
    double v0 = (x <= 0.0) ? (x - mu) : (-sigma / (x + mu));
    beta = 2.0 * v0 * v0 / (sigma + v0 * v0);
    v[1] = y / v0;
    v[2] = z / v0;
    return {v, beta};
}

// ======================================================================
// Compute Householder reflector for 2-vector [x, y]^T
// ======================================================================
static std::pair<std::vector<double>, double> house2(double x, double y) {
    std::vector<double> v(2, 0.0);
    v[0] = 1.0;
    double beta = 0.0;
    if (y == 0.0) return {v, 0.0};
    double mu = std::sqrt(x * x + y * y);
    double v0 = (x <= 0.0) ? (x - mu) : (-y * y / (x + mu));
    v[1] = y / v0;
    beta = 2.0 / (1.0 + v[1] * v[1]);
    return {v, beta};
}

// ======================================================================
// Apply P = I - βvv^T from LEFT: rows [r0..r0+nv-1], columns [c0..c1-1]
// ======================================================================
static void applyLeft(Matrix& H, std::size_t r0, std::size_t c0, std::size_t c1,
                      const std::vector<double>& v, double beta) {
    if (beta == 0.0) return;
    std::size_t nv = v.size(), nc = c1 - c0;
    for (std::size_t j = 0; j < nc; ++j) {
        double sum = 0.0;
        for (std::size_t i = 0; i < nv; ++i) sum += H(r0 + i, c0 + j) * v[i];
        sum *= beta;
        for (std::size_t i = 0; i < nv; ++i) H(r0 + i, c0 + j) -= v[i] * sum;
    }
}

// ======================================================================
// Apply P from RIGHT: rows [r0..r1-1], columns [c0..c0+nv-1]
// ======================================================================
static void applyRight(Matrix& H, std::size_t r0, std::size_t r1,
                       std::size_t c0, const std::vector<double>& v, double beta) {
    if (beta == 0.0) return;
    std::size_t nv = v.size(), nr = r1 - r0;
    for (std::size_t i = 0; i < nr; ++i) {
        double sum = 0.0;
        for (std::size_t j = 0; j < nv; ++j) sum += H(r0 + i, c0 + j) * v[j];
        sum *= beta;
        for (std::size_t j = 0; j < nv; ++j) H(r0 + i, c0 + j) -= sum * v[j];
    }
}

// ======================================================================
// Update Q: Q[:, c0:c0+nv-1] = Q[:, c0:c0+nv-1] * P
// ======================================================================
static void updateQ(Matrix& Q, std::size_t c0,
                    const std::vector<double>& v, double beta) {
    if (beta == 0.0) return;
    std::size_t nv = v.size(), n = Q.lines();
    for (std::size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (std::size_t j = 0; j < nv; ++j) sum += Q(i, c0 + j) * v[j];
        sum *= beta;
        for (std::size_t j = 0; j < nv; ++j) Q(i, c0 + j) -= sum * v[j];
    }
}

// ======================================================================
// Francis double-shift QR step for an irreducible Hessenberg block
// ======================================================================
static void francisQRStep(Matrix& H, std::size_t first, std::size_t blk, Matrix& Q) {
    if (blk <= 1) return;

    // Handle 2×2 block separately — use explicit single-shift QR
    if (blk == 2) {
        // Wilkinson shift: eigenvalue of bottom-right 2×2 closer to d
        // For a symmetric matrix, b=c, disc >= 0
        double a = H(first, first), b = H(first, first + 1);
        double c = H(first + 1, first), d = H(first + 1, first + 1);
        double tr = a + d;
        double det = a * d - b * c;
        double disc = tr * tr - 4.0 * det;
        double e1, e2;
        if (disc >= 0) {
            double sd = std::sqrt(disc);
            e1 = (tr + sd) / 2.0;
            e2 = (tr - sd) / 2.0;
        } else {
            e1 = tr / 2.0;
            e2 = tr / 2.0;
        }
        // Pick eigenvalue closer to d
        double mu = (std::abs(e1 - d) < std::abs(e2 - d)) ? e1 : e2;
        double x = a - mu;
        double y = c;
        auto [v, beta] = house2(x, y);
        if (beta != 0.0) {
            applyLeft(H, first, first, first + 2, v, beta);
            applyRight(H, first, first + 2, first, v, beta);
            updateQ(Q, first, v, beta);
        }
        return;
    }

    // blk >= 3: Francis double-shift step
    std::size_t btm = first + blk - 2;
    double s = H(btm, btm) + H(btm + 1, btm + 1);
    double t = H(btm, btm) * H(btm + 1, btm + 1) - H(btm, btm + 1) * H(btm + 1, btm);

    double x = H(first, first) * H(first, first)
               + H(first, first + 1) * H(first + 1, first)
               - s * H(first, first) + t;
    double y = H(first + 1, first) * (H(first, first) + H(first + 1, first + 1) - s);
    double z = H(first + 1, first) * H(first + 2, first + 1);

    // If bulge is all zeros, use small perturbation
    // (handles companion matrices where eigenvalues lie on the unit circle)
    if (x == 0.0 && y == 0.0 && z == 0.0) {
        double eps = std::numeric_limits<double>::epsilon();
        z = std::abs(H(first + 1, first)) * std::sqrt(eps);
        if (z == 0.0) z = eps;
    }

    for (std::size_t k = 0; k + 2 < blk; ++k) {
        auto [v, beta] = house3(x, y, z);

        // LEFT: rows k..k+2, ALL columns of active block
        applyLeft(H, first + k, first, first + blk, v, beta);

        // RIGHT: rows 0..min(k+4,blk)-1, columns k..k+2
        std::size_t r = std::min(k + 4, blk);
        applyRight(H, first, first + r, first + k, v, beta);

        updateQ(Q, first + k, v, beta);

        if (k + 3 < blk) {
            x = H(first + k + 1, first + k);
            y = H(first + k + 2, first + k);
            z = H(first + k + 3, first + k);
        } else {
            x = H(first + blk - 2, first + blk - 3);
            y = H(first + blk - 1, first + blk - 3);
        }
    }

    // Final 2×2 step
    if (blk >= 2) {
        auto [v, beta] = house2(x, y);
        if (beta != 0.0) {
            std::size_t c0 = first + blk - 2;
            // LEFT: rows blk-2, blk-1, columns blk-3..blk-1 (3 cols)
            // for blk>=3, also cover the bulge column at (blk-1, blk-3)
            std::size_t lc = (blk > 2) ? (c0 - 1) : c0;
            applyLeft(H, c0, lc, first + blk, v, beta);
            applyRight(H, first, first + blk, c0, v, beta);
            updateQ(Q, c0, v, beta);
        }
    }
}

// ======================================================================
// Full real Schur decomposition
// ======================================================================
SchurResult realSchurDecomposition(const Matrix& A, double /*tol*/, std::size_t maxIter) {
    if (A.lines() != A.rows())
        throw std::invalid_argument("realSchurDecomposition: A must be square");

    std::size_t n = A.lines();
    Matrix H, Q;
    hessenbergReduction(A, H, Q);

    std::size_t active = n;
    std::size_t iter = 0;

    while (active > 1 && iter < maxIter) {
        iter++;

        // Zero out negligible subdiagonals
        for (std::size_t i = active - 1; i > 0; --i) {
            double thresh = (std::abs(H(i - 1, i - 1)) + std::abs(H(i, i))) * EPS * 100;
            if (std::abs(H(i, i - 1)) <= thresh)
                H(i, i - 1) = 0.0;
        }

        // Find top of irreducible block
        std::size_t first = 0;
        for (std::size_t i = active - 1; i > 0; --i)
            if (H(i, i - 1) == 0.0) { first = i; break; }

        std::size_t blk = active - first;
        if (blk <= 1) { active = first; continue; }

        francisQRStep(H, first, blk, Q);
    }

    // Extract eigenvalues from quasi-triangular form
    std::vector<double> eigenvalues;
    for (std::size_t i = 0; i < n; ++i) {
        if (i + 1 < n && std::abs(H(i + 1, i)) > EPS * (std::abs(H(i, i)) + std::abs(H(i + 1, i + 1))) * 100) {
            double a = H(i, i), b = H(i, i + 1);
            double c = H(i + 1, i), d = H(i + 1, i + 1);
            double tr = a + d;
            double det = a * d - b * c;
            double disc = tr * tr - 4.0 * det;
            if (disc >= 0) {
                double sd = std::sqrt(disc);
                eigenvalues.push_back((tr - sd) / 2.0);
                eigenvalues.push_back((tr + sd) / 2.0);
            } else {
                eigenvalues.push_back(tr / 2.0);
                eigenvalues.push_back(tr / 2.0);
            }
            i++;
        } else {
            eigenvalues.push_back(H(i, i));
        }
    }

    SchurResult result;
    result.T = H;
    result.Q = Q;
    result.eigenvalues = eigenvalues;
    result.converged = (active <= 1 || iter < maxIter);
    result.iterations = iter;
    return result;
}

}  // namespace numalg
