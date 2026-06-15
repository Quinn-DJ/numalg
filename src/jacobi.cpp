#include "jacobi.hpp"
#include "operations.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>

namespace numalg {

namespace {

constexpr double EPS = std::numeric_limits<double>::epsilon();

// Compute the Frobenius off-diagonal norm: sqrt(2 * sum_{i<j} A(i,j)^2)
double offNorm(const Matrix& A) {
    std::size_t n = A.lines();
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            sum += A(i, j) * A(i, j);
        }
    }
    return std::sqrt(2.0 * sum);
}

// Apply one Jacobi rotation J(p,q) to symmetric A and accumulate Q.
// On input, A is symmetric. The rotation zeroes A(p,q) = A(q,p).
// Updates:
//   A ← J^T A J   (J = Givens rotation in plane (p,q) with angle θ)
//   Q ← Q J
void jacobiRotation(Matrix& A, Matrix& Q,
                    std::size_t p, std::size_t q) {
    std::size_t n = A.lines();

    // Compute rotation angle
    double app = A(p, p);
    double aqq = A(q, q);
    double apq = A(p, q);

    if (apq == 0.0) return;

    double tau = (aqq - app) / (2.0 * apq);
    double t;

    if (tau >= 0.0) {
        t = 1.0 / (tau + std::sqrt(1.0 + tau * tau));
    } else {
        t = -1.0 / (-tau + std::sqrt(1.0 + tau * tau));
    }

    double c = 1.0 / std::sqrt(1.0 + t * t);
    double s = t * c;

    // Update diagonal entries
    double app_new = app - t * apq;
    double aqq_new = aqq + t * apq;
    A(p, p) = app_new;
    A(q, q) = aqq_new;
    A(p, q) = 0.0;
    A(q, p) = 0.0;

    // Update off-diagonal entries involving p or q
    for (std::size_t k = 0; k < n; ++k) {
        if (k == p || k == q) continue;

        double apk = A(p, k);
        double aqk = A(q, k);

        A(p, k) = c * apk - s * aqk;
        A(k, p) = A(p, k);

        A(q, k) = s * apk + c * aqk;
        A(k, q) = A(q, k);
    }

    // Accumulate eigenvectors: Q ← Q * J
    // J(p,q) = Givens rotation: J(p,p)=c, J(p,q)=s, J(q,p)=-s, J(q,q)=c
    // (Q * J)(i, p) = Q(i,p)*c + Q(i,q)*(-s) = c*Q(i,p) - s*Q(i,q)
    // (Q * J)(i, q) = Q(i,p)*s  + Q(i,q)*c    = s*Q(i,p) + c*Q(i,q)
    for (std::size_t i = 0; i < n; ++i) {
        double qip = Q(i, p);
        double qiq = Q(i, q);
        Q(i, p) = c * qip - s * qiq;
        Q(i, q) = s * qip + c * qiq;
    }
}

// Build a JacobiResult from the converged A and Q matrices
JacobiResult makeResult(const Matrix& A, const Matrix& Q,
                        std::size_t sweeps, bool converged) {
    JacobiResult result;
    std::size_t n = A.lines();

    // Copy eigenvalues onto the diagonal
    result.eigenvalues = Matrix(n, n);
    result.eigenvalues_list.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        result.eigenvalues(i, i) = A(i, i);
        result.eigenvalues_list[i] = A(i, i);
    }

    result.Q = Q;
    result.sweeps = sweeps;
    result.converged = converged;
    return result;
}

}  // anonymous namespace

// ======================================================================
// Classic Jacobi: always pick the largest |A(p,q)|
// ======================================================================
JacobiResult classicJacobi(const Matrix& A,
                           double tol,
                           std::size_t maxSweeps) {
    std::size_t n = A.lines();

    Matrix Q(n, n);
    for (std::size_t i = 0; i < n; ++i) Q(i, i) = 1.0;

    Matrix Ak = A;  // working copy

    double initialOff = offNorm(Ak);
    double target = tol * initialOff;

    std::size_t sweeps = 0;
    bool converged = false;

    while (sweeps < maxSweeps) {
        double curOff = offNorm(Ak);
        if (curOff <= target) {
            converged = true;
            break;
        }

        // Find the largest off-diagonal |A(p,q)|
        double maxVal = 0.0;
        std::size_t p = 0, q = 1;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                double absVal = std::abs(Ak(i, j));
                if (absVal > maxVal) {
                    maxVal = absVal;
                    p = i;
                    q = j;
                }
            }
        }

        if (maxVal < target) {
            converged = true;
            break;
        }

        jacobiRotation(Ak, Q, p, q);
        sweeps++;
    }

    return makeResult(Ak, Q, sweeps, converged);
}

// ======================================================================
// Cyclic Jacobi: sweep through upper triangle row-wise
// ======================================================================
JacobiResult cyclicJacobi(const Matrix& A,
                          double tol,
                          std::size_t maxSweeps) {
    std::size_t n = A.lines();

    Matrix Q(n, n);
    for (std::size_t i = 0; i < n; ++i) Q(i, i) = 1.0;

    Matrix Ak = A;

    double initialOff = offNorm(Ak);
    double target = tol * initialOff;

    std::size_t sweeps = 0;
    bool converged = false;

    while (sweeps < maxSweeps) {
        double curOff = offNorm(Ak);
        if (curOff <= target) {
            converged = true;
            break;
        }

        // One sweep: rotate through all (p,q) pairs in row-wise order
        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                jacobiRotation(Ak, Q, p, q);
            }
        }

        sweeps++;
    }

    return makeResult(Ak, Q, sweeps, converged);
}

// ======================================================================
// Threshold Jacobi: skip small off-diagonals within each sweep
// ======================================================================
JacobiResult thresholdJacobi(const Matrix& A,
                             double tol,
                             std::size_t maxSweeps) {
    std::size_t n = A.lines();

    Matrix Q(n, n);
    for (std::size_t i = 0; i < n; ++i) Q(i, i) = 1.0;

    Matrix Ak = A;

    double initialOff = offNorm(Ak);
    double target = tol * initialOff;

    std::size_t sweeps = 0;
    bool converged = false;
    double threshold = initialOff / static_cast<double>(n);

    while (sweeps < maxSweeps) {
        double curOff = offNorm(Ak);
        if (curOff <= target) {
            converged = true;
            break;
        }

        // Reduce threshold each sweep
        threshold /= static_cast<double>(n);

        std::size_t rotationsInSweep = 0;

        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                if (std::abs(Ak(p, q)) > threshold) {
                    jacobiRotation(Ak, Q, p, q);
                    rotationsInSweep++;
                }
            }
        }

        // If no rotations were done, the matrix is close enough to diagonal
        if (rotationsInSweep == 0) {
            converged = true;
            break;
        }

        sweeps++;
    }

    return makeResult(Ak, Q, sweeps, converged);
}

}  // namespace numalg
