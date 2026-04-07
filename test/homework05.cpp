/**
 * Homework 5 - Accuracy estimation of computed solutions
 *
 * Problem 1: Estimate the infinity-norm condition number of 5~20 order Hilbert matrices
 *   cond(A)_\infty = ||A||_\infty * ||A^{-1}||_\infty
 *   where ||A^{-1}||_\infty is estimated via Hager's method
 *
 * Problem 2: For special lower-triangular matrix A_n (n=5~30), generate random x,
 *   compute b=Ax, solve with partial-pivot Gauss elimination,
 *   estimate relative error of the computed solution and compare with true relative error
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>

#include "matrix.hpp"
#include "operations.hpp"
#include "gauss.hpp"
#include "inverse_norm_estimate.hpp"

// Helper: construct Hilbert matrix H(i,j) = 1/(i+j+1)
numalg::Matrix hilbert(std::size_t n) {
    numalg::Matrix H(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            H(i, j) = 1.0 / (i + j + 1.0);
        }
    }
    return H;
}

// Helper: construct A_n matrix from Problem 2
// Diagonal = 1, column n = 1, strict lower triangle = -1, rest = 0
numalg::Matrix make_An(std::size_t n) {
    numalg::Matrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) {
                A(i, j) = 1.0;
            } else if (j == n - 1) {
                A(i, j) = 1.0;
            } else if (i > j) {
                A(i, j) = -1.0;
            }
            // else: 0 (already initialized)
        }
    }
    return A;
}

// Helper: random n-vector with components in [-1, 1]
numalg::Matrix random_vector(std::size_t n) {
    numalg::Matrix x(n, 1);
    for (std::size_t i = 0; i < n; ++i) {
        x(i, 0) = 2.0 * rand() / RAND_MAX - 1.0;
    }
    return x;
}

// Helper: vector infinity norm
double vec_inf_norm(const numalg::Matrix& v) {
    double mx = 0.0;
    for (std::size_t i = 0; i < v.lines(); ++i) {
        mx = std::max(mx, std::abs(v(i, 0)));
    }
    return mx;
}

int main() {
    srand(static_cast<unsigned int>(time(nullptr)));
    std::cout << std::scientific << std::setprecision(6);

    // =============================================
    // Problem 1: Infinity-norm condition number of Hilbert matrix
    // =============================================
    // cond(A)_\infty = ||A||_\infty * ||A^{-1}||_\infty
    // Let B = A^{-T}, then estimate_inverse_norm_inf(B) = ||A^{-1}||_\infty
    std::cout << "===== Problem 1: Hilbert matrix condition number (inf-norm) =====\n";
    std::cout << "  n       ||H||_inf     ||H^{-1}||_inf  cond(H)_inf\n";
    for (std::size_t n = 5; n <= 20; ++n) {
        numalg::Matrix H = hilbert(n);

        // ||H||_\infty: use Matrix::normInf()
        double norm_H = H.normInf();

        // ||H^{-1}||_\infty: compute inverse then take inf-norm
        // H^{-T} = (H^{-1})^T; since H is symmetric, H^{-1} is also symmetric
        numalg::Matrix H_inv = gaussInverse(H);
        double norm_H_inv = H_inv.normInf();

        double cond = norm_H * norm_H_inv;
        std::cout << "  " << std::setw(2) << n << "   "
                  << std::setw(12) << norm_H << "  "
                  << std::setw(14) << norm_H_inv << "  "
                  << std::setw(14) << cond << "\n";
    }

    std::cout << "\n";

    // =============================================
    // Problem 2: Accuracy estimation for special matrix A_n
    // =============================================
    // For n=5~30:
    //   1. Generate random x
    //   2. Compute b = A_n * x
    //   3. Solve with partial-pivot Gauss elimination -> x_tilde
    //   4. Estimate relative error:
    //      - v_hat = estimate_inverse_norm_inf(A_n^T)  (i.e. ||A^{-1}||_\infty)
    //      - mu = ||A||_\infty
    //      - gamma = ||r||_\infty, r = b - A*x_tilde
    //      - beta = ||b||_\infty
    //      - rho = v_hat * mu * gamma / beta  (estimated relative error)
    //   5. True relative error: ||x - x_tilde||_\infty / ||x||_\infty
    std::cout << "===== Problem 2: A_n computed solution accuracy =====\n";
    std::cout << "   n   Est. rel. error  True rel. error\n";
    for (std::size_t n = 5; n <= 30; ++n) {
        numalg::Matrix A = make_An(n);

        // Generate random true solution x
        numalg::Matrix x = random_vector(n);

        // Compute b = A * x
        numalg::Matrix b = A * x;

        // Solve with partial-pivot Gauss elimination
        numalg::Matrix x_tilde = PgaussSolve(A, b);

        // Residual r = b - A * x_tilde
        numalg::Matrix r = b - A * x_tilde;

        // Various norms
        double v_hat = numalg::estimate_inverse_norm_inf(A.transpose());
        double mu = A.normInf();
        double gamma = vec_inf_norm(r);
        double beta = vec_inf_norm(b);

        // Estimated relative error
        double rho = (beta != 0.0) ? v_hat * mu * gamma / beta : 0.0;

        // True relative error
        double err = vec_inf_norm(x - x_tilde) / vec_inf_norm(x);

        std::cout << "  " << std::setw(2) << n << "   "
                  << std::setw(14) << rho << "  "
                  << std::setw(14) << err << "\n";
    }

    return 0;
}
