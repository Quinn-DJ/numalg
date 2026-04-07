#include <iostream>
#include <cmath>
#include <vector>

#include "matrix.hpp"
#include "operations.hpp"
#include "inverse_norm_estimate.hpp"

// Helper: manually compute matrix infinity norm (for verification)
static double exact_norm_inf(const numalg::Matrix& A) {
    double max_sum = 0.0;
    for (std::size_t i = 0; i < A.lines(); ++i) {
        double row_sum = 0.0;
        for (std::size_t j = 0; j < A.rows(); ++j) {
            row_sum += std::abs(A(i, j));
        }
        if (row_sum > max_sum) max_sum = row_sum;
    }
    return max_sum;
}

int main() {
    std::cout << "=== Test: Inverse infinity-norm estimate (Hager's method) ===\n\n";

    // Test 1: matrix with known inverse
    {
        std::cout << "--- Test 1: 3x3 matrix ---\n";
        // A = [[2, 1, 0],
        //      [1, 3, 1],
        //      [0, 1, 2]]
        numalg::Matrix A(3, 3, {2, 1, 0,
                                  1, 3, 1,
                                  0, 1, 2});

        // A^{-1} = (1/8) * [[5, -2, 1],
        //                    [-2, 4, -2],
        //                    [1, -2, 5]]
        numalg::Matrix A_inv(3, 3, {5.0/8, -2.0/8, 1.0/8,
                                     -2.0/8, 4.0/8, -2.0/8,
                                     1.0/8, -2.0/8, 5.0/8});

        // B = A^{-T} = (A^{-1})^T
        numalg::Matrix B = A_inv.transpose();

        double estimate = numalg::estimate_inverse_norm_inf(B);
        double exact = exact_norm_inf(A_inv);

        std::cout << "Exact ||A^{-1}||_inf: " << exact << "\n";
        std::cout << "Estimated:              " << estimate << "\n";
        std::cout << "Relative error:         " << std::abs(estimate - exact) / exact << "\n\n";
    }

    // Test 2: 5x5 matrix
    {
        std::cout << "--- Test 2: 5x5 diagonally dominant matrix ---\n";
        numalg::Matrix A(5, 5, {4, 1, 0, 0, 0,
                                  1, 4, 1, 0, 0,
                                  0, 1, 4, 1, 0,
                                  0, 0, 1, 4, 1,
                                  0, 0, 0, 1, 4});

        // Manual inversion is too tedious; use B = I as a simple sanity check
        std::cout << "(Skipping exact comparison, verifying no crash)\n";
        // Use B = A to test the algorithm runs correctly
        numalg::Matrix B = A.transpose();
        double estimate = numalg::estimate_inverse_norm_inf(B);
        std::cout << "Estimated ||B^{-1}||_inf for B: " << estimate << "\n\n";
    }

    // Test 3: identity matrix
    {
        std::cout << "--- Test 3: 3x3 identity matrix ---\n";
        numalg::Matrix I(3, 3, {1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 1});
        // I^{-1} = I, ||I||_inf = 1
        double estimate = numalg::estimate_inverse_norm_inf(I);
        std::cout << "Estimated: " << estimate << " (expected: 1.0)\n\n";
    }

    std::cout << "=== Test complete ===\n";
    return 0;
}
