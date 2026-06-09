#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "hessenberg.hpp"
#include "qr_algorithm.hpp"
#include "operations.hpp"

// Check if a matrix is upper Hessenberg (within tolerance)
bool isHessenberg(const numalg::Matrix& H, double tol = 1e-12) {
    std::size_t n = H.lines();
    for (std::size_t i = 2; i < n; ++i)
        for (std::size_t j = 0; j + 2 <= i; ++j)  // j < i-1, i.e. j <= i-2
            if (std::abs(H(i, j)) > tol) return false;
    return true;
}

// Check orthogonality: ||Q^T Q - I||_inf
double orthogonalityError(const numalg::Matrix& Q) {
    std::size_t n = Q.lines();
    double maxErr = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < n; ++k)
                sum += Q(k, i) * Q(k, j);
            if (i == j) sum -= 1.0;
            maxErr = std::max(maxErr, std::abs(sum));
        }
    }
    return maxErr;
}

// Check similarity: ||Q^T A Q - T||_inf
double schurSimilarityError(const numalg::Matrix& A, const numalg::Matrix& Q, const numalg::Matrix& T) {
    numalg::Matrix QtAQ = numalg::tr(Q) * A * Q;
    double maxErr = 0.0;
    for (std::size_t i = 0; i < A.lines(); ++i)
        for (std::size_t j = 0; j < A.lines(); ++j)
            maxErr = std::max(maxErr, std::abs(QtAQ(i, j) - T(i, j)));
    return maxErr;
}

void runTests() {
    std::ofstream out("./output/homework14.txt");
    if (!out.is_open()) {
        std::cerr << "Failed to open output file.\n";
        return;
    }

    out << std::fixed << std::setprecision(15);
    out << "============================================================\n";
    out << "Homework 14: Implicit QR Algorithm for Real Schur Decomposition\n";
    out << "============================================================\n\n";

    // ==================================================================
    // Test 1: Householder reduction to Hessenberg
    // ==================================================================
    {
        out << "--- Test 1: Householder Reduction to Hessenberg Form ---\n\n";

        numalg::Matrix A(5, 5, {
            4.0,  1.0,  2.0,  3.0,  1.0,
            1.0,  3.0,  1.0,  2.0,  2.0,
            2.0,  1.0,  5.0,  1.0,  3.0,
            3.0,  2.0,  1.0,  2.0,  1.0,
            1.0,  2.0,  3.0,  1.0,  4.0
        });

        numalg::Matrix H, Q;
        numalg::hessenbergReduction(A, H, Q);

        out << "  Hessenberg form H =\n";
        for (std::size_t i = 0; i < 5; ++i) {
            out << "    ";
            for (std::size_t j = 0; j < 5; ++j)
                out << std::setw(14) << H(i, j);
            out << "\n";
        }
        out << "\n";

        double orthoErr = orthogonalityError(Q);
        out << "  Orthogonality error ||Q^T Q - I||:        " << std::scientific << orthoErr << "\n";

        double simErr = schurSimilarityError(A, Q, H);
        out << "  Similarity error ||Q^T A Q - H||_inf:     " << std::scientific << simErr << "\n";

        out << "  Is Hessenberg:                            " << (isHessenberg(H) ? "Yes" : "No") << "\n";
        out << "\n";
    }

    // ==================================================================
    // Test 2: 3×3 symmetric matrix
    // ==================================================================
    {
        out << "--- Test 2: 3×3 Symmetric Matrix ---\n\n";

        // A = [3 1 1; 1 3 1; 1 1 3]
        // Eigenvalues: λ = 5, 2, 2
        numalg::Matrix A(3, 3, {
            3.0, 1.0, 1.0,
            1.0, 3.0, 1.0,
            1.0, 1.0, 3.0
        });

        auto result = numalg::realSchurDecomposition(A);
        double orthoErr = orthogonalityError(result.Q);
        double simErr = schurSimilarityError(A, result.Q, result.T);

        out << "  Converged:  " << (result.converged ? "Yes" : "No")
            << "  (iterations: " << result.iterations << ")\n\n";

        out << "  Eigenvalues:\n";
        std::vector<double> expected = {5.0, 2.0, 2.0};
        std::vector<double> computed = result.eigenvalues;
        std::sort(computed.begin(), computed.end(), std::greater<double>());
        std::sort(expected.begin(), expected.end(), std::greater<double>());
        for (std::size_t i = 0; i < 3; ++i)
            out << "    λ[" << i << "] = " << std::setw(18) << computed[i]
                << "  (expected " << expected[i] << ")"
                << "  err = " << std::scientific << std::abs(computed[i] - expected[i]) << "\n";
        out << "\n";

        out << "  Orthogonality error:  " << std::scientific << orthoErr << "\n";
        out << "  Similarity error:     " << std::scientific << simErr << "\n";
        out << "\n";
    }

    // ==================================================================
    // Test 3: 4×4 Poisson matrix
    // ==================================================================
    {
        out << "--- Test 3: 4×4 Poisson (tridiagonal 2,-1) ---\n\n";

        numalg::Matrix A(4, 4, {
            2.0, -1.0,  0.0,  0.0,
            -1.0,  2.0, -1.0,  0.0,
            0.0, -1.0,  2.0, -1.0,
            0.0,  0.0, -1.0,  2.0
        });

        auto result = numalg::realSchurDecomposition(A);
        double orthoErr = orthogonalityError(result.Q);
        double simErr = schurSimilarityError(A, result.Q, result.T);

        out << "  Converged:  " << (result.converged ? "Yes" : "No")
            << "  (iterations: " << result.iterations << ")\n\n";

        // λ_k = 2 - 2cos(kπ/5), k = 1..4
        std::vector<double> expected(4);
        for (int k = 1; k <= 4; ++k)
            expected[k-1] = 2.0 - 2.0 * std::cos(k * M_PI / 5.0);

        out << "  Expected eigenvalues:\n";
        for (std::size_t i = 0; i < 4; ++i)
            out << "    λ[" << i << "] = " << std::setw(18) << expected[i] << "\n";

        std::vector<double> computed = result.eigenvalues;
        std::sort(computed.begin(), computed.end());
        std::sort(expected.begin(), expected.end());

        out << "\n  Computed eigenvalues (sorted):\n";
        for (std::size_t i = 0; i < 4; ++i)
            out << "    λ[" << i << "] = " << std::setw(18) << computed[i]
                << "  err = " << std::scientific << std::abs(computed[i] - expected[i]) << "\n";
        out << "\n";

        out << "  Orthogonality error:  " << std::scientific << orthoErr << "\n";
        out << "  Similarity error:     " << std::scientific << simErr << "\n";
        out << "\n";
    }

    // ==================================================================
    // Test 4: 8×8 Hilbert-like matrix
    // ==================================================================
    {
        out << "--- Test 4: 8×8 Symmetric Cauchy Matrix A_ij = 1/(1+i+j) ---\n\n";

        std::size_t n = 8;
        numalg::Matrix A(n, n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                A(i, j) = 1.0 / (1.0 + i + j);

        auto result = numalg::realSchurDecomposition(A);
        double orthoErr = orthogonalityError(result.Q);
        double simErr = schurSimilarityError(A, result.Q, result.T);

        out << "  Converged:  " << (result.converged ? "Yes" : "No")
            << "  (iterations: " << result.iterations << ")\n\n";

        std::vector<double> evals = result.eigenvalues;
        std::sort(evals.begin(), evals.end(), std::greater<double>());

        out << "  Computed eigenvalues (sorted descending):\n";
        for (std::size_t i = 0; i < evals.size(); ++i)
            out << "    λ[" << i << "] = " << std::scientific << std::setw(18) << evals[i] << "\n";
        out << "\n";

        out << "  Orthogonality error:  " << std::scientific << orthoErr << "\n";
        out << "  Similarity error:     " << std::scientific << simErr << "\n";
        out << "\n";
    }

    // ==================================================================
    // Test 5: 5×5 symmetric matrix
    // ==================================================================
    {
        out << "--- Test 5: 5×5 Symmetric Random Matrix ---\n\n";

        numalg::Matrix A(5, 5, {
            10.0,  2.0,  3.0,  1.0,  4.0,
             2.0,  8.0,  1.0,  5.0,  2.0,
             3.0,  1.0, 12.0,  2.0,  3.0,
             1.0,  5.0,  2.0,  6.0,  1.0,
             4.0,  2.0,  3.0,  1.0,  9.0
        });

        auto result = numalg::realSchurDecomposition(A);
        double orthoErr = orthogonalityError(result.Q);
        double simErr = schurSimilarityError(A, result.Q, result.T);

        out << "  Converged:  " << (result.converged ? "Yes" : "No")
            << "  (iterations: " << result.iterations << ")\n\n";

        std::vector<double> evals = result.eigenvalues;
        std::sort(evals.begin(), evals.end(), std::greater<double>());

        double trace = 0;
        for (std::size_t i = 0; i < 5; ++i) trace += evals[i];
        double original_trace = 10 + 8 + 12 + 6 + 9;

        out << "  Computed eigenvalues (sorted descending):\n";
        for (std::size_t i = 0; i < evals.size(); ++i)
            out << "    λ[" << i << "] = " << std::setw(18) << evals[i] << "\n";
        out << "\n";

        out << "  Trace(Λ) = " << trace << "  (original trace = " << original_trace << ")\n";
        out << "  Orthogonality error:  " << std::scientific << orthoErr << "\n";
        out << "  Similarity error:     " << std::scientific << simErr << "\n";
        out << "\n";
    }

    // ==================================================================
    // Summary
    // ==================================================================
    out << "============================================================\n";
    out << "Results Summary\n";
    out << "============================================================\n\n";

    out << "  All tests completed.\n";
    out << "\n";

    out.close();
}

int main() {
    std::filesystem::create_directories("./output");
    runTests();
    std::cout << "Results written to ./output/homework14.txt" << std::endl;
    return 0;
}
