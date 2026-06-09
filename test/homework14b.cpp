#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "sym_eigen.hpp"
#include "operations.hpp"

// Check if a matrix is tridiagonal (within tolerance)
bool isTridiagonal(const numalg::Matrix& T, double tol = 1e-12) {
    std::size_t n = T.lines();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if ((i > j + 1 || j > i + 1) && std::abs(T(i, j)) > tol) return false;
    return true;
}

bool isDiagonal(const numalg::Matrix& T, double tol = 1e-12) {
    std::size_t n = T.lines();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (i != j && std::abs(T(i, j)) > tol) return false;
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
double similarityError(const numalg::Matrix& A, const numalg::Matrix& Q, const numalg::Matrix& T) {
    numalg::Matrix QtAQ = numalg::tr(Q) * A * Q;
    double maxErr = 0.0;
    for (std::size_t i = 0; i < A.lines(); ++i)
        for (std::size_t j = 0; j < A.lines(); ++j)
            maxErr = std::max(maxErr, std::abs(QtAQ(i, j) - T(i, j)));
    return maxErr;
}

void runTests() {
    std::ofstream out("./output/homework14b.txt");
    if (!out.is_open()) {
        std::cerr << "Failed to open output file.\n";
        return;
    }

    out << std::scientific << std::setprecision(15);
    out << "============================================================\n";
    out << "Homework 14b: Implicit Symmetric QR Algorithm\n";
    out << "============================================================\n\n";

    // ==================================================================
    // Test 1: Tridiagonal reduction (Algorithm 7.2.1)
    // ==================================================================
    {
        out << "--- Test 1: Tridiagonal Reduction (Algorithm 7.2.1) ---\n\n";

        numalg::Matrix A(5, 5, {
            4.0,  1.0,  2.0,  3.0,  1.0,
            1.0,  3.0,  1.0,  2.0,  2.0,
            2.0,  1.0,  5.0,  1.0,  3.0,
            3.0,  2.0,  1.0,  2.0,  1.0,
            1.0,  2.0,  3.0,  1.0,  4.0
        });

        numalg::Matrix T, Q;
        numalg::tridiagonalReduction(A, T, Q);

        out << "  Tridiagonal form T:\n";
        for (std::size_t i = 0; i < 5; ++i) {
            out << "    ";
            for (std::size_t j = 0; j < 5; ++j)
                out << std::setw(16) << T(i, j);
            out << "\n";
        }
        out << "\n";

        double orthoErr = orthogonalityError(Q);
        double simErr = similarityError(A, Q, T);

        out << "  Orthogonality error ||Q^T Q - I||:   " << orthoErr << "\n";
        out << "  Similarity error ||Q^T A Q - T||_inf: " << simErr << "\n";
        out << "  Is tridiagonal: " << (isTridiagonal(T) ? "Yes" : "No") << "\n";
        out << "\n";
    }

    // ==================================================================
    // Test 2: 3×3 symmetric matrix (full spectral decomposition)
    // ==================================================================
    {
        out << "--- Test 2: 3×3 Symmetric ─ Full Spectral Decomposition ---\n\n";

        numalg::Matrix A(3, 3, {
            3.0, 1.0, 1.0,
            1.0, 3.0, 1.0,
            1.0, 1.0, 3.0
        });

        auto result = numalg::symmetricEigensolver(A);

        std::vector<double> expected = {5.0, 2.0, 2.0};
        std::vector<double> computed = result.eigenvalues;
        std::sort(computed.begin(), computed.end(), std::greater<double>());

        out << "  Converged: " << (result.converged ? "Yes" : "No")
            << "  (iterations: " << result.iterations << ")\n\n";

        out << "  Eigenvalues:\n";
        for (std::size_t i = 0; i < 3; ++i)
            out << "    λ[" << i << "] = " << std::setw(18) << computed[i]
                << "  (expected " << expected[i] << ")"
                << "  err = " << std::abs(computed[i] - expected[i]) << "\n";
        out << "\n";

        double orthoErr = orthogonalityError(result.Q);
        double simErr = similarityError(A, result.Q, result.T);
        bool diag = isDiagonal(result.T);

        out << "  Orthogonality error:  " << orthoErr << "\n";
        out << "  Similarity error:     " << simErr << "\n";
        out << "  Is diagonal:          " << (diag ? "Yes" : "No") << "\n";
        out << "\n";
    }

    // ==================================================================
    // Test 3: 4×4 Poisson matrix
    // ==================================================================
    {
        out << "--- Test 3: 4×4 Poisson ─ Full Spectral Decomposition ---\n\n";

        numalg::Matrix A(4, 4, {
            2.0, -1.0,  0.0,  0.0,
            -1.0,  2.0, -1.0,  0.0,
            0.0, -1.0,  2.0, -1.0,
            0.0,  0.0, -1.0,  2.0
        });

        auto result = numalg::symmetricEigensolver(A);

        std::vector<double> expected(4);
        for (int k = 1; k <= 4; ++k)
            expected[k-1] = 2.0 - 2.0 * std::cos(k * M_PI / 5.0);

        std::vector<double> computed = result.eigenvalues;
        std::sort(computed.begin(), computed.end());
        std::sort(expected.begin(), expected.end());

        out << "  Converged: " << (result.converged ? "Yes" : "No")
            << "  (iterations: " << result.iterations << ")\n\n";

        out << "  Expected eigenvalues:\n";
        for (std::size_t i = 0; i < 4; ++i)
            out << "    λ[" << i << "] = " << std::setw(18) << expected[i] << "\n";

        out << "\n  Computed eigenvalues:\n";
        double maxErr = 0.0;
        for (std::size_t i = 0; i < 4; ++i) {
            double err = std::abs(computed[i] - expected[i]);
            maxErr = std::max(maxErr, err);
            out << "    λ[" << i << "] = " << std::setw(18) << computed[i]
                << "  err = " << err << "\n";
        }
        out << "\n";

        double orthoErr = orthogonalityError(result.Q);
        double simErr = similarityError(A, result.Q, result.T);

        out << "  Max eigenvalue error:  " << maxErr << "\n";
        out << "  Orthogonality error:   " << orthoErr << "\n";
        out << "  Similarity error:      " << simErr << "\n";
        out << "  Is diagonal:           " << (isDiagonal(result.T) ? "Yes" : "No") << "\n";
        out << "\n";
    }

    // ==================================================================
    // Test 4: 8×8 Cauchy matrix
    // ==================================================================
    {
        out << "--- Test 4: 8×8 Symmetric Cauchy Matrix ---\n\n";

        std::size_t n = 8;
        numalg::Matrix A(n, n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                A(i, j) = 1.0 / (1.0 + i + j);

        auto result = numalg::symmetricEigensolver(A);

        std::vector<double> evals = result.eigenvalues;
        std::sort(evals.begin(), evals.end(), std::greater<double>());

        double trace = 0;
        for (auto v : evals) trace += v;

        out << "  Converged: " << (result.converged ? "Yes" : "No")
            << "  (iterations: " << result.iterations << ")\n\n";

        out << "  Eigenvalues (sorted descending):\n";
        for (std::size_t i = 0; i < evals.size(); ++i)
            out << "    λ[" << i << "] = " << evals[i] << "\n";
        out << "\n";

        double orthoErr = orthogonalityError(result.Q);
        double simErr = similarityError(A, result.Q, result.T);

        out << "  Sum(Λ) = " << trace << "  (trace(A) = " << n << ")\n";
        out << "  Orthogonality error:  " << orthoErr << "\n";
        out << "  Similarity error:     " << simErr << "\n";
        out << "  Is diagonal:          " << (isDiagonal(result.T) ? "Yes" : "No") << "\n";
        out << "\n";
    }

    // ==================================================================
    // Test 5: 5×5 random symmetric matrix
    // ==================================================================
    {
        out << "--- Test 5: 5×5 Random Symmetric Matrix ---\n\n";

        numalg::Matrix A(5, 5, {
            10.0,  2.0,  3.0,  1.0,  4.0,
             2.0,  8.0,  1.0,  5.0,  2.0,
             3.0,  1.0, 12.0,  2.0,  3.0,
             1.0,  5.0,  2.0,  6.0,  1.0,
             4.0,  2.0,  3.0,  1.0,  9.0
        });

        auto result = numalg::symmetricEigensolver(A);

        std::vector<double> evals = result.eigenvalues;
        std::sort(evals.begin(), evals.end(), std::greater<double>());

        double traceA = 10.0 + 8.0 + 12.0 + 6.0 + 9.0;
        double traceE = 0;
        for (auto v : evals) traceE += v;

        out << "  Converged: " << (result.converged ? "Yes" : "No")
            << "  (iterations: " << result.iterations << ")\n\n";

        out << "  Eigenvalues (sorted descending):\n";
        for (std::size_t i = 0; i < evals.size(); ++i)
            out << "    λ[" << i << "] = " << std::setw(18) << evals[i] << "\n";
        out << "\n";

        double orthoErr = orthogonalityError(result.Q);
        double simErr = similarityError(A, result.Q, result.T);

        out << "  Trace(A) = " << traceA << "  trace(Λ) = " << traceE << "\n";
        out << "  Orthogonality error:  " << orthoErr << "\n";
        out << "  Similarity error:     " << simErr << "\n";
        out << "  Is diagonal:          " << (isDiagonal(result.T) ? "Yes" : "No") << "\n";
        out << "\n";
    }

    // ==================================================================
    // Summary
    // ==================================================================
    out << "============================================================\n";
    out << "Results Summary\n";
    out << "============================================================\n\n";
    out << "  All tests completed.\n\n";

    out.close();
}

int main() {
    std::filesystem::create_directories("./output");
    runTests();
    std::cout << "Results written to ./output/homework14b.txt" << std::endl;
    return 0;
}
