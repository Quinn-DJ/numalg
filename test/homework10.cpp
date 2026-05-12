#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>

#include "matrix.hpp"
#include "operations.hpp"
#include "solve.hpp"
#include "gauss.hpp"
#include "iterative.hpp"

// Infinity-norm of a vector
double vecInfNorm(const numalg::Matrix& v) {
    double d = 0.0;
    for (std::size_t i = 0; i < v.lines(); ++i)
        d = std::max(d, std::abs(v(i)));
    return d;
}

// Print iteration result summary
void printResult(std::ofstream& out, const std::string& method,
                 const numalg::IterResult& r, const numalg::Matrix& x_exact) {
    out << "  " << method << ":\n";
    out << "    Converged:  " << (r.converged ? "Yes" : "No") << "\n";
    out << "    Iterations: " << r.iterations << "\n";
    if (x_exact.lines() > 0) {
        double err = vecInfNorm(r.x - x_exact);
        out << "    ||x - x*||_inf: " << std::scientific << err << "\n";
    }
    out << "\n";
}

// Print residual history (first few + last)
void printHistory(std::ofstream& out, const std::string& method,
                  const std::vector<double>& hist) {
    out << "  " << method << " residual history:\n    ";
    std::size_t show = std::min(hist.size(), std::size_t(10));
    for (std::size_t i = 0; i < show; ++i)
        out << std::scientific << hist[i] << " ";
    if (hist.size() > 10)
        out << "... " << hist.back();
    out << "\n\n";
}

// ============================================
// Problem 1: 3x3 strictly diagonally dominant system
//   [10  -1  -2] [x1]   [  7.2]
//   [-1  11  -1] [x2] = [  7.3]
//   [-2  -1  10] [x3]   [  8.4]
// ============================================
void prob1(std::ofstream& out) {
    out << "========== Problem 1: 3x3 diagonally dominant system ==========\n";
    out << "  A = [[10, -1, -2], [-1, 11, -1], [-2, -1, 10]]\n";
    out << "  b = [7.2, 7.3, 8.4]^T\n";

    numalg::Matrix A(3, 3, {10, -1, -2, -1, 11, -1, -2, -1, 10});
    numalg::Matrix b(3, {7.2, 7.3, 8.4});
    numalg::Matrix x_exact = PgaussSolve(A, b);

    out << std::fixed << std::setprecision(6);
    out << "  Exact: x = [" << x_exact(0) << ", " << x_exact(1) << ", " << x_exact(2) << "]^T\n\n";

    auto r_jacobi = numalg::jacobiSolve(A, b);
    auto r_gs     = numalg::gaussSeidelSolve(A, b);
    auto r_sor    = numalg::sorSolve(A, b, 1.1);

    printResult(out, "Jacobi", r_jacobi, x_exact);
    printResult(out, "Gauss-Seidel", r_gs, x_exact);
    printResult(out, "SOR (omega=1.1)", r_sor, x_exact);

    // Compare convergence speed
    out << "  Convergence speed (iterations to tol=1e-10):\n";
    out << "    Jacobi:       " << r_jacobi.iterations << "\n";
    out << "    Gauss-Seidel: " << r_gs.iterations << "\n";
    out << "    SOR(1.1):     " << r_sor.iterations << "\n";
    out << "    SOR(1.2):     " << numalg::sorSolve(A, b, 1.2).iterations << "\n";
    out << "    SOR(1.5):     " << numalg::sorSolve(A, b, 1.5).iterations << "\n\n";

    printHistory(out, "Jacobi", r_jacobi.residual_history);
    printHistory(out, "Gauss-Seidel", r_gs.residual_history);
}

// ============================================
// Problem 2: 5x5 strictly diagonally dominant system
//   A is tridiagonal with 4 on diagonal, -1 on sub/super diagonals
//   b = ones
// ============================================
void prob2(std::ofstream& out) {
    out << "========== Problem 2: 5x5 tridiagonal system ==========\n";
    out << "  A: diagonal=4, sub/super=-1\n";
    out << "  b = [1, 1, 1, 1, 1]^T\n\n";

    int n = 5;
    numalg::Matrix A(n, n);
    for (int i = 0; i < n; ++i) A(i, i) = 4;
    for (int i = 0; i < n - 1; ++i) {
        A(i, i + 1) = -1;
        A(i + 1, i) = -1;
    }
    numalg::Matrix b(n, {1, 1, 1, 1, 1});

    // Exact solution via Gaussian elimination
    numalg::Matrix x_exact = PgaussSolve(A, b);

    auto r_jacobi = numalg::jacobiSolve(A, b);
    auto r_gs     = numalg::gaussSeidelSolve(A, b);
    auto r_sor    = numalg::sorSolve(A, b, 1.3);

    printResult(out, "Jacobi", r_jacobi, x_exact);
    printResult(out, "Gauss-Seidel", r_gs, x_exact);
    printResult(out, "SOR (omega=1.3)", r_sor, x_exact);

    out << "  Convergence speed:\n";
    out << "    Jacobi:       " << r_jacobi.iterations << "\n";
    out << "    Gauss-Seidel: " << r_gs.iterations << "\n";
    out << "    SOR(1.0):     " << numalg::sorSolve(A, b, 1.0).iterations << "\n";
    out << "    SOR(1.2):     " << numalg::sorSolve(A, b, 1.2).iterations << "\n";
    out << "    SOR(1.3):     " << r_sor.iterations << "\n";
    out << "    SOR(1.5):     " << numalg::sorSolve(A, b, 1.5).iterations << "\n";
    out << "    SOR(1.7):     " << numalg::sorSolve(A, b, 1.7).iterations << "\n\n";

    printHistory(out, "Jacobi", r_jacobi.residual_history);
    printHistory(out, "Gauss-Seidel", r_gs.residual_history);
    printHistory(out, "SOR(1.3)", r_sor.residual_history);
}

// ============================================
// Problem 3: Large 100x100 diagonally dominant tridiagonal system
//   Compare convergence speed of three methods at different scales
// ============================================
void prob3(std::ofstream& out) {
    out << "========== Problem 3: 100x100 tridiagonal system ==========\n";
    out << "  A: diagonal=10, sub/super=1 (strongly diagonally dominant)\n";
    out << "  b = [1, 2, 3, ..., 100]^T\n\n";

    int n = 100;
    numalg::Matrix A(n, n);
    for (int i = 0; i < n; ++i) A(i, i) = 10;
    for (int i = 0; i < n - 1; ++i) {
        A(i, i + 1) = 1;
        A(i + 1, i) = 1;
    }
    numalg::Matrix b(n);
    for (int i = 0; i < n; ++i) b(i) = i + 1;

    auto r_jacobi = numalg::jacobiSolve(A, b, 1e-12);
    auto r_gs     = numalg::gaussSeidelSolve(A, b, 1e-12);
    auto r_sor    = numalg::sorSolve(A, b, 1.2, 1e-12);

    out << std::scientific << std::setprecision(6);
    printResult(out, "Jacobi", r_jacobi, numalg::Matrix(0, 0));
    printResult(out, "Gauss-Seidel", r_gs, numalg::Matrix(0, 0));
    printResult(out, "SOR (omega=1.2)", r_sor, numalg::Matrix(0, 0));

    out << "  Convergence speed (tol=1e-12):\n";
    out << "    Jacobi:       " << r_jacobi.iterations << "\n";
    out << "    Gauss-Seidel: " << r_gs.iterations << "\n";
    out << "    SOR (1.1):     " << numalg::sorSolve(A, b, 1.1, 1e-12).iterations << "\n";
    out << "    SOR (1.2):     " << r_sor.iterations << "\n";
    out << "    SOR (1.5):     " << numalg::sorSolve(A, b, 1.5, 1e-12).iterations << "\n\n";

    printHistory(out, "Jacobi", r_jacobi.residual_history);
    printHistory(out, "Gauss-Seidel", r_gs.residual_history);
    printHistory(out, "SOR(1.2)", r_sor.residual_history);
}

// ============================================
// Problem 4: Non-diagonally dominant system
//   Shows that iterative methods may diverge
//   [1  2] [x1]   [3]
//   [3  4] [x2] = [7]
// ============================================
void prob4(std::ofstream& out) {
    out << "========== Problem 4: Non-diagonally dominant system ==========\n";
    out << "  A = [[1, 2], [3, 4]] (not diagonally dominant)\n";
    out << "  b = [3, 7]^T\n\n";

    numalg::Matrix A(2, 2, {1, 2, 3, 4});
    numalg::Matrix b(2, {3, 7});

    auto r_jacobi = numalg::jacobiSolve(A, b, 1e-10, 1000);
    auto r_gs     = numalg::gaussSeidelSolve(A, b, 1e-10, 1000);

    out << "  Jacobi:\n";
    out << "    Converged: " << (r_jacobi.converged ? "Yes" : "No") << "\n";
    out << "    Iterations: " << r_jacobi.iterations << "\n";
    if (r_jacobi.converged)
        out << "    x = [" << r_jacobi.x(0) << ", " << r_jacobi.x(1) << "]\n";

    out << "  Gauss-Seidel:\n";
    out << "    Converged: " << (r_gs.converged ? "Yes" : "No") << "\n";
    out << "    Iterations: " << r_gs.iterations << "\n";
    if (r_gs.converged)
        out << "    x = [" << r_gs.x(0) << ", " << r_gs.x(1) << "]\n";
    out << "\n";
}

int main() {
    std::filesystem::create_directories("./output");
    std::ofstream out("./output/homework10.txt");
    if (!out.is_open()) {
        std::cerr << "Failed to open output file: ./output/homework10.txt" << std::endl;
        return 1;
    }

    prob1(out);
    prob2(out);
    prob3(out);
    prob4(out);

    std::cout << "Results written to ./output/homework10.txt" << std::endl;
    return 0;
}
