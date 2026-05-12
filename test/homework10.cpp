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

// Exact solution of the BVP: y(x) = (1-a)/(1-e^{-1/eps}) * (1 - e^{-x/eps}) + a*x
double exactSolution(double x, double eps, double a) {
    double coeff = (1.0 - a) / (1.0 - std::exp(-1.0 / eps));
    return coeff * (1.0 - std::exp(-x / eps)) + a * x;
}

// Build the (n-1)x(n-1) tridiagonal coefficient matrix A
numalg::Matrix buildMatrix(std::size_t n, double eps) {
    double h = 1.0 / n;
    std::size_t m = n - 1;  // interior points
    numalg::Matrix A(m, m);

    double diag = -(2.0 * eps + h);
    double upper = eps + h;
    double lower = eps;

    for (std::size_t i = 0; i < m; ++i) {
        A(i, i) = diag;
        if (i > 0) A(i, i - 1) = lower;
        if (i < m - 1) A(i, i + 1) = upper;
    }
    return A;
}

// Build the RHS vector b (n-1 elements)
numalg::Matrix buildRHS(std::size_t n, double eps, double a) {
    double h = 1.0 / n;
    std::size_t m = n - 1;
    numalg::Matrix b(m);

    for (std::size_t i = 0; i < m; ++i) {
        b(i) = a * h * h;
    }
    // Boundary conditions
    b(0) -= eps * 0.0;        // y(0) = 0, so no contribution
    b(m - 1) -= (eps + h) * 1.0;  // y(1) = 1

    return b;
}

// Compute max error ||x_iter - x_exact||_inf
double maxError(const numalg::Matrix& x_iter,
                std::size_t n, double eps, double a) {
    double err = 0.0;
    for (std::size_t i = 0; i < n - 1; ++i) {
        double x = (i + 1) / static_cast<double>(n);
        double exact = exactSolution(x, eps, a);
        err = std::max(err, std::abs(x_iter(i) - exact));
    }
    return err;
}

// Solve and report for one (eps, a, n) configuration
void solveAndReport(std::ofstream& out, double eps, double a, std::size_t n) {
    out << "eps = " << std::scientific << std::setprecision(4) << eps << ", "
        << "a = " << std::fixed << std::setprecision(2) << a << ", "
        << "n = " << n << "\n\n";

    numalg::Matrix A = buildMatrix(n, eps);
    numalg::Matrix b = buildRHS(n, eps, a);

    // Jacobi
    auto r_jacobi = numalg::jacobiSolve(A, b, 1e-10, 100000);
    // Gauss-Seidel
    auto r_gs = numalg::gaussSeidelSolve(A, b, 1e-10, 100000);
    // SOR with a few omega values
    auto r_sor1 = numalg::sorSolve(A, b, 1.2, 1e-10, 100000);
    auto r_sor2 = numalg::sorSolve(A, b, 1.5, 1e-10, 100000);
    auto r_sor3 = numalg::sorSolve(A, b, 1.8, 1e-10, 100000);

    double err_jacobi = maxError(r_jacobi.x, n, eps, a);
    double err_gs     = maxError(r_gs.x, n, eps, a);
    double err_sor1   = maxError(r_sor1.x, n, eps, a);
    double err_sor2   = maxError(r_sor2.x, n, eps, a);
    double err_sor3   = maxError(r_sor3.x, n, eps, a);

    out << std::scientific << std::setprecision(6);
    out << "  Method        Converged  Iterations   ||x - x*||_inf\n";
    out << "  ----------    --------   ----------   ---------------\n";
    out << "  Jacobi        " << (r_jacobi.converged ? "Yes" : "No ")
        << "       " << std::setw(10) << r_jacobi.iterations
        << "   " << err_jacobi << "\n";
    out << "  Gauss-Seidel  " << (r_gs.converged ? "Yes" : "No ")
        << "       " << std::setw(10) << r_gs.iterations
        << "   " << err_gs << "\n";
    out << "  SOR(1.2)      " << (r_sor1.converged ? "Yes" : "No ")
        << "       " << std::setw(10) << r_sor1.iterations
        << "   " << err_sor1 << "\n";
    out << "  SOR(1.5)      " << (r_sor2.converged ? "Yes" : "No ")
        << "       " << std::setw(10) << r_sor2.iterations
        << "   " << err_sor2 << "\n";
    out << "  SOR(1.8)      " << (r_sor3.converged ? "Yes" : "No ")
        << "       " << std::setw(10) << r_sor3.iterations
        << "   " << err_sor3 << "\n\n";

    // Print first few and last few solution values
    out << std::fixed << std::setprecision(6);
    out << "  Sample solutions (interior points):\n";
    out << "  i       x_i         Jacobi      G-S         SOR(1.5)    Exact\n";
    out << "  -----   --------    --------    --------    --------    --------\n";

    std::vector<int> indices;
    if (n <= 20) {
        for (std::size_t i = 0; i < n - 1; ++i) indices.push_back(i);
    } else {
        for (int i = 0; i < 5; ++i) indices.push_back(i);
        for (int i = static_cast<int>(n) - 7; i < static_cast<int>(n) - 1; ++i)
            indices.push_back(i);
    }

    for (int idx : indices) {
        double xi = (idx + 1) / static_cast<double>(n);
        double ex = exactSolution(xi, eps, a);
        out << "  " << std::setw(5) << idx << "   "
            << std::setw(10) << xi << "  "
            << std::setw(10) << r_jacobi.x(idx) << "  "
            << std::setw(10) << r_gs.x(idx) << "  "
            << std::setw(10) << r_sor2.x(idx) << "  "
            << std::setw(10) << ex << "\n";
    }
    out << "\n" << std::string(60, '-') << "\n\n";
}

int main() {
    std::filesystem::create_directories("./output");
    std::ofstream out("./output/homework10.txt");
    if (!out.is_open()) {
        std::cerr << "Failed to open output file: ./output/homework10.txt" << std::endl;
        return 1;
    }

    out << "Week 10: Iterative methods for BVP\n";
    out << "BVP: eps*y'' + y' = a, y(0)=0, y(1)=1\n";
    out << std::string(60, '-') << "\n\n";

    double a = 0.5;
    std::size_t n = 100;

    double eps_values[] = {1.0, 0.1, 0.01, 0.0001};

    for (double eps : eps_values) {
        solveAndReport(out, eps, a, n);
    }

    std::cout << "Results written to ./output/homework10.txt" << std::endl;
    return 0;
}
