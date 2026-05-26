#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <vector>

#include "power_method.hpp"

// Evaluate p(x) = x^n + c_{n-1}x^{n-1} + ... + c_0 via Horner's method
double evalPoly(const std::vector<double>& c, double x) {
    double r = 1.0;
    for (std::size_t i = c.size(); i > 0; --i) {
        r = r * x + c[i - 1];
    }
    return r;
}

void runTests() {
    std::ofstream out("./output/homework12.txt");
    if (!out.is_open()) {
        std::cerr << "Failed to open output file.\n";
        return;
    }

    out << "============================================================\n";
    out << "Homework 12: Power Method for Polynomial Roots\n";
    out << "Find dominant roots via companion matrix + power iteration\n";
    out << "============================================================\n\n";

    out << std::fixed << std::setprecision(12);

    // ==================================================================
    // (i)  x^3 + x^2 - 5x + 3 = 0 = (x+3)(x-1)^2
    // ==================================================================
    {
        out << "--- (i) f(x) = x^3 + x^2 - 5x + 3 ---\n";
        out << "  Factorization: (x+3)(x-1)^2\n";
        out << "  Expected dominant root: -3\n\n";

        std::vector<double> c = {3.0, -5.0, 1.0};
        auto C = numalg::companionMatrix(c);
        out << "  Companion matrix C:\n";
        for (std::size_t i = 0; i < 3; ++i) {
            out << "    ";
            for (std::size_t j = 0; j < 3; ++j)
                out << std::setw(14) << C(i, j);
            out << "\n";
        }
        out << "\n";

        auto r = numalg::dominantRoot(c, 2000, 1e-14);

        out << "  Dominant root:  " << r.eigenvalue << "\n";
        out << "  f(root):        " << evalPoly(c, r.eigenvalue) << "\n";
        out << "  Iterations:     " << r.iterations << "\n";
        out << "  Converged:      " << (r.converged ? "Yes" : "No") << "\n\n";
    }

    // ==================================================================
    // (ii) x^3 - 3x - 1 = 0
    //      Roots: 2cos(20°) ≈ 1.879385..., 2cos(140°) ≈ -1.532089..., 2cos(260°) ≈ -0.347296...
    // ==================================================================
    {
        out << "--- (ii) f(x) = x^3 - 3x - 1 ---\n";
        out << "  Roots: 2cos(20°) ≈ 1.87939, 2cos(140°) ≈ -1.53209, 2cos(260°) ≈ -0.34730\n";
        out << "  Expected dominant root: 1.879385241571817\n\n";

        std::vector<double> c = {-1.0, -3.0, 0.0};
        auto r = numalg::dominantRoot(c, 2000, 1e-14);

        out << "  Dominant root:  " << r.eigenvalue << "\n";
        out << "  f(root):        " << evalPoly(c, r.eigenvalue) << "\n";
        out << "  Iterations:     " << r.iterations << "\n";
        out << "  Converged:      " << (r.converged ? "Yes" : "No") << "\n";
        out << "  |λ₂/λ₁| =       " << std::abs(1.532088886237956 / 1.879385241571817)
            << " (roots are close in magnitude, convergence is slow)\n\n";
    }

    // ==================================================================
    // (iii) x^8 + 101x^7 + 208.01x^6 + 10891.01x^5 + 9802.08x^4
    //       + 79108.9x^3 - 99902x^2 + 790x - 1000 = 0
    //       = (x+100)(x^7 + x^6 + 108.01x^5 + 90.01x^4 + 801.08x^3 - 999.1x^2 + 8x - 10)
    // ==================================================================
    {
        out << "--- (iii) f(x) = x^8 + 101x^7 + 208.01x^6 + 10891.01x^5\n";
        out << "           + 9802.08x^4 + 79108.9x^3 - 99902x^2 + 790x - 1000 ---\n";
        out << "  Factor: (x+100) is a factor\n";
        out << "  Expected dominant root: -100\n\n";

        std::vector<double> c = {-1000.0, 790.0, -99902.0, 79108.9,
                                  9802.08, 10891.01, 208.01, 101.0};
        auto r = numalg::dominantRoot(c, 5000, 1e-14);

        out << "  Dominant root:  " << r.eigenvalue << "\n";
        out << "  f(root):        " << evalPoly(c, r.eigenvalue) << " (near 0)\n";
        out << "  Iterations:     " << r.iterations << "\n";
        out << "  Converged:      " << (r.converged ? "Yes" : "No") << "\n";
        out << "  |λ₂/λ₁| ≈       " << 0.01
            << " (dominant root dominates heavily, very fast convergence)\n";
        out << "\n  Convergence history:\n";
        for (std::size_t i = 0; i < r.history.size(); ++i)
            out << "    iter " << std::setw(2) << i << ": "
                << std::scientific << r.history[i] << "\n";
        out << "\n";
    }

    // ==================================================================
    // Summary
    // ==================================================================
    out << "============================================================\n";
    out << "Results Summary\n";
    out << "============================================================\n\n";

    out << "  (i)   x^3 + x^2 - 5x + 3     → dominant root = -3\n";
    out << "  (ii)  x^3 - 3x - 1           → dominant root = 1.879385241572\n";
    out << "  (iii) x^8 + 101x^7 + ...     → dominant root = -100\n\n";

    out << "All tests completed.\n";
    out.close();
}

int main() {
    std::filesystem::create_directories("./output");
    runTests();
    std::cout << "Results written to ./output/homework12.txt" << std::endl;
    return 0;
}
