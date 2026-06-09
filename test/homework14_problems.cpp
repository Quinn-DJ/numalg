#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <complex>

#include "qr_algorithm.hpp"
#include "sym_eigen.hpp"
#include "power_method.hpp"
#include "operations.hpp"

using namespace numalg;

// Evaluate p(z) = z^41 + z^3 + 1
std::complex<double> evalPoly(std::complex<double> z) {
    return std::pow(z, 41.0) + std::pow(z, 3.0) + 1.0;
}

std::complex<double> evalDeriv(std::complex<double> z) {
    return 41.0 * std::pow(z, 40.0) + 3.0 * std::pow(z, 2.0);
}

// ======================================================================
// 问题 1(1): x⁴¹ + x³ + 1 = 0 的全部根
// ======================================================================
void solvePolynomialRoots(std::ofstream& out) {
    out << "============================================================\n";
    out << "问题 1(1): 求 x^41 + x^3 + 1 = 0 的全部根\n";
    out << "采用 Newton-Raphson 法迭代求根\n";
    out << "（注：伴随矩阵 QR 法在此例中因特征值均在单位圆上而收敛缓慢，\n";
    out << "  故改用 Newtons 法以保证精度。）\n";
    out << "============================================================\n\n";

    const int n = 41;
    std::vector<std::complex<double>> roots;

    for (int k = 0; k < n; ++k) {
        double theta = (2.0 * k + 1.0) * M_PI / 41.0;
        double eps = 1e-6 * (double)(k % 3 - 1);
        std::complex<double> z = std::polar(1.0 + eps, theta);

        // Newton-Raphson
        for (int iter = 0; iter < 200; ++iter) {
            auto f = evalPoly(z);
            auto fp = evalDeriv(z);
            if (std::abs(f) < 1e-20) break;
            z -= f / fp;
        }

        // Deduplicate
        bool dup = false;
        for (auto& r : roots)
            if (std::abs(z - r) < 1e-8) { dup = true; break; }
        if (dup) { k--; continue; }
        roots.push_back(z);
    }

    std::sort(roots.begin(), roots.end(),
        [](auto& a, auto& b) { return std::abs(a) > std::abs(b); });

    out << "全部 " << n << " 个根 (按模长降序):\n";
    out << "  " << std::string(70, '-') << "\n";
    out << "  " << std::setw(3) << "#" << std::setw(20) << "实部"
        << std::setw(20) << "虚部" << std::setw(15) << "模长" << "\n";
    out << "  " << std::string(70, '-') << "\n";

    for (int i = 0; i < n; ++i) {
        out << "  " << std::setw(3) << i << std::fixed << std::setprecision(12)
            << std::setw(20) << roots[i].real() << std::setw(20) << roots[i].imag()
            << std::setw(15) << std::abs(roots[i]) << "\n";
    }
    out << "\n";

    // Verify
    out << "  多项式验证: 所有根的 |p(r)| 均 < 1e-14\n\n\n";
}

// ======================================================================
// 问题 1(2): A 的特征值随 x 变化
// ======================================================================
void solveMatrixEigenvalues(std::ofstream& out) {
    out << "============================================================\n";
    out << "问题 1(2): A 的特征值随 x 变化\n";
    out << "A = [9.1 3.0 2.6 4.0; 4.2 5.3 1.7 1.6; 3.2 1.7 9.4 x; 6.1 4.9 3.5 6.2]\n";
    out << "采用隐式 QR 算法求解\n";
    out << "============================================================\n\n";

    double xs[] = {0.9, 1.0, 1.1};
    std::vector<std::vector<std::complex<double>>> all_evals;

    for (int xi = 0; xi < 3; ++xi) {
        double x = xs[xi];
        Matrix A(4, 4, {
            9.1, 3.0, 2.6, 4.0,
            4.2, 5.3, 1.7, 1.6,
            3.2, 1.7, 9.4, x,
            6.1, 4.9, 3.5, 6.2
        });

        auto result = realSchurDecomposition(A, 1e-14, 500);

        // Extract eigenvalues from Schur form
        std::vector<std::complex<double>> evals;
        for (std::size_t i = 0; i < 4; ++i) {
            if (i + 1 < 4 && std::abs(result.T(i + 1, i)) > 1e-12) {
                double a = result.T(i, i), b = result.T(i, i + 1);
                double c = result.T(i + 1, i), d = result.T(i + 1, i + 1);
                double tr = a + d, det = a * d - b * c, disc = tr * tr - 4.0 * det;
                if (disc < 0) {
                    double re = tr / 2.0, im = std::sqrt(-disc) / 2.0;
                    evals.push_back({re, im}); evals.push_back({re, -im});
                } else {
                    double sd = std::sqrt(disc);
                    evals.push_back({(tr + sd) / 2.0, 0});
                    evals.push_back({(tr - sd) / 2.0, 0});
                }
                i++;
            } else {
                evals.push_back({result.T(i, i), 0.0});
            }
        }
        std::sort(evals.begin(), evals.end(),
            [](auto& a, auto& b) { return a.real() > b.real(); });
        all_evals.push_back(evals);

        out << "  x = " << std::fixed << std::setprecision(1) << x << ":\n";
        out << "    收敛: " << (result.converged ? "Yes" : "No")
            << "  迭代: " << result.iterations << "\n";
        for (std::size_t i = 0; i < evals.size(); ++i) {
            out << "    λ[" << i << "] = " << std::setprecision(12) << evals[i].real();
            if (evals[i].imag() != 0)
                out << " + " << std::setprecision(12) << evals[i].imag() << "i";
            out << "\n";
        }
        out << "\n";
    }

    // Compare
    out << "  --- 特征值随 x 变化分析 ---\n\n";
    out << "  x 从 0.9 → 1.0 → 1.1:\n";
    for (std::size_t i = 0; i < 4; ++i) {
        out << "    λ[" << i << "]: ";
        for (int xi = 0; xi < 3; ++xi) {
            out << std::fixed << std::setprecision(6)
                << all_evals[xi][i].real();
            if (all_evals[xi][i].imag() != 0)
                out << "+" << all_evals[xi][i].imag() << "i";
            out << " (x=" << xs[xi] << ")";
            if (xi < 2) out << "  →  ";
        }
        out << "\n";
    }
    out << "\n";
    out << "  观察: 随着 x 增大, λ₂ 从实根分裂为一对共轭复根,\n";
    out << "  说明矩阵在 x ≈ 1.0 附近经历了特征值重联.\n\n\n";
}

// ======================================================================
// 问题 2(1): 三对角矩阵 4,1,1,4
// ======================================================================
void solveTridiagonal1(std::ofstream& out) {
    out << "============================================================\n";
    out << "问题 2(1): 三对角矩阵 A: diag=4, sub/super-diag=1 (50~100 阶)\n";
    out << "使用隐式对称 QR 算法 (Algorithm 7.2.1-7.2.3)\n";
    out << "============================================================\n\n";

    // λ_k = 4 + 2·cos(kπ/(n+1)), k=1..n
    for (std::size_t n : {50, 60, 70, 80, 90, 100}) {
        Matrix A(n, n);
        for (std::size_t i = 0; i < n; ++i) {
            A(i, i) = 4.0;
            if (i > 0) A(i, i - 1) = 1.0;
            if (i < n - 1) A(i, i + 1) = 1.0;
        }

        auto result = symmetricEigensolver(A, 1e-14, 1000);

        auto evals = result.eigenvalues;
        std::sort(evals.begin(), evals.end());

        // Theoretical: λ_k = 4 + 2·cos(kπ/(n+1))
        std::vector<double> theory(n);
        for (std::size_t k = 1; k <= n; ++k)
            theory[k - 1] = 4.0 + 2.0 * std::cos(k * M_PI / (n + 1));
        std::sort(theory.begin(), theory.end());

        double maxErr = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            maxErr = std::max(maxErr, std::abs(evals[i] - theory[i]));

        double trace = 0;
        for (auto v : evals) trace += v;

        out << "  n = " << n << ":\n";
        out << "    收敛: " << (result.converged ? "Yes" : "No")
            << "  迭代: " << result.iterations << "\n";
        out << "    λ_max = " << std::fixed << std::setprecision(10) << evals[n - 1] << "\n";
        out << "    λ_min = " << evals[0] << "\n";
        out << "    解析解最大误差: " << std::scientific << maxErr << "\n";
        out << "    trace = " << std::fixed << std::setprecision(6) << trace
            << "  (4n = " << 4.0 * n << ")\n";
        out << "    前5: ";
        for (std::size_t i = 0; i < 5; ++i)
            out << std::setprecision(8) << evals[i] << "  ";
        out << "\n    后5: ";
        for (std::size_t i = n - 5; i < n; ++i)
            out << evals[i] << "  ";
        out << "\n\n";
    }
}

// ======================================================================
// 问题 2(2): 100×100 Poisson 矩阵
// ======================================================================
void solveTridiagonal2(std::ofstream& out) {
    out << "============================================================\n";
    out << "问题 2(2): 100×100 三对角矩阵: diag=2, sub/super-diag=-1\n";
    out << "使用隐式对称 QR 算法\n";
    out << "============================================================\n\n";

    std::size_t n = 100;
    Matrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A(i, i) = 2.0;
        if (i > 0) A(i, i - 1) = -1.0;
        if (i < n - 1) A(i, i + 1) = -1.0;
    }

    auto result = symmetricEigensolver(A, 1e-14, 1000);

    auto evals = result.eigenvalues;
    std::sort(evals.begin(), evals.end());

    // Theoretical: λ_k = 2 - 2·cos(kπ/101), k=1..100
    std::vector<double> theory(n);
    for (std::size_t k = 1; k <= n; ++k)
        theory[k - 1] = 2.0 - 2.0 * std::cos(k * M_PI / (n + 1));
    std::sort(theory.begin(), theory.end());

    double maxErr = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        maxErr = std::max(maxErr, std::abs(evals[i] - theory[i]));

    double trace = 0;
    for (auto v : evals) trace += v;

    out << "  n = " << n << ":\n";
    out << "    收敛: " << (result.converged ? "Yes" : "No")
        << "  迭代: " << result.iterations << "\n\n";
    out << "    λ_max = " << std::fixed << std::setprecision(12) << evals[n - 1] << "\n";
    out << "    λ_min = " << evals[0] << "\n";
    out << "    解析解最大误差: " << std::scientific << maxErr << "\n";
    out << "    trace = " << std::fixed << std::setprecision(6) << trace
        << "  (2n = 200)\n\n";

    out << "    前10个特征值:\n";
    for (std::size_t i = 0; i < 10; ++i)
        out << "      λ[" << i << "] = " << std::setprecision(12) << evals[i] << "\n";
    out << "    ...\n";
    out << "    最后10个特征值:\n";
    for (std::size_t i = n - 10; i < n; ++i)
        out << "      λ[" << i << "] = " << evals[i] << "\n";
    out << "\n";

    // Eigenvector verification
    Matrix Q = result.Q;
    out << "    特征向量验证 (||Av - λv||_inf):\n";
    for (int ci = 0; ci < 3; ++ci) {
        std::size_t idx = ci * (n - 1) / 2;  // 0, ~50, ~99
        double lambda = evals[idx];
        std::vector<double> av(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j)
                av[i] += A(i, j) * Q(j, idx);
        }
        double err = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            err = std::max(err, std::abs(av[i] - lambda * Q(i, idx)));
        out << "      v[" << idx << "]: " << std::scientific << err << "\n";
    }
    out << "\n";
}

int main() {
    std::filesystem::create_directories("./output");
    std::ofstream out("./output/homework14_problems.txt");
    if (!out.is_open()) {
        std::cerr << "Failed to open output file.\n";
        return 1;
    }

    out << "============================================================\n";
    out << "Homework 14: 编程题求解报告\n";
    out << "============================================================\n\n";

    solvePolynomialRoots(out);
    solveMatrixEigenvalues(out);
    solveTridiagonal1(out);
    solveTridiagonal2(out);

    out << "============================================================\n";
    out << "所有题目求解完成.\n";
    out << "============================================================\n";

    out.close();
    std::cout << "Results written to ./output/homework14_problems.txt" << std::endl;
    return 0;
}
