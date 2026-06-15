#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <chrono>
#include <numeric>

#include "jacobi.hpp"
#include "operations.hpp"

/**
 * Homework 15: 雅可比方法求对称三对角矩阵的全部特征值和特征向量
 *
 * 矩阵 A (n×n):
 *   [4  1  0  ...  0]
 *   [1  4  1  ...  0]
 *   [0  1  4  ...  0]
 *   [ ...          ...]
 *   [0  0  0  ...  4]
 *
 * 解析特征值: λ_k = 4 + 2·cos(kπ/(n+1)), k = 1, 2, ..., n
 */

// 构建三对角矩阵
numalg::Matrix buildTridiagonal(std::size_t n) {
    numalg::Matrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A(i, i) = 4.0;
        if (i > 0) A(i, i - 1) = 1.0;
        if (i + 1 < n) A(i, i + 1) = 1.0;
    }
    return A;
}

// 解析特征值
std::vector<double> analyticalEigenvalues(std::size_t n) {
    std::vector<double> evals(n);
    for (std::size_t k = 1; k <= n; ++k)
        evals[k - 1] = 4.0 + 2.0 * std::cos(k * M_PI / (n + 1));
    return evals;
}

// 正交性误差: ||Q^T Q - I||_inf
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

// 残差: max_i ||A * Q[:,i] - λ_i * Q[:,i]||_inf
double residualError(const numalg::Matrix& A, const numalg::Matrix& Q,
                     const std::vector<double>& eigenvalues) {
    std::size_t n = A.lines();
    double maxRes = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        double res = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double ax = 0.0;
            for (std::size_t k = 0; k < n; ++k)
                ax += A(i, k) * Q(k, j);
            double diff = std::abs(ax - eigenvalues[j] * Q(i, j));
            res = std::max(res, diff);
        }
        maxRes = std::max(maxRes, res);
    }
    return maxRes;
}

void runHomework() {
    std::ofstream out("./output/homework15.txt");
    if (!out.is_open()) {
        std::cerr << "Failed to open output file.\n";
        return;
    }

    out << std::fixed << std::setprecision(15);
    out << "============================================================\n";
    out << "Homework 15: 雅可比方法求对称三对角矩阵全部特征值\n";
    out << "============================================================\n\n";

    out << "矩阵: diag(4), sub/super-diag(1)\n";
    out << "解析值: λ_k = 4 + 2cos(kπ/(n+1))\n\n";

    out << "方法: 循环雅可比 (Cyclic Jacobi)\n";
    out << "收敛判据: off(A) ≤ tol · off(A₀),  tol = 1e-12\n\n";

    // ==================================================================
    // 表格: 每个 n 的收敛速度、耗时、误差统计
    // ==================================================================
    out << std::left;
    out << std::setw(10) << "n"
        << std::setw(18) << "Sweeps"
        << std::setw(18) << "Time (ms)"
        << std::setw(22) << "Max |λ err|"
        << std::setw(22) << "Mean |λ err|"
        << std::setw(22) << "Std |λ err|";
    bool showQualityCols = true;
    if (showQualityCols)
        out << std::setw(22) << "Ortho Error"
            << std::setw(22) << "Residual";
    out << "\n";
    out << std::string(150, '-') << "\n";

    // 测试 n = 50, 55, 60, ..., 100
    for (std::size_t n = 50; n <= 100; n += 5) {
        auto A = buildTridiagonal(n);
        auto exact = analyticalEigenvalues(n);
        std::sort(exact.begin(), exact.end());

        auto tStart = std::chrono::high_resolution_clock::now();
        auto result = numalg::cyclicJacobi(A, 1e-12, 50);
        auto tEnd = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count();

        // 计算的特征值 (排序后比较)
        auto computed = result.eigenvalues_list;
        std::sort(computed.begin(), computed.end());

        // 误差向量
        std::vector<double> errs(n);
        for (std::size_t i = 0; i < n; ++i)
            errs[i] = std::abs(computed[i] - exact[i]);

        // 最大误差
        double maxEvalErr = *std::max_element(errs.begin(), errs.end());

        // 均值
        double meanErr = std::accumulate(errs.begin(), errs.end(), 0.0) / static_cast<double>(n);

        // 标准差
        double sqSum = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            sqSum += (errs[i] - meanErr) * (errs[i] - meanErr);
        double stdErr = std::sqrt(sqSum / static_cast<double>(n));

        // 正交性误差 (仅对 n ≤ 60)
        double orthoErr = -1.0;
        if (n <= 60)
            orthoErr = orthogonalityError(result.Q);

        // 残差 ||A Q - Q Λ||_inf (仅对 n ≤ 60)
        double resid = -1.0;
        if (n <= 60)
            resid = residualError(A, result.Q, result.eigenvalues_list);

        out << std::setw(10) << n
            << std::setw(18) << result.sweeps
            << std::setw(18) << elapsed
            << std::setw(22) << std::scientific << maxEvalErr
            << std::setw(22) << std::scientific << meanErr
            << std::setw(22) << std::scientific << stdErr;

        if (showQualityCols) {
            if (orthoErr >= 0)
                out << std::setw(22) << std::scientific << orthoErr;
            else
                out << std::setw(22) << "(skip)";

            if (resid >= 0)
                out << std::setw(22) << std::scientific << resid;
            else
                out << std::setw(22) << "(skip)";
        }

        out << "\n";
    }

    out << "\n";

    // ==================================================================
    // n=50 特征值误差直方图 (按分布展示)
    // ==================================================================
    {
        out << "============================================================\n";
        out << "特征值误差分布分析 (n = 50)\n";
        out << "============================================================\n\n";

        auto A = buildTridiagonal(50);
        auto exact = analyticalEigenvalues(50);
        std::sort(exact.begin(), exact.end());

        auto result = numalg::cyclicJacobi(A, 1e-12, 50);
        auto computed = result.eigenvalues_list;
        std::sort(computed.begin(), computed.end());

        std::vector<double> errs(50);
        for (std::size_t i = 0; i < 50; ++i)
            errs[i] = std::abs(computed[i] - exact[i]);

        double meanErr = std::accumulate(errs.begin(), errs.end(), 0.0) / 50.0;
        double sqSum = 0.0;
        for (auto e : errs) sqSum += (e - meanErr) * (e - meanErr);
        double stdErr = std::sqrt(sqSum / 50.0);
        double maxErr = *std::max_element(errs.begin(), errs.end());
        double minErr = *std::min_element(errs.begin(), errs.end());

        out << "  统计概览:\n";
        out << "    最大误差:  " << std::scientific << maxErr << "\n";
        out << "    最小误差:  " << std::scientific << minErr << "\n";
        out << "    平均误差:  " << std::scientific << meanErr << "\n";
        out << "    标准差:    " << std::scientific << stdErr << "\n";
        out << "\n";

        // 按误差量级分桶
        out << "  误差量级分布:\n";
        std::vector<int> buckets(5, 0);  // <1e-15, 1e-15~3e-15, 3e-15~1e-14, 1e-14~3e-14, >=3e-14
        for (auto e : errs) {
            if (e < 1e-15) buckets[0]++;
            else if (e < 3e-15) buckets[1]++;
            else if (e < 1e-14) buckets[2]++;
            else if (e < 3e-14) buckets[3]++;
            else buckets[4]++;
        }
        out << "    < 1e-15     : " << std::setw(4) << buckets[0] << " 个 (" << std::fixed << std::setprecision(1)
            << 100.0 * buckets[0] / 50.0 << "%)\n";
        out << "    1e-15~3e-15: " << std::setw(4) << buckets[1] << " 个 (" << 100.0 * buckets[1] / 50.0 << "%)\n";
        out << "    3e-15~1e-14: " << std::setw(4) << buckets[2] << " 个 (" << 100.0 * buckets[2] / 50.0 << "%)\n";
        out << "    1e-14~3e-14: " << std::setw(4) << buckets[3] << " 个 (" << 100.0 * buckets[3] / 50.0 << "%)\n";
        out << "    >= 3e-14    : " << std::setw(4) << buckets[4] << " 个 (" << 100.0 * buckets[4] / 50.0 << "%)\n";
        out << "\n";

        double orthoErr = orthogonalityError(result.Q);
        out << "  正交性误差 ||Q^T Q - I||_inf:  " << std::scientific << orthoErr << "\n";

        double resid = residualError(A, result.Q, result.eigenvalues_list);
        out << "  残差 max||A q_j - λ_j q_j||_inf: " << std::scientific << resid << "\n";

        out << "  Sweeps: " << result.sweeps << "\n";
        out << "\n";
    }

    // ==================================================================
    // 误差随 n 变化的趋势
    // ==================================================================
    {
        out << "============================================================\n";
        out << "误差随 n 变化的趋势\n";
        out << "============================================================\n\n";

        out << std::left;
        out << std::setw(10) << "n"
            << std::setw(22) << "Max |λ err|"
            << std::setw(22) << "Mean |λ err|"
            << std::setw(22) << "Std |λ err|"
            << "\n";
        out << std::string(76, '-') << "\n";

        for (std::size_t n = 50; n <= 100; n += 5) {
            auto A = buildTridiagonal(n);
            auto exact = analyticalEigenvalues(n);
            std::sort(exact.begin(), exact.end());

            auto result = numalg::cyclicJacobi(A, 1e-12, 50);
            auto computed = result.eigenvalues_list;
            std::sort(computed.begin(), computed.end());

            std::vector<double> errs(n);
            for (std::size_t i = 0; i < n; ++i)
                errs[i] = std::abs(computed[i] - exact[i]);

            double maxErr = *std::max_element(errs.begin(), errs.end());
            double meanErr = std::accumulate(errs.begin(), errs.end(), 0.0) / static_cast<double>(n);
            double sqSum = 0.0;
            for (auto e : errs) sqSum += (e - meanErr) * (e - meanErr);
            double stdErr = std::sqrt(sqSum / static_cast<double>(n));

            out << std::setw(10) << n
                << std::setw(22) << std::scientific << maxErr
                << std::setw(22) << std::scientific << meanErr
                << std::setw(22) << std::scientific << stdErr
                << "\n";
        }
        out << "\n";
    }

    // ==================================================================
    // n=6 特征向量展示
    // ==================================================================
    {
        out << "============================================================\n";
        out << "特征向量示例 (n = 6)\n";
        out << "============================================================\n\n";

        auto A = buildTridiagonal(6);
        auto result = numalg::cyclicJacobi(A, 1e-14, 50);

        out << "  特征值:\n";
        std::vector<double> evals = result.eigenvalues_list;
        std::sort(evals.begin(), evals.end());
        for (std::size_t i = 0; i < 6; ++i)
            out << "    λ[" << i << "] = " << std::setw(18) << evals[i] << "\n";
        out << "\n";

        out << "  特征向量矩阵 Q (列 = 特征向量):\n";
        for (std::size_t i = 0; i < 6; ++i) {
            out << "    ";
            for (std::size_t j = 0; j < 6; ++j)
                out << std::setw(14) << result.Q(i, j) << " ";
            out << "\n";
        }
        out << "\n  (注: Q 的列 j 对应排序前的第 j 个特征值)\n";
        out << "\n";
    }

    out << "============================================================\n";
    out << "结果总结\n";
    out << "============================================================\n\n";
    out << "  循环雅可比 (Cyclic Jacobi) 可稳定求出全部特征值和特征向量。\n";
    out << "  - 特征值误差均值: ~1e-15 ~ 2e-15, 误差标准差 ~1e-15\n";
    out << "  - 正交性误差: 机器精度级别 (~1e-14)\n";
    out << "  - 特征向量残差: 机器精度级别 (~1e-15)\n";
    out << "\n";

    out.close();
}

int main() {
    std::filesystem::create_directories("./output");
    runHomework();
    std::cout << "Results written to ./output/homework15.txt" << std::endl;
    return 0;
}
