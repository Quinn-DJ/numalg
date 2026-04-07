#include <iostream>
#include <cmath>
#include <vector>

#include "matrix.hpp"
#include "operations.hpp"
#include "inverse_norm_estimate.hpp"

// 辅助：手动计算矩阵无穷范数（用于对比验证）
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
    std::cout << "=== 测试：矩阵逆的无穷范数估计（优化法/Hager's method） ===\n\n";

    // 测试1：已知逆的矩阵
    {
        std::cout << "--- 测试1：3x3 矩阵 ---\n";
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

        std::cout << "A^{-1} 的精确无穷范数: " << exact << "\n";
        std::cout << "估计值:                " << estimate << "\n";
        std::cout << "相对误差:              " << std::abs(estimate - exact) / exact << "\n\n";
    }

    // 测试2：5x5 矩阵
    {
        std::cout << "--- 测试2：5x5 对角占优矩阵 ---\n";
        numalg::Matrix A(5, 5, {4, 1, 0, 0, 0,
                                  1, 4, 1, 0, 0,
                                  0, 1, 4, 1, 0,
                                  0, 0, 1, 4, 1,
                                  0, 0, 0, 1, 4});

        // 手动求逆太麻烦，直接用 B = I（相当于 A = I 的情况）作为简单验证
        std::cout << "(跳过精确对比，仅验证不崩溃)\n";
        // 用 B = A 本身测试算法能正常运行
        numalg::Matrix B = A.transpose();
        double estimate = numalg::estimate_inverse_norm_inf(B);
        std::cout << "B 的 ||B^{-1}||_inf 估计值: " << estimate << "\n\n";
    }

    // 测试3：单位矩阵
    {
        std::cout << "--- 测试3：3x3 单位矩阵 ---\n";
        numalg::Matrix I(3, 3, {1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 1});
        // I^{-1} = I, ||I||_inf = 1
        double estimate = numalg::estimate_inverse_norm_inf(I);
        std::cout << "估计值: " << estimate << " (期望: 1.0)\n\n";
    }

    std::cout << "=== 测试完成 ===\n";
    return 0;
}
