/**
 * Homework 5 - 解的精度估计
 *
 * 题目1: 估计 5~20 阶 Hilbert 矩阵的无穷范数条件数
 *   条件数 cond(A)_inf = ||A||_inf * ||A^{-1}||_inf
 *   其中 ||A^{-1}||_inf 通过 Hager 优化法估计
 *
 * 题目2: 对特殊下三角矩阵 A_n（n=5~30），随机生成 x，计算 b=Ax，
 *   用列主元 Gauss 消元法求解，估计计算解的相对误差并与真实相对误差比较
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

// 辅助：构造 Hilbert 矩阵 H(i,j) = 1/(i+j+1)
numalg::Matrix hilbert(std::size_t n) {
    numalg::Matrix H(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            H(i, j) = 1.0 / (i + j + 1.0);
        }
    }
    return H;
}

// 辅助：构造题目2中的 A_n 矩阵
// 对角线为1，第n列为1，下三角为-1，其余为0
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

// 辅助：随机生成 n 维向量，分量在 [-1, 1]
numalg::Matrix random_vector(std::size_t n) {
    numalg::Matrix x(n, 1);
    for (std::size_t i = 0; i < n; ++i) {
        x(i, 0) = 2.0 * rand() / RAND_MAX - 1.0;
    }
    return x;
}

// 辅助：向量无穷范数
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
    // 题目1: Hilbert 矩阵的无穷范数条件数
    // =============================================
    // cond(A)_inf = ||A||_inf * ||A^{-1}||_inf
    // 令 B = A^{-T}，则 estimate_inverse_norm_inf(B) = ||A^{-1}||_inf
    std::cout << "===== 题目1: Hilbert 矩阵无穷范数条件数 =====\n";
    std::cout << "  n       ||H||_inf     ||H^{-1}||_inf  cond(H)_inf\n";
    for (std::size_t n = 5; n <= 20; ++n) {
        numalg::Matrix H = hilbert(n);

        // ||H||_inf：直接用矩阵类的方法
        double norm_H = H.normInf();

        // ||H^{-1}||_inf：先求逆，再直接算无穷范数
        // H^{-T} = (H^{-1})^T，由于 H 对称，H^{-1} 也对称
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
    // 题目2: 特殊矩阵 A_n 的计算解精度估计
    // =============================================
    // 对 n=5~30:
    //   1. 随机生成 x
    //   2. 计算 b = A_n * x
    //   3. 用列主元 Gauss 消元法求解 x_tilde
    //   4. 估计相对误差:
    //      - v_hat = estimate_inverse_norm_inf(A_n^T)  (即 ||A^{-1}||_inf)
    //      - mu = ||A||_inf
    //      - gamma = ||r||_inf, r = b - A*x_tilde
    //      - beta = ||b||_inf
    //      - rho = v_hat * mu * gamma / beta  (估计的相对误差)
    //   5. 真实相对误差: ||x - x_tilde||_inf / ||x||_inf
    std::cout << "===== 题目2: A_n 计算解精度估计 =====\n";
    std::cout << "   n   估计相对误差    真实相对误差\n";
    for (std::size_t n = 5; n <= 30; ++n) {
        numalg::Matrix A = make_An(n);

        // 随机生成真实解 x
        numalg::Matrix x = random_vector(n);

        // 计算 b = A * x
        numalg::Matrix b = A * x;

        // 列主元 Gauss 消元法求解
        numalg::Matrix x_tilde = PgaussSolve(A, b);

        // 残差 r = b - A * x_tilde
        numalg::Matrix r = b - A * x_tilde;

        // 各范数
        double v_hat = numalg::estimate_inverse_norm_inf(A.transpose());
        double mu = A.normInf();
        double gamma = vec_inf_norm(r);
        double beta = vec_inf_norm(b);

        // 估计的相对误差
        double rho = (beta != 0.0) ? v_hat * mu * gamma / beta : 0.0;

        // 真实相对误差
        double err = vec_inf_norm(x - x_tilde) / vec_inf_norm(x);

        std::cout << "  " << std::setw(2) << n << "   "
                  << std::setw(14) << rho << "  "
                  << std::setw(14) << err << "\n";
    }

    return 0;
}
