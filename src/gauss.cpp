#include "matrix.hpp"
#include "operations.hpp"
#include "solve.hpp"

#include <fstream>
#include <filesystem>

// 不选主元的高斯消元法
numalg::Matrix gaussSolve(const numalg::Matrix& A,
                          const numalg::Matrix& b) {
    if (A.lines() != A.rows()) {
        throw std::invalid_argument("matrix A must be square");
    }
    if (A.lines() != b.lines()) {
        throw std::invalid_argument("matrix A and b must have the same number of lines");
    }
    std::size_t n = A.lines();
    numalg::Matrix U = A;  // 上三角矩阵
    numalg::Matrix y = b;  // 中间结果
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t i = k + 1; i < n; ++i) {
            double factor = U(i, k) / U(k, k);
            for (std::size_t j = k; j < n; ++j) {
                U(i, j) -= factor * U(k, j);
            }
            y(i, 0) -= factor * y(k, 0);
        }
    }
    return solve_upper_triangular(U, y);
}

// 列主元的高斯消元法
numalg::Matrix PgaussSolve(const numalg::Matrix& A,
                           const numalg::Matrix& b) {
    if (A.lines() != A.rows()) {
        throw std::invalid_argument("matrix A must be square");
    }
    if (A.lines() != b.lines()) {
        throw std::invalid_argument("matrix A and b must have the same number of lines");
    }
    std::size_t n = A.lines();
    numalg::Matrix U = A;
    numalg::Matrix y = b;
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t max_row = k;
        for (std::size_t i = k + 1; i < n; ++i) {
            if (abs(U(i, k)) > abs(U(max_row, k))) {
                max_row = i;
            }
        }
        if (max_row != k) {
            // 交换行
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(U(k, j), U(max_row, j));
            }
            std::swap(y(k, 0), y(max_row, 0));
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            double factor = U(i, k) / U(k, k);
            for (std::size_t j = k; j < n; ++j) {
                U(i, j) -= factor * U(k, j);
            }
            y(i, 0) -= factor * y(k, 0);
        }
    }
    return solve_upper_triangular(U, y);
}

// 平方根法 aka Cholesky 分解
numalg::Matrix choleskySolve(const numalg::Matrix& A,
                             const numalg::Matrix& b) {
    if (A.lines() != A.rows()) {
        throw std::invalid_argument("matrix A must be square");
    }
    if (A.lines() != b.lines()) {
        throw std::invalid_argument("matrix A and b must have the same number of lines");
    }
    
    /*
    for k = 1 : n
        A(k, k) = sqrt(A(k, k))
        A(k+1:n, k) = A(k+1:n, k) / A(k, k)
        for j = k+1 : n
            A(j:n, j) = A(j:n, j) - A(j:n, k) * A(j, k)
        end
    end
    so L = A 的下三角部分，A = L L^T
    */
    std::size_t n = A.lines();
    numalg::Matrix L = A;  // 下三角矩阵
    for (std::size_t k = 0; k < n; ++k) {
        L(k, k) = sqrt(L(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            L(i, k) /= L(k, k);
        }
        for (std::size_t j = k + 1; j < n; ++j) {
            for (std::size_t i = j; i < n; ++i) {
                L(i, j) -= L(i, k) * L(j, k);
            }
        }
    }
    // L = A 的下三角部分
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            L(i, j) = 0.0;
        }
    }
    // 求解 Ly = b
    numalg::Matrix y = solve_lower_triangular(L, b);
    // 求解 L^T x = y
    return solve_upper_triangular(L.transpose(), y);
}

// 改进的平方根法 let A = L D L^T
numalg::Matrix modifiedCholeskySolve(const numalg::Matrix& A,
                                     const numalg::Matrix& b) {
    if (A.lines() != A.rows()) {
        throw std::invalid_argument("matrix A must be square");
    }
    if (A.lines() != b.lines()) {
        throw std::invalid_argument("matrix A and b must have the same number of lines");
    }

    /*
    for j = 1 : n
        for i = 1 : j-1
            v(i) = A(j, i)A(i, i)
        end
        A(j, j) = A(j, j) - A(j, 1:j-1) * v(1:j-1)
        A(j+1:n, j) = (A(j+1:n, j) - A(j+1:n, 1:j-1) * v(1:j-1)) / A(j, j)
    end
    so L 的严格下三角部分 = A 的严格下三角部分，L 的对角元素 = 1，D = A 的对角线部分，A = L D L^T
    */
    std::size_t n = A.lines();
    numalg::Matrix L = A;  // 下三角矩阵
    for (std::size_t j = 0; j < n; ++j) {
        numalg::Matrix v(n);
        for (std::size_t i = 0; i < j; ++i) {
            v(i) = L(j, i) * L(i, i);
        }
        for (std::size_t i = 0; i < j; ++i) {
            L(j, j) -= L(j, i) * v(i);
        }
        for (std::size_t i = j + 1; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t k = 0; k < j; ++k) {
                sum += L(i, k) * v(k);
            }
            L(i, j) = (L(i, j) - sum) / L(j, j);
        }
    }
    // L 的严格下三角部分 = A 的严格下三角部分，L 的对角元素 = 1，D = A 的对角线部分
    numalg::Matrix D(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        D(i, i) = L(i, i);
        L(i, i) = 1.0;
        for (std::size_t j = i + 1; j < n; ++j) {
            L(i, j) = 0.0;
        }
    }
    // 求解 Ly = b
    numalg::Matrix y = solve_lower_triangular(L, b);
    // 求解 Dz = y
    for (std::size_t i = 0; i < n; ++i) {
        y(i) /= D(i, i);
    }
    // 求解 L^T x = z
    return solve_upper_triangular(L.transpose(), y);
}