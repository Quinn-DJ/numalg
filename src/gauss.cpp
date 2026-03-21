#include "matrix.hpp"
#include "operations.hpp"
#include "solve.hpp"

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

    std::size_t n = A.lines();
    numalg::Matrix L(n, n);  // 下三角矩阵
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < j; ++k) {
                sum += L(i, k) * L(j, k);
            }
            if (i == j) {
                L(i, j) = std::sqrt(A(i, j) - sum);
            } else {
                L(i, j) = (A(i, j) - sum) / L(j, j);
            }
        }
    }
    // 求解 Ly = b
    numalg::Matrix y = solve_lower_triangular(L, b);
    // 求解 L^T x = y
    return solve_upper_triangular(L.transpose(), y);
}

// 上三角平方根法 A = U^T U
numalg::Matrix upperCholeskySolve(const numalg::Matrix& A,
                                  const numalg::Matrix& b) {
    if (A.lines() != A.rows()) {
        throw std::invalid_argument("matrix A must be square");
    }
    if (A.lines() != b.lines()) {
        throw std::invalid_argument("matrix A and b must have the same number of lines");
    }

    std::size_t n = A.lines();
    numalg::Matrix U(n, n);  // 上三角矩阵
    for (std::size_t j = 0; j < n; ++j) {
        // 计算对角元素 U(j,j)
        double sum = 0.0;
        for (std::size_t k = 0; k < j; ++k) {
            sum += U(k, j) * U(k, j);
        }
        U(j, j) = std::sqrt(A(j, j) - sum);

        // 计算 U(j, i)，i > j，写入第 j 行，内存连续
        for (std::size_t i = j + 1; i < n; ++i) {
            double s = 0.0;
            for (std::size_t k = 0; k < j; ++k) {
                s += U(k, j) * U(k, i);
            }
            U(j, i) = (A(j, i) - s) / U(j, j);
        }
    }
    // 求解 U^T y = b（下三角矩阵求解）
    numalg::Matrix y = solve_lower_triangular(U.transpose(), b);
    // 求解 U x = y（上三角矩阵求解）
    return solve_upper_triangular(U, y);
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

    std::size_t n = A.lines();
    numalg::Matrix L(n, n);  // 下三角矩阵
    numalg::Matrix D(n, n);  // 对角矩阵
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < j; ++k) {
                sum += L(i, k) * D(k, k) * L(j, k);
            }
            if (i == j) {
                D(i, i) = A(i, i) - sum;
                L(i, j) = 1.0;
            } else {
                L(i, j) = (A(i, j) - sum) / D(j, j);
            }
        }
    }
    // 求解 Ly = b
    numalg::Matrix y = solve_lower_triangular(L, b);
    // 求解 (D L^T) x = y
    // 由于 D 是对角矩阵，所以先让 y 除以 D 的对角元素，再求解 L^T x = y 以减少计算量
    for (std::size_t i = 0; i < n; ++i) {
        y(i, 0) /= D(i, i);
    } 
    return solve_upper_triangular(L.transpose(), y);
}