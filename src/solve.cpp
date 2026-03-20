#include "matrix.hpp"
#include "operations.hpp"

// 下三角矩阵求解
numalg::Matrix solve_lower_triangular(const numalg::Matrix& L,
                                      const numalg::Matrix& b) {
    if (L.is_lower_triangular() == false) {
        throw std::invalid_argument("matrix L must be lower triangular");
    }
    numalg::Matrix x(L.lines(), 1);
    for (std::size_t i = 0; i < L.lines(); ++i) {
        double sum = 0.0;
        for (std::size_t j = 0; j < i; ++j) {
            sum += L(i, j) * x(j, 0);
        }
        x(i, 0) = (b(i, 0) - sum) / L(i, i);
    }
    return x;
}

// 上三角矩阵求解
numalg::Matrix solve_upper_triangular(const numalg::Matrix& U,
                                      const numalg::Matrix& b) {
    if (U.is_upper_triangular() == false) {
        throw std::invalid_argument("matrix U must be upper triangular");
    }
    numalg::Matrix x(U.lines(), 1);
    for (std::size_t i = U.lines(); i-- > 0;) {
        double sum = 0.0;
        for (std::size_t j = i + 1; j < U.rows(); ++j) {
            sum += U(i, j) * x(j, 0);
        }
        x(i, 0) = (b(i, 0) - sum) / U(i, i);
    }
    return x;
}