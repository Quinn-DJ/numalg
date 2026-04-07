#ifndef NUMALG_GAUSS_HPP
#define NUMALG_GAUSS_HPP

#include "matrix.hpp"

/**
 * @brief 列主元 Gauss 消元法求矩阵逆
 * 对 [A | I] 做列主元 Gauss 消元，得到 [I | A^{-1}]
 *
 * @param A 可逆方阵
 * @return numalg::Matrix A 的逆矩阵
 */
numalg::Matrix gaussInverse(const numalg::Matrix& A);

/**
 * @brief 不进行主元选择的高斯消元法
 * 
 * @param A 系数矩阵
 * @param b 常数项矩阵
 * @return numalg::Matrix 
 */
numalg::Matrix gaussSolve(const numalg::Matrix& A,
                          const numalg::Matrix& b);

/**
 * @brief 列主元高斯消元法
 * 
 * @param A 系数矩阵
 * @param b 常数项矩阵
 * @return numalg::Matrix 
 */
numalg::Matrix PgaussSolve(const numalg::Matrix& A,
                           const numalg::Matrix& b);

/**
 * @brief 平方根法 aka Cholesky 分解
 * 
 * @param A 系数矩阵，必须是对称正定矩阵
 * @param b 常数项矩阵
 * @return numalg::Matrix 
 */
numalg::Matrix choleskySolve(const numalg::Matrix& A,
                             const numalg::Matrix& b);

/**
 * @brief 改进的平方根法
 * 避免了开方
 * @param A 系数矩阵，必须是对称正定矩阵
 * @param b 常数项矩阵
 * @return numalg::Matrix
 */
numalg::Matrix modifiedCholeskySolve(const numalg::Matrix& A,
                                     const numalg::Matrix& b);

#endif  // NUMALG_GAUSS_HPP
