#ifndef NUMALG_GAUSS_HPP
#define NUMALG_GAUSS_HPP

#include "matrix.hpp"
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
 * @brief 上三角平方根法 (A = U^T U)
 * 分解为上三角矩阵 U，内层循环按行访问，对行主序存储更缓存友好
 * @param A 系数矩阵，必须是对称正定矩阵
 * @param b 常数项矩阵
 * @return numalg::Matrix
 */
numalg::Matrix upperCholeskySolve(const numalg::Matrix& A,
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