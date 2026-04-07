#ifndef NUMALG_GAUSS_HPP
#define NUMALG_GAUSS_HPP

#include "matrix.hpp"

/**
 * @brief Matrix inverse via Gaussian elimination with partial pivoting
 * Perform partial-pivot Gauss elimination on [A | I] to obtain [I | A^{-1}]
 *
 * @param A invertible square matrix
 * @return numalg::Matrix the inverse of A
 */
numalg::Matrix gaussInverse(const numalg::Matrix& A);

/**
 * @brief Gaussian elimination without pivoting
 * 
 * @param A coefficient matrix
 * @param b constant terms
 * @return numalg::Matrix 
 */
numalg::Matrix gaussSolve(const numalg::Matrix& A,
                          const numalg::Matrix& b);

/**
 * @brief Gaussian elimination with partial (column) pivoting
 * 
 * @param A coefficient matrix
 * @param b constant terms
 * @return numalg::Matrix 
 */
numalg::Matrix PgaussSolve(const numalg::Matrix& A,
                           const numalg::Matrix& b);

/**
 * @brief Cholesky decomposition: A = L L^T
 * 
 * @param A coefficient matrix (must be symmetric positive definite)
 * @param b constant terms
 * @return numalg::Matrix 
 */
numalg::Matrix choleskySolve(const numalg::Matrix& A,
                             const numalg::Matrix& b);

/**
 * @brief Modified Cholesky (LDL^T) decomposition
 * Avoids square root operations
 * @param A coefficient matrix (must be symmetric positive definite)
 * @param b constant terms
 * @return numalg::Matrix
 */
numalg::Matrix modifiedCholeskySolve(const numalg::Matrix& A,
                                     const numalg::Matrix& b);

#endif  // NUMALG_GAUSS_HPP
