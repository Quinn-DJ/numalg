#include "matrix.hpp"
#include "operations.hpp"
#include "solve.hpp"


// Gaussian elimination without pivoting
numalg::Matrix gaussSolve(const numalg::Matrix& A,
                          const numalg::Matrix& b) {
    if (A.lines() != A.rows()) {
        throw std::invalid_argument("matrix A must be square");
    }
    if (A.lines() != b.lines()) {
        throw std::invalid_argument("matrix A and b must have the same number of lines");
    }
    std::size_t n = A.lines();
    numalg::Matrix U = A;  // upper triangular matrix
    numalg::Matrix y = b;  // intermediate result
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

// Gaussian elimination with partial (column) pivoting
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
            if (std::abs(U(i, k)) > std::abs(U(max_row, k))) {
                max_row = i;
            }
        }
        if (max_row != k) {
            // swap rows
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

// Cholesky decomposition: A = L L^T
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
    so L = lower triangular part of A, A = L L^T
    */
    std::size_t n = A.lines();
    numalg::Matrix L = A;  // will become lower triangular matrix
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
    // Extract lower triangular part of L
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            L(i, j) = 0.0;
        }
    }
    // Solve Ly = b
    numalg::Matrix y = solve_lower_triangular(L, b);
    // Solve L^T x = y
    return solve_upper_triangular(L.transpose(), y);
}

// Modified Cholesky (LDL^T) decomposition: A = L D L^T
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
    so strict lower part of L = strict lower part of A, diag(L) = 1, D = diag(A), A = L D L^T
    */
    std::size_t n = A.lines();
    numalg::Matrix L = A;  // will store L and D
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
    // Extract L (unit lower triangular) and D (diagonal)
    numalg::Matrix D(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        D(i, i) = L(i, i);
        L(i, i) = 1.0;
        for (std::size_t j = i + 1; j < n; ++j) {
            L(i, j) = 0.0;
        }
    }
    // Solve Ly = b
    numalg::Matrix y = solve_lower_triangular(L, b);
    // Solve Dz = y
    for (std::size_t i = 0; i < n; ++i) {
        y(i) /= D(i, i);
    }
    // Solve L^T x = z
    return solve_upper_triangular(L.transpose(), y);
}
// Matrix inverse via Gaussian elimination with partial pivoting
// Augmented matrix [A | I] -> [I | A^{-1}]
numalg::Matrix gaussInverse(const numalg::Matrix& A) {
    if (A.lines() != A.rows()) {
        throw std::invalid_argument("matrix A must be square");
    }
    std::size_t n = A.lines();

    // Build augmented matrix [A | I], n rows x 2n columns
    numalg::Matrix aug(n, 2 * n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            aug(i, j) = A(i, j);
        }
        aug(i, n + i) = 1.0;
    }

    // Forward elimination with partial pivoting
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t max_row = k;
        for (std::size_t i = k + 1; i < n; ++i) {
            if (std::abs(aug(i, k)) > std::abs(aug(max_row, k))) {
                max_row = i;
            }
        }
        if (max_row != k) {
            for (std::size_t j = 0; j < 2 * n; ++j) {
                std::swap(aug(k, j), aug(max_row, j));
            }
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            double factor = aug(i, k) / aug(k, k);
            for (std::size_t j = k; j < 2 * n; ++j) {
                aug(i, j) -= factor * aug(k, j);
            }
        }
    }

    // Back substitution
    for (std::size_t i = n; i-- > 0;) {
        double pivot = aug(i, i);
        for (std::size_t j = 0; j < 2 * n; ++j) {
            aug(i, j) /= pivot;
        }
        for (std::size_t k = 0; k < i; ++k) {
            double factor = aug(k, i);
            for (std::size_t j = 0; j < 2 * n; ++j) {
                aug(k, j) -= factor * aug(i, j);
            }
        }
    }

    // Extract inverse from the right half
    numalg::Matrix inv(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            inv(i, j) = aug(i, n + j);
        }
    }
    return inv;
}
