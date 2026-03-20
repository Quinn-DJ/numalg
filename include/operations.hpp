#ifndef NUMALG_OPERATIONS_HPP
#define NUMALG_OPERATIONS_HPP

#include "matrix.hpp"

namespace numalg {

inline Matrix operator+(const Matrix& a, const Matrix& b) {
    if (a.lines() != b.lines() || a.rows() != b.rows()) {
        throw std::invalid_argument("matrix size mismatch for addition");
    }
    Matrix result(a.lines(), a.rows());
    for (std::size_t i = 0; i < a.lines(); ++i) {
        for (std::size_t j = 0; j < a.rows(); ++j) {
            result(i, j) = a(i, j) + b(i, j);
        }
    }
    return result;
}

inline Matrix operator-(const Matrix& a, const Matrix& b) {
    if (a.lines() != b.lines() || a.rows() != b.rows()) {
        throw std::invalid_argument("matrix size mismatch for subtraction");
    }
    Matrix result(a.lines(), a.rows());
    for (std::size_t i = 0; i < a.lines(); ++i) {
        for (std::size_t j = 0; j < a.rows(); ++j) {
            result(i, j) = a(i, j) - b(i, j);
        }
    }
    return result;
}

inline Matrix operator*(const Matrix& a, const Matrix& b) {
    if (a.rows() != b.lines()) {
        throw std::invalid_argument("matrix size mismatch for multiplication");
    }
    Matrix result(a.lines(), b.rows());
    for (std::size_t i = 0; i < a.lines(); ++i) {
        for (std::size_t j = 0; j < b.rows(); ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < a.rows(); ++k) {
                sum += a(i, k) * b(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

inline Matrix operator*(double scalar, const Matrix& a) {
    Matrix result(a.lines(), a.rows());
    for (std::size_t i = 0; i < a.lines(); ++i) {
        for (std::size_t j = 0; j < a.rows(); ++j) {
            result(i, j) = scalar * a(i, j);
        }
    }
    return result;
}

inline Matrix tr(const Matrix& a) {
    Matrix result(a.rows(), a.lines());
    for (std::size_t i = 0; i < a.lines(); ++i) {
        for (std::size_t j = 0; j < a.rows(); ++j) {
            result(j, i) = a(i, j);
        }
    }
    return result;
}

inline bool operator==(const Matrix& a, const Matrix& b) {
    if (a.lines() != b.lines() || a.rows() != b.rows()) {
        return false;
    }
    for (std::size_t i = 0; i < a.lines(); ++i) {
        for (std::size_t j = 0; j < a.rows(); ++j) {
            if (abs(a(i, j) - b(i, j)) > 1e-9) {
                return false;
            }
        }
    }
    return true;
}

}  // namespace numalg

#endif  // NUMALG_OPERATIONS_HPP
