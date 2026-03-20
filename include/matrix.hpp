#ifndef NUMALG_MATRIX_HPP
#define NUMALG_MATRIX_HPP

#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>

namespace numalg {

class Matrix {
public:

    Matrix(std::size_t lines, std::size_t rows)
        : lines_(lines), rows_(rows), data_(lines * rows, 0.0) {}

    Matrix(std::size_t lines, std::size_t rows, std::initializer_list<double> init)
        : lines_(lines), rows_(rows), data_(init) {
        if (data_.size() != lines_ * rows_) {
            throw std::invalid_argument("initializer_list size mismatch");
        }
    }

    Matrix(std::size_t lines)
        : lines_(lines), rows_((std::size_t)1), data_(lines, 0.0) {}

    Matrix(std::size_t lines, std::initializer_list<double> init)
        : lines_(lines), rows_((std::size_t)1), data_(init) {
        if (data_.size() != lines_) {
            throw std::invalid_argument("initializer_list size mismatch");
        }
    }

    std::size_t lines() const { return lines_; }
    std::size_t rows() const { return rows_; }

    double& operator()(std::size_t i, std::size_t j) { return data_[index(i, j)]; }
    double operator()(std::size_t i, std::size_t j) const { return data_[index(i, j)]; }
    double& operator()(std::size_t i) { return data_[index(i, 0)]; }
    double operator()(std::size_t i) const { return data_[index(i, 0)]; }

    double& at(std::size_t i, std::size_t j) {
        check_bounds(i, j);
        return data_[index(i, j)];
    }

    double at(std::size_t i, std::size_t j) const {
        check_bounds(i, j);
        return data_[index(i, j)];
    }

    double& at(std::size_t i) {
        if (rows_ != 1) {
            throw std::invalid_argument("matrix is not a vector");
        }
        check_bounds(i, 0);
        return data_[i];
    }

    double at(std::size_t i) const {
        if (rows_ != 1) {
            throw std::invalid_argument("matrix is not a vector");
        }
        check_bounds(i, 0);
        return data_[i];
    }

    const std::vector<double>& data() const { return data_; }

    // 下三角矩阵检测
    bool is_lower_triangular() const {
        if (lines_ != rows_) {
            return false;
        }
        for (std::size_t i = 0; i < lines_; ++i) {
            for (std::size_t j = i + 1; j < rows_; ++j) {
                if (abs((*this)(i, j)) > 1e-9) {
                    return false;
                }
            }
        }
        return true;
    }

    // 上三角矩阵检测
    bool is_upper_triangular() const {
        if (lines_ != rows_) {
            return false;
        }
        for (std::size_t i = 1; i < lines_; ++i) {
            for (std::size_t j = 0; j < i; ++j) {
                if (abs((*this)(i, j)) > 1e-9) {
                    return false;
                }
            }
        }
        return true;
    }

    void print() const {
        for (std::size_t i = 0; i < lines_; ++i) {
            for (std::size_t j = 0; j < rows_; ++j) {
                std::cout << (*this)(i, j) << " ";
            }
            std::cout << "\n";
        }
    }

    // transpose
    Matrix transpose() const {
        Matrix result(rows_, lines_);
        for (std::size_t i = 0; i < lines_; ++i) {
            for (std::size_t j = 0; j < rows_; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

private:
    std::size_t lines_;
    std::size_t rows_;
    std::vector<double> data_;

    std::size_t index(std::size_t i, std::size_t j) const {
        return i * rows_ + j;
    }

    void check_bounds(std::size_t i, std::size_t j) const {
        if (i >= lines_ || j >= rows_) {
            throw std::out_of_range("matrix index out of range");
        }
    }
};

}  // namespace numalg

#endif  // NUMALG_MATRIX_HPP
