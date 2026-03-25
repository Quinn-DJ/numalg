#include <iostream>

#include "matrix.hpp"
#include "operations.hpp"
#include "solve.hpp"
#include "gauss.hpp"

int main() {
    // 随机数种子
    srand(static_cast<unsigned int>(time(nullptr)));

    // 测试矩阵类
    numalg::Matrix A(10, 10);
    for (std::size_t i = 0; i < A.lines(); ++i) {
        for (std::size_t j = 0; j < A.rows(); ++j) {
            A.at(i, j) = rand() % 1000;  // 填充随机数
        }
    }
    A.print();

    numalg::Matrix b(10, 1);
    for (std::size_t i = 0; i < b.lines(); ++i) {
        b.at(i, 0) = rand() % 1000;  // 填充随机数
    }
    b.print();

    // 测试高斯消元法
    try {
        numalg::Matrix x = gaussSolve(A, b);
        std::cout << "Solution:\n";
        // for (std::size_t i = 0; i < x.lines(); ++i) {
        //     std::cout << x(i, 0) << "\n";
        // }
        if (A*x == b) {
            std::cout << "Verification passed: A * x == b\n";
        } else {
            std::cout << "Verification failed: A * x != b\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    // 测试列主元的高斯消元法
    try {
        numalg::Matrix x_pivot = PgaussSolve(A, b);
        std::cout << "Solution with pivoting:\n";
        // for (std::size_t i = 0; i < x_pivot.lines(); ++i) {
        //     std::cout << x_pivot(i, 0) << "\n";
        // }
        if (A*x_pivot == b) {
            std::cout << "Verification passed: A * x_pivot == b\n";
        } else {
            std::cout << "Verification failed: A * x_pivot != b\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    // 生成一个对称正定矩阵 A_symm
    numalg::Matrix A_symm(10, 10);
    for (std::size_t i = 0; i < A_symm.lines(); ++i) {
        A_symm.at(i, i) = 10;
        if (i > 0) {
            A_symm.at(i, i - 1) = 1;
            A_symm.at(i - 1, i) = 1;
        }
    }
    A_symm.print();

    // 测试平方根法
    try {
        numalg::Matrix x_cholesky = choleskySolve(A_symm, b);
        std::cout << "Solution with Cholesky:\n";
        // for (std::size_t i = 0; i < x_cholesky.lines(); ++i) {
        //     std::cout << x_cholesky(i, 0) << "\n";
        // }
        if (A_symm*x_cholesky == b) {
            std::cout << "Verification passed: A_symm * x_cholesky == b\n";
        } else {
            std::cout << "Verification failed: A_symm * x_cholesky != b\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    // 测试改进的平方根法
    try {
        numalg::Matrix x_modified_cholesky = modifiedCholeskySolve(A_symm, b);
        std::cout << "Solution with modified Cholesky:\n";
        // for (std::size_t i = 0; i < x_modified_cholesky.lines(); ++i) {
        //     std::cout << x_modified_cholesky(i, 0) << "\n";
        // }
        if (A_symm*x_modified_cholesky == b) {
            std::cout << "Verification passed: A_symm * x_modified_cholesky == b\n";
        } else {
            std::cout << "Verification failed: A_symm * x_modified_cholesky != b\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
    return 0;
}