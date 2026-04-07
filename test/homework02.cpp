#include <iostream>

#include "matrix.hpp"
#include "operations.hpp"
#include "solve.hpp"
#include "gauss.hpp"

int main() {
    // definition of the test matrix A and vector b
    numalg::Matrix A(84, 84);
    for (std::size_t i = 0; i < A.lines(); ++i) {
        A.at(i, i) = 6;
    }
    for (std::size_t i = 0; i < A.lines() - 1; ++i) {
        A.at(i, i + 1) = 1;
        A.at(i + 1, i) = 8;
    }

    numalg::Matrix b(84);
    for (std::size_t i = 0; i < b.lines(); ++i) {
        b.at(i) = 15;
    }
    b.at(0) = 7;
    b.at(83) = 14;

    numalg::Matrix x_1 = gaussSolve(A, b);
    numalg::Matrix x_2 = PgaussSolve(A, b);
    std::cout << "Gauss elimination x_1:\n";
    x_1.print();
    std::cout << "Gauss elimination with pivoting x_2:\n";
    x_2.print();
    return 0;
}