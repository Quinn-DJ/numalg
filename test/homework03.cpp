#include <iostream>
#include <fstream>
#include <filesystem>

#include "matrix.hpp"
#include "operations.hpp"
#include "solve.hpp"
#include "gauss.hpp"

void prob1(std::ofstream& outfile) {
    numalg::Matrix A(100, 100);
    for (int i = 0; i < 100; ++i) {
        A(i, i) = 10;
        if (i > 0) {
            A(i, i - 1) = 1;
        }
        if (i < 99) {
            A(i, i + 1) = 1;
        }
    }
    std::initializer_list<double> b_values = {
         34, -55, 46, 75, 4, 26, -58, 30, 39, 81, 74, 0, 48, 6, -38, 18, 27, 39, -27, -97, -31, 56, -75, 98, 72, 31, -59, -7, -71, 76, 84, -60, 42, 55, 84, -96, -96, 53, -67, 19, 0, 80, -31, 94, 84, 81, -40, -95, -7, 10, -16, -23, 64, -71, -84, -45, -40, 74, 42, 96, 100, -19, -72, 40, -13, -83, 89, 34, 28, 30, -57, 17, -58, -3, -95, -27, 88, -7, -68, -88, -18, -62, -21, -26, -38, -15, 78, 64, -4, 58, -27, 2, -46, 67, 57, -1, -59, 73, -77, 63
    };
    numalg::Matrix b = numalg::Matrix(100, b_values);

    outfile << "Prob1 Vector b:\n";
    for (int i = 0; i < 100; ++i) {
        outfile << b(i) << " ";
    }
    outfile << "\n";

    numalg::Matrix x_cholesky = choleskySolve(A, b);
    outfile << "Prob1 Cholesky solution:\n";
    for (int i = 0; i < 100; ++i) {
        outfile << x_cholesky(i) << " ";
    }
    outfile << "\n";

    numalg::Matrix x_modified_cholesky = modifiedCholeskySolve(A, b);
    outfile << "Prob1 Modified Cholesky solution:\n";
    for (int i = 0; i < 100; ++i) {
        outfile << x_modified_cholesky(i) << " ";
    }
    outfile << "\n";

    numalg::Matrix x_gauss = gaussSolve(A, b);
    outfile << "Prob1 Gaussian elimination solution:\n";
    for (int i = 0; i < 100; ++i) {
        outfile << x_gauss(i) << " ";
    }
    outfile << "\n";

    numalg::Matrix x_pGauss = PgaussSolve(A, b);
    outfile << "Prob1 Partial pivoting Gaussian elimination solution:\n";
    for (int i = 0; i < 100; ++i) {
        outfile << x_pGauss(i) << " ";
    }
    outfile << "\n";
}

void prob2(std::ofstream& outfile) {
    int n = 20;
    numalg::Matrix A(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i, j) = 1.0 / (i + j + 1.0); // the index starts from 0!
        }
    }
    numalg::Matrix b(n);
    for (int i = 0; i < n; ++i) {
        double sum_i = 0.0;
        for (int j = 0; j < n; ++j) {
            sum_i += 1.0 / (i + j + 1.0); // the index starts from 0!
        }
        b(i) = sum_i;
    }

    numalg::Matrix x_cholesky = choleskySolve(A, b);
    outfile << "Prob2 Cholesky solution:\n";
    for (int i = 0; i < n; ++i) {
        outfile << x_cholesky(i) << " ";
    }
    outfile << "\n";

    numalg::Matrix x_modified_cholesky = modifiedCholeskySolve(A, b);
    outfile << "Prob2 Modified Cholesky solution:\n";
    for (int i = 0; i < n; ++i) {
        outfile << x_modified_cholesky(i) << " ";
    }
    outfile << "\n";

    numalg::Matrix x_gauss = gaussSolve(A, b);
    outfile << "Prob2 Gaussian elimination solution:\n";
    for (int i = 0; i < n; ++i) {
        outfile << x_gauss(i) << " ";
    }
    outfile << "\n";

    numalg::Matrix x_pGauss = PgaussSolve(A, b);
    outfile << "Prob2 Partial pivoting Gaussian elimination solution:\n";
    for (int i = 0; i < n; ++i) {
        outfile << x_pGauss(i) << " ";
    }
    outfile << "\n";
}

int main() {
    std::filesystem::create_directories("./output");
    std::ofstream outfile("./output/homework03.txt");
    if (!outfile.is_open()) {
        std::cerr << "Failed to open output file: ./output/homework03.txt" << std::endl;
        return 1;
    }
    prob1(outfile);
    prob2(outfile);
    return 0;
}