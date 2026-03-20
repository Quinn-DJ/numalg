// example 1.3.2 on page 32
#include <iostream>
#include <ctime>
#include <fstream>
#include <filesystem>

#include "matrix.hpp"
#include "operations.hpp"
#include "solve.hpp"
#include "gauss.hpp"

/**
 * @brief 随机生成 n 维列向量 b, 其元素是服从区间 [0, 1] 上的均匀分布的随机数
 * 
 * @return numalg::Matrix
 */
numalg::Matrix generate_vector_b(std::size_t n) {
    numalg::Matrix b(n);
    for (std::size_t i = 0; i < n; ++i) {
        b(i, 0) = static_cast<double>(rand()) / RAND_MAX;
    }
    return b;
}

/**
 * @brief 随机生成 n * n 矩阵 A, A = L L^T, 其中 L 是一个随机生成的下三角矩阵, L 的元素是服从区间 [1, 2] 上的均匀分布的随机数
 * 
 * @return numalg::Matrix
 */
numalg::Matrix generate_matrix_A(std::size_t n) {
    numalg::Matrix L(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i > j) {
                L(i, j) = 0;
            } else {
                L(i, j) = 1 + static_cast<double>(rand()) / RAND_MAX;
            }
        }
    }
    return L * L.transpose();
}

int main(int argc, char* argv[]) {
    if (!(argc == 1 || argc == 4)) {
        std::cerr << "Usage: " << argv[0] << " [START_N] [END_N] [STEP]" << std::endl;
        return 1;
    }
    // the output file is ./output/time_consuming.csv
    // n, time for gaussSolve, time for PgaussSolve, time for choleskySolve, time for modifiedCholeskySolve
    // Ensure output directory exists before opening the file.
    std::filesystem::create_directories("./output");
    std::ofstream outfile("./output/time_consuming.csv");
    if (!outfile.is_open()) {
        std::cerr << "Failed to open output file: ./output/time_consuming.csv" << std::endl;
        return 1;
    }
    outfile << "n, gaussSolve, PgaussSolve, choleskySolve, modifiedCholeskySolve\n";

    // argv[0] is the program name
    // argv[1] is the starting n, argv[2] is the ending n, argv[3] is the step
    int start_n = argv[1] ? std::stoi(argv[1]) : 10;
    int end_n = argv[2] ? std::stoi(argv[2]) : 500;
    int step = argv[3] ? std::stoi(argv[3]) : 10;
    // 设置随机数种子
    srand(time(0));
    for (int n = start_n; n <= end_n; n += step) {
        std::cout << "Testing n = " << n << "..." << std::endl;
        numalg::Matrix A = generate_matrix_A(n);
        numalg::Matrix b = generate_vector_b(n);
        outfile << n;

        clock_t start = clock();
        gaussSolve(A, b);
        clock_t end = clock();
        double time_gauss = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        outfile << ", " << time_gauss;

        start = clock();
        PgaussSolve(A, b);
        end = clock();
        double time_pgauss = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        outfile << ", " << time_pgauss;

        start = clock();
        choleskySolve(A, b);
        end = clock();
        double time_cholesky = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        outfile << ", " << time_cholesky;

        start = clock();
        modifiedCholeskySolve(A, b);
        end = clock();
        double time_modified_cholesky = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        outfile << ", " << time_modified_cholesky;

        outfile << "\n";
    }
    return 0;
}