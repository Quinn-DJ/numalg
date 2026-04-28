#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>

#include "matrix.hpp"
#include "operations.hpp"
#include "solve.hpp"
#include "gauss.hpp"
#include "qr.hpp"

// Frobenius norm
double frobNorm(const numalg::Matrix& M) {
    double s = 0.0;
    for (std::size_t i = 0; i < M.lines(); ++i)
        for (std::size_t j = 0; j < M.rows(); ++j)
            s += M(i, j) * M(i, j);
    return std::sqrt(s);
}

// Max absolute difference between two vectors
double maxDiff(const numalg::Matrix& a, const numalg::Matrix& b) {
    double d = 0.0;
    for (std::size_t i = 0; i < a.lines(); ++i)
        for (std::size_t j = 0; j < a.rows(); ++j)
            d = std::max(d, std::abs(a(i, j) - b(i, j)));
    return d;
}

// ============================================
// Part 1: Compare QR solve with previous methods
// ============================================

void part1_prob_h02(std::ofstream& out) {
    // homework02: 84x84 tridiagonal system
    out << "========== Part 1: homework02 problem (84x84 tridiagonal) ==========\n";

    numalg::Matrix A(84, 84);
    for (std::size_t i = 0; i < 84; ++i) A(i, i) = 6;
    for (std::size_t i = 0; i < 83; ++i) {
        A(i, i + 1) = 1;
        A(i + 1, i) = 8;
    }

    numalg::Matrix b(84);
    for (std::size_t i = 0; i < 84; ++i) b(i) = 15;
    b(0) = 7; b(83) = 14;

    numalg::Matrix x_gauss = gaussSolve(A, b);
    numalg::Matrix x_pgauss = PgaussSolve(A, b);
    numalg::Matrix x_qr = qrSolve(A, b);

    out << std::scientific << std::setprecision(6);
    out << "Max diff (Gauss vs QR):      " << maxDiff(x_gauss, x_qr) << "\n";
    out << "Max diff (Pgauss vs QR):     " << maxDiff(x_pgauss, x_qr) << "\n";
    out << "Residual ||Ax - b|| (QR):    " << frobNorm(A * x_qr - b) << "\n\n";
}

void part1_prob1(std::ofstream& out) {
    // homework03 prob1: 100x100 tridiagonal system
    out << "========== Part 1: homework03 prob1 (100x100 tridiagonal) ==========\n";

    numalg::Matrix A(100, 100);
    for (int i = 0; i < 100; ++i) {
        A(i, i) = 10;
        if (i > 0) A(i, i - 1) = 1;
        if (i < 99) A(i, i + 1) = 1;
    }
    std::initializer_list<double> b_values = {
         34, -55, 46, 75, 4, 26, -58, 30, 39, 81, 74, 0, 48, 6, -38, 18, 27, 39, -27, -97, -31, 56, -75, 98, 72, 31, -59, -7, -71, 76, 84, -60, 42, 55, 84, -96, -96, 53, -67, 19, 0, 80, -31, 94, 84, 81, -40, -95, -7, 10, -16, -23, 64, -71, -84, -45, -40, 74, 42, 96, 100, -19, -72, 40, -13, -83, 89, 34, 28, 30, -57, 17, -58, -3, -95, -27, 88, -7, -68, -88, -18, -62, -21, -26, -38, -15, 78, 64, -4, 58, -27, 2, -46, 67, 57, -1, -59, 73, -77, 63
    };
    numalg::Matrix b = numalg::Matrix(100, b_values);

    numalg::Matrix x_cholesky = choleskySolve(A, b);
    numalg::Matrix x_gauss = gaussSolve(A, b);
    numalg::Matrix x_qr = qrSolve(A, b);

    out << std::scientific << std::setprecision(6);
    out << "Max diff (Cholesky vs QR):   " << maxDiff(x_cholesky, x_qr) << "\n";
    out << "Max diff (Gauss vs QR):      " << maxDiff(x_gauss, x_qr) << "\n";
    out << "Residual ||Ax - b|| (QR):    " << frobNorm(A * x_qr - b) << "\n\n";
}

void part1_prob2(std::ofstream& out) {
    // homework03 prob2: 20x20 Hilbert matrix system
    out << "========== Part 1: homework03 prob2 (20x20 Hilbert) ==========\n";

    int n = 20;
    numalg::Matrix A(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A(i, j) = 1.0 / (i + j + 1.0);

    numalg::Matrix b(n);
    for (int i = 0; i < n; ++i) {
        double sum_i = 0.0;
        for (int j = 0; j < n; ++j) sum_i += 1.0 / (i + j + 1.0);
        b(i) = sum_i;
    }

    numalg::Matrix x_cholesky = choleskySolve(A, b);
    numalg::Matrix x_gauss = gaussSolve(A, b);
    numalg::Matrix x_pgauss = PgaussSolve(A, b);
    numalg::Matrix x_qr = qrSolve(A, b);

    // Exact solution is [1, 1, ..., 1]
    numalg::Matrix x_exact(n, 1);
    for (int i = 0; i < n; ++i) x_exact(i, 0) = 1.0;

    out << std::scientific << std::setprecision(6);
    out << "Max diff (Cholesky vs exact):  " << maxDiff(x_cholesky, x_exact) << "\n";
    out << "Max diff (Gauss vs exact):     " << maxDiff(x_gauss, x_exact) << "\n";
    out << "Max diff (Pgauss vs exact):    " << maxDiff(x_pgauss, x_exact) << "\n";
    out << "Max diff (QR vs exact):        " << maxDiff(x_qr, x_exact) << "\n";
    out << "Residual ||Ax - b|| (QR):      " << frobNorm(A * x_qr - b) << "\n\n";
}

// ============================================
// Part 2: Quadratic polynomial least squares fit
// ============================================

void part2(std::ofstream& out) {
    out << "========== Part 2: Quadratic polynomial least squares ==========\n";

    double t[] = {-1.0, -0.75, -0.5, 0.0, 0.25, 0.5, 0.75};
    double y[] = {1.00, 0.8125, 0.75, 1.00, 1.3125, 1.75, 2.3125};
    int m = 7;

    // Build Vandermonde matrix: A = [1, t, t^2]
    numalg::Matrix A(m, 3);
    numalg::Matrix b(m, 1);
    for (int i = 0; i < m; ++i) {
        A(i, 0) = 1.0;
        A(i, 1) = t[i];
        A(i, 2) = t[i] * t[i];
        b(i, 0) = y[i];
    }

    numalg::Matrix x = qrLeastSquares(A, b);

    out << std::fixed << std::setprecision(6);
    out << "Fitting: y = ax^2 + bx + c\n";
    out << "a = " << x(2, 0) << "\n";
    out << "b = " << x(1, 0) << "\n";
    out << "c = " << x(0, 0) << "\n";

    // Compute residual ||Ax - b||
    double residual = frobNorm(A * x - b);
    out << std::scientific << std::setprecision(6);
    out << "Residual ||Ax - b||_2 = " << residual << "\n";

    // Print fitted values vs actual
    out << std::fixed << std::setprecision(4);
    out << "\nComparison:\n";
    out << "  t_i    y_i       fitted    residual\n";
    for (int i = 0; i < m; ++i) {
        double fitted = x(0, 0) + x(1, 0) * t[i] + x(2, 0) * t[i] * t[i];
        out << "  " << std::setw(6) << t[i] << "  " << std::setw(7) << y[i]
            << "  " << std::setw(9) << fitted
            << "  " << std::setw(9) << (fitted - y[i]) << "\n";
    }
    out << "\n";
}

// ============================================
// Part 3: Real estate linear model
// ============================================

void part3(std::ofstream& out) {
    out << "========== Part 3: Real estate least squares model ==========\n";

    double y_data[] = {
        25.9, 29.5, 27.9, 25.9, 29.9, 29.9, 30.9,
        28.9, 84.9, 82.9, 35.9, 31.5, 31.0, 30.9,
        30.0, 28.9, 36.9, 41.9, 40.5, 43.9, 37.5,
        37.9, 44.5, 37.9, 38.9, 36.9, 45.8, 41.0
    };
    int m = 28;

    // Model: y = x_0 + a_1*x_1 + a_2*x_2 + ... + a_11*x_11
    // Design matrix: column 0 = 1 (intercept), columns 1..11 = features
    numalg::Matrix b(m, 1);
    for (int i = 0; i < m; ++i) b(i, 0) = y_data[i];

    numalg::Matrix A(m, 12);
    for (int i = 0; i < m; ++i) A(i, 0) = 1.0;

    // Feature data: a_1=tax, a_2=bathrooms, a_3=lot_area, a_4=living_area,
    // a_5=garages, a_6=rooms, a_7=bedrooms, a_8=age,
    // a_9=building_type, a_10=floor_plan, a_11=fireplaces
    double features[28][11] = {
        { 4.9176,  1.0,  3.4720,  0.9980, 1.0,  7, 4, 42, 3, 1, 0},
        { 5.0208,  1.0,  3.5310,  1.5000, 2.0,  7, 4, 62, 1, 1, 0},
        { 4.5429,  1.0,  2.2750,  1.1750, 1.0,  6, 3, 40, 2, 1, 0},
        { 4.5573,  1.0,  4.0500,  1.2320, 1.0,  6, 3, 54, 4, 1, 0},
        { 5.0597,  1.0,  4.4550,  1.1210, 1.0,  6, 3, 42, 3, 1, 0},
        { 3.8910,  1.0,  4.4550,  0.9880, 1.0,  6, 3, 56, 2, 1, 0},
        { 5.8980,  1.0,  5.8500,  1.2400, 1.0,  7, 3, 51, 2, 1, 1},
        { 5.6039,  1.0,  9.5200,  1.5010, 0.0,  6, 3, 32, 1, 1, 0},
        {15.4202,  2.5,  9.8000,  3.4200, 2.0, 10, 5, 42, 2, 1, 1},
        {14.4598,  2.5, 19.8000,  3.7000, 2.0,  9, 5, 14, 4, 1, 1},
        { 5.8282,  1.0,  6.1350,  1.2250, 2.0,  6, 3, 32, 1, 1, 0},
        { 5.3003,  1.0,  4.9883,  1.5520, 1.0,  6, 3, 30, 1, 2, 0},
        { 6.2712,  1.0,  5.5200,  0.9750, 1.0,  5, 2, 30, 1, 2, 0},
        { 5.9592,  1.0,  6.6660,  1.1210, 2.0,  6, 3, 32, 2, 1, 0},
        { 5.0500,  1.0,  5.0000,  1.0200, 0.0,  5, 2, 46, 4, 1, 0},
        { 5.6039,  1.0,  9.5200,  1.5010, 0.0,  6, 3, 32, 1, 1, 0},
        { 8.2464,  1.5,  5.1500,  1.6640, 2.0,  8, 4, 50, 4, 1, 0},
        { 6.6969,  1.5,  6.0920,  1.4880, 1.5,  7, 3, 22, 1, 1, 1},
        { 7.7841,  1.5,  7.1020,  1.3760, 1.0,  6, 3, 17, 2, 1, 0},
        { 9.0384,  1.0,  7.8000,  1.5000, 1.5,  7, 3, 23, 3, 3, 0},
        { 5.9894,  1.0,  5.5200,  1.2560, 2.0,  6, 3, 40, 4, 1, 1},
        { 7.5422,  1.5,  4.0000,  1.6900, 1.0,  6, 3, 22, 1, 1, 0},
        { 8.7951,  1.5,  9.8900,  1.8200, 2.0,  8, 4, 50, 1, 1, 1},
        { 6.0931,  1.5,  6.7265,  1.6520, 1.0,  6, 3, 44, 4, 1, 0},
        { 8.3607,  1.5,  9.1500,  1.7770, 2.0,  8, 4, 48, 1, 1, 1},
        { 8.1400,  1.0,  8.0000,  1.5040, 2.0,  7, 3,  3, 1, 3, 0},
        { 9.1416,  1.5,  7.3262,  1.8310, 1.5,  8, 4, 31, 4, 1, 0},
        {12.0000,  1.5,  5.0000,  1.2000, 2.0,  6, 3, 30, 3, 1, 1},
    };
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < 11; ++j)
            A(i, j + 1) = features[i][j];

    numalg::Matrix x = qrLeastSquares(A, b);

    out << std::fixed << std::setprecision(6);
    out << "Model: y = x_0 + x_1*a_1 + x_2*a_2 + ... + x_11*a_11\n";
    out << "x_0 (intercept) = " << x(0, 0) << "\n";
    out << "x_1  (tax)        = " << x(1, 0) << "\n";
    out << "x_2  (bathrooms)  = " << x(2, 0) << "\n";
    out << "x_3  (lot_area)   = " << x(3, 0) << "\n";
    out << "x_4  (living_area)= " << x(4, 0) << "\n";
    out << "x_5  (garages)    = " << x(5, 0) << "\n";
    out << "x_6  (rooms)      = " << x(6, 0) << "\n";
    out << "x_7  (bedrooms)   = " << x(7, 0) << "\n";
    out << "x_8  (age)        = " << x(8, 0) << "\n";
    out << "x_9  (bldg_type)  = " << x(9, 0) << "\n";
    out << "x_10 (floor_plan) = " << x(10, 0) << "\n";
    out << "x_11 (fireplaces) = " << x(11, 0) << "\n";

    double residual = frobNorm(A * x - b);
    out << std::scientific << std::setprecision(6);
    out << "Residual ||Ax - b||_2 = " << residual << "\n\n";
}

int main() {
    std::filesystem::create_directories("./output");
    std::ofstream out("./output/homework08.txt");
    if (!out.is_open()) {
        std::cerr << "Failed to open output file: ./output/homework08.txt" << std::endl;
        return 1;
    }

    part1_prob_h02(out);
    part1_prob1(out);
    part1_prob2(out);
    part2(out);
    part3(out);

    std::cout << "Results written to ./output/homework08.txt" << std::endl;
    return 0;
}
