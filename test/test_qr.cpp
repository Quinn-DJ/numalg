#include "matrix.hpp"
#include "qr.hpp"
#include "operations.hpp"
#include <iostream>
#include <cmath>
#include <cassert>

// Frobenius norm of a matrix
double frobNorm(const numalg::Matrix& M) {
    double s = 0.0;
    for (std::size_t i = 0; i < M.lines(); ++i)
        for (std::size_t j = 0; j < M.rows(); ++j)
            s += M(i, j) * M(i, j);
    return std::sqrt(s);
}

void test_house() {
    std::cout << "=== test_house ===" << std::endl;

    // Test: house(x) should give (I - beta*v*v^T)*x = ||x|| * e_1
    std::vector<double> x = {3.0, 4.0};
    auto [v, beta] = numalg::house(x);
    std::size_t n = x.size();

    // Apply H*x = x - beta*v*(v^T*x)
    double vTx = 0.0;
    for (std::size_t i = 0; i < n; ++i) vTx += v[i] * x[i];
    std::vector<double> Hx(n);
    for (std::size_t i = 0; i < n; ++i) Hx[i] = x[i] - beta * v[i] * vTx;

    double norm_x = std::sqrt(3.0 * 3.0 + 4.0 * 4.0);
    std::cout << "  ||Hx|| = " << std::sqrt(Hx[0]*Hx[0] + Hx[1]*Hx[1])
              << ", expected ||x|| = " << norm_x << std::endl;
    std::cout << "  Hx[1] = " << Hx[1] << " (should be ~0)" << std::endl;
    assert(std::abs(Hx[1]) < 1e-12);
    assert(std::abs(std::abs(Hx[0]) - norm_x) < 1e-12);  // |Hx[0]| = ||x||
    std::cout << "  PASSED" << std::endl;
}

void test_qr_decomposition() {
    std::cout << "=== test_qr_decomposition ===" << std::endl;

    // Square matrix
    numalg::Matrix A(3, 3);
    A(0,0)=12; A(0,1)=-51; A(0,2)=4;
    A(1,0)=6;  A(1,1)=167; A(1,2)=-68;
    A(2,0)=-4; A(2,1)=24;  A(2,2)=-41;

    auto [Q, R] = numalg::qrDecomposition(A);

    // Check Q*R ≈ A
    numalg::Matrix QR = Q * R;
    double err = frobNorm(QR - A);
    std::cout << "  ||QR - A|| = " << err << std::endl;
    assert(err < 1e-10);

    // Check Q is orthogonal: Q^T*Q ≈ I
    numalg::Matrix QtQ = numalg::tr(Q) * Q;
    numalg::Matrix I3(3, 3);
    for (std::size_t i = 0; i < 3; ++i) I3(i, i) = 1.0;
    double orth_err = frobNorm(QtQ - I3);
    std::cout << "  ||Q^T Q - I|| = " << orth_err << std::endl;
    assert(orth_err < 1e-10);

    // Check R is upper triangular
    for (std::size_t i = 1; i < 3; ++i)
        for (std::size_t j = 0; j < i; ++j)
            assert(std::abs(R(i, j)) < 1e-10);

    std::cout << "  PASSED" << std::endl;
}

void test_qr_rectangular() {
    std::cout << "=== test_qr_rectangular ===" << std::endl;

    // Tall matrix (4x2)
    numalg::Matrix A(4, 2);
    A(0,0)=1; A(0,1)=1;
    A(1,0)=1; A(1,1)=2;
    A(2,0)=1; A(2,1)=3;
    A(3,0)=1; A(3,1)=4;

    auto [Q, R] = numalg::qrDecomposition(A);

    // Check Q*R ≈ A
    double err = frobNorm(Q * R - A);
    std::cout << "  ||QR - A|| = " << err << std::endl;
    assert(err < 1e-10);

    // Q should be 4x4 orthogonal
    assert(Q.lines() == 4 && Q.rows() == 4);
    numalg::Matrix I4(4, 4);
    for (std::size_t i = 0; i < 4; ++i) I4(i, i) = 1.0;
    double orth_err = frobNorm(numalg::tr(Q) * Q - I4);
    std::cout << "  ||Q^T Q - I|| = " << orth_err << std::endl;
    assert(orth_err < 1e-10);

    std::cout << "  PASSED" << std::endl;
}

void test_qr_solve() {
    std::cout << "=== test_qr_solve ===" << std::endl;

    // Solve Ax = b
    numalg::Matrix A(3, 3);
    A(0,0)=2;  A(0,1)=1;  A(0,2)=-1;
    A(1,0)=-3; A(1,1)=-1; A(1,2)=2;
    A(2,0)=-2; A(2,1)=1;  A(2,2)=2;

    numalg::Matrix b(3, 1);
    b(0,0) = 8;
    b(1,0) = -11;
    b(2,0) = -3;

    numalg::Matrix x = numalg::qrSolve(A, b);

    // Expected: x = [2, 3, -1]
    std::cout << "  x = [" << x(0,0) << ", " << x(1,0) << ", " << x(2,0) << "]" << std::endl;
    assert(std::abs(x(0,0) - 2.0) < 1e-10);
    assert(std::abs(x(1,0) - 3.0) < 1e-10);
    assert(std::abs(x(2,0) + 1.0) < 1e-10);

    // Verify A*x = b
    numalg::Matrix Ax = A * x;
    double err = frobNorm(Ax - b);
    std::cout << "  ||Ax - b|| = " << err << std::endl;
    assert(err < 1e-10);

    std::cout << "  PASSED" << std::endl;
}

void test_qr_least_squares() {
    std::cout << "=== test_qr_least_squares ===" << std::endl;

    // Overdetermined system: fit y = c0 + c1*x to points (1,2), (2,3), (3,5)
    numalg::Matrix A(3, 2);
    A(0,0)=1; A(0,1)=1;  // x=1
    A(1,0)=1; A(1,1)=2;  // x=2
    A(2,0)=1; A(2,1)=3;  // x=3

    numalg::Matrix b(3, 1);
    b(0,0)=2;
    b(1,0)=3;
    b(2,0)=5;

    numalg::Matrix x = numalg::qrLeastSquares(A, b);

    std::cout << "  coeffs = [" << x(0,0) << ", " << x(1,0) << "]" << std::endl;
    // Expected: ~[0.333, 1.5]
    assert(std::abs(x(0,0) - 1.0/3.0) < 1e-10);
    assert(std::abs(x(1,0) - 1.5) < 1e-10);

    std::cout << "  PASSED" << std::endl;
}

void test_random_system() {
    std::cout << "=== test_random_system ===" << std::endl;

    // 10x10 diagonally dominant system for numerical stability
    std::srand(42);
    std::size_t n = 10;
    numalg::Matrix A(n, n);
    numalg::Matrix b(n, 1);
    numalg::Matrix x_true(n, 1);

    for (std::size_t i = 0; i < n; ++i) {
        x_true(i, 0) = (double)i;
        for (std::size_t j = 0; j < n; ++j) {
            A(i, j) = (std::rand() % 100) / 100.0;
        }
        A(i, i) += n;  // diagonally dominant
        b(i, 0) = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            b(i, 0) += A(i, j) * x_true(j, 0);
        }
    }

    numalg::Matrix I_n(n, n);
    for (std::size_t i = 0; i < n; ++i) I_n(i, i) = 1.0;
    auto [Q, R] = numalg::qrDecomposition(A);
    std::cout << "  ||QR - A|| = " << frobNorm(Q * R - A) << std::endl;
    std::cout << "  ||Q^T Q - I|| = " << frobNorm(numalg::tr(Q) * Q - I_n) << std::endl;
    numalg::Matrix x = numalg::qrSolve(A, b);
    double residual = frobNorm(A * x - b);
    std::cout << "  ||Ax - b|| = " << residual << std::endl;
    assert(residual < 1e-9);
    std::cout << "  PASSED" << std::endl;
}

int main() {
    test_house();
    test_qr_decomposition();
    test_qr_rectangular();
    test_qr_solve();
    test_qr_least_squares();
    test_random_system();

    std::cout << "\n✅ All tests passed!" << std::endl;
    return 0;
}
