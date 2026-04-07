#include "inverse_norm_estimate.hpp"
#include <vector>
#include <cmath>

namespace numalg {

// 向量无穷范数
static double vec_norm_inf(const std::vector<double>& v) {
    double max_val = 0.0;
    for (auto x : v) {
        double ax = std::abs(x);
        if (ax > max_val) max_val = ax;
    }
    return max_val;
}

// 向量1范数
static double vec_norm_1(const std::vector<double>& v) {
    double sum = 0.0;
    for (auto x : v) sum += std::abs(x);
    return sum;
}

// w = B * x
static std::vector<double> mat_vec_mul(const Matrix& B, const std::vector<double>& x) {
    std::size_t n = x.size();
    std::vector<double> w(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            w[i] += B(i, j) * x[j];
        }
    }
    return w;
}

// z = B^T * v
static std::vector<double> matT_vec_mul(const Matrix& B, const std::vector<double>& v) {
    std::size_t n = v.size();
    std::vector<double> z(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            z[i] += B(j, i) * v[j];
        }
    }
    return z;
}

// 点积 z^T * x
static double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

double estimate_inverse_norm_inf(const Matrix& B) {
    std::size_t n = B.lines();
    if (n != B.rows()) {
        throw std::invalid_argument("B must be a square matrix");
    }

    // x = e_1 (初始为单位向量)
    std::vector<double> x(n, 0.0);
    x[0] = 1.0;

    while (true) {
        // w = B * x
        std::vector<double> w = mat_vec_mul(B, x);

        // v = sign(w)
        std::vector<double> v(n);
        for (std::size_t i = 0; i < n; ++i) {
            v[i] = (w[i] >= 0.0) ? 1.0 : -1.0;
        }

        // z = B^T * v
        std::vector<double> z = matT_vec_mul(B, v);

        if (vec_norm_inf(z) <= dot(z, x)) {
            // 收敛：返回 ||w||_1
            return vec_norm_1(w);
        } else {
            // x = e_j，其中 |z_j| = ||z||_∞
            double z_inf = vec_norm_inf(z);
            std::fill(x.begin(), x.end(), 0.0);
            for (std::size_t j = 0; j < n; ++j) {
                if (std::abs(z[j]) == z_inf) {
                    x[j] = 1.0;
                    break;
                }
            }
        }
    }
}

}  // namespace numalg
