#ifndef NUMALG_SOLVE_HPP
#define NUMALG_SOLVE_HPP

#include "matrix.hpp"

numalg::Matrix solve_lower_triangular(const numalg::Matrix& L, const numalg::Matrix& b);
numalg::Matrix solve_upper_triangular(const numalg::Matrix& U, const numalg::Matrix& b);

#endif  // NUMALG_SOLVE_HPP