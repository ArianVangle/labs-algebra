#ifndef SLAE_LU_HPP
#define SLAE_LU_HPP

#include "types.hpp"
#include <utility>

namespace slae {
    std::pair<Matrix, Matrix> luDecomposition(const Matrix& A);
    Vector forwardSubstitution(const Matrix& L, const Vector& b);
    Vector backwardSubstitution(const Matrix& U, const Vector& y);
    Vector solveLU(const Matrix& L, const Matrix& U, const Vector& b);
}

#endif