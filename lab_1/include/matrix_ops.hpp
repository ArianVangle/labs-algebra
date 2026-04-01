#ifndef SLAE_MATRIX_OPS_HPP
#define SLAE_MATRIX_OPS_HPP

#include "types.hpp"

namespace slae {
    Vector matVecMult(const Matrix& A, const Vector& x);
    double computeResidual(const Matrix& A, const Vector& x, const Vector& b);
    double computeRelativeError(const Vector& approx, const Vector& exact);
}

#endif