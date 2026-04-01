#include "matrix_ops.hpp"
#include "utils.hpp"
#include <cmath>

namespace slae {

Vector matVecMult(const Matrix& A, const Vector& x) {
    int n = A.size();
    Vector res(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            res[i] += A[i][j] * x[j];
    return res;
}

double computeResidual(const Matrix& A, const Vector& x, const Vector& b) {
    Vector Ax = matVecMult(A, x);
    double sum = 0.0;
    for (size_t i = 0; i < b.size(); ++i) {
        double d = Ax[i] - b[i];
        sum += d * d;
    }
    return std::sqrt(sum);
}

double computeRelativeError(const Vector& approx, const Vector& exact) {
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < approx.size(); ++i) {
        double d = approx[i] - exact[i];
        num += d * d;
        den += exact[i] * exact[i];
    }
    return std::sqrt(num) / std::sqrt(den);
}

} // namespace slae