#include "lu.hpp"
#include <cmath>

namespace slae {

std::pair<Matrix, Matrix> luDecomposition(const Matrix& A) {
    int n = A.size();
    Matrix L(n, Vector(n, 0.0));
    Matrix U(n, Vector(n, 0.0));
    
    for (int i = 0; i < n; ++i) {
        for (int k = i; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) sum += L[i][j] * U[j][k];
            U[i][k] = A[i][k] - sum;
        }
        L[i][i] = 1.0;
        for (int k = i + 1; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) sum += L[k][j] * U[j][i];
            L[k][i] = (std::abs(U[i][i]) < 1e-12) ? 0 : (A[k][i] - sum) / U[i][i];
        }
    }
    return {L, U};
}

Vector forwardSubstitution(const Matrix& L, const Vector& b) {
    int n = L.size();
    Vector y(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) sum += L[i][j] * y[j];
        y[i] = (b[i] - sum) / L[i][i];
    }
    return y;
}

Vector backwardSubstitution(const Matrix& U, const Vector& y) {
    int n = U.size();
    Vector x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += U[i][j] * x[j];
        x[i] = (std::abs(U[i][i]) < 1e-12) ? 0 : (y[i] - sum) / U[i][i];
    }
    return x;
}

Vector solveLU(const Matrix& L, const Matrix& U, const Vector& b) {
    return backwardSubstitution(U, forwardSubstitution(L, b));
}

} // namespace slae