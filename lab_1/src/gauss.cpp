#include "gauss.hpp"
#include <cmath>
#include <algorithm>

namespace slae {

Vector gaussianEliminationNoPivot(Matrix A, Vector b) {
    int n = A.size();
    for (int k = 0; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(A[k][k]) < 1e-12) continue;
            double f = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) A[i][j] -= f * A[k][j];
            b[i] -= f * b[k];
        }
    }
    Vector x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += A[i][j] * x[j];
        x[i] = (std::abs(A[i][i]) < 1e-12) ? 0 : (b[i] - sum) / A[i][i];
    }
    return x;
}

Vector gaussianEliminationPivot(Matrix A, Vector b) {
    int n = A.size();
    for (int k = 0; k < n - 1; ++k) {
        int maxRow = k;
        double maxVal = std::abs(A[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(A[i][k]) > maxVal) {
                maxVal = std::abs(A[i][k]);
                maxRow = i;
            }
        }
        if (maxRow != k) {
            std::swap(A[k], A[maxRow]);
            std::swap(b[k], b[maxRow]);
        }
        if (std::abs(A[k][k]) < 1e-12) continue;
        for (int i = k + 1; i < n; ++i) {
            double f = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) A[i][j] -= f * A[k][j];
            b[i] -= f * b[k];
        }
    }
    Vector x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += A[i][j] * x[j];
        x[i] = (std::abs(A[i][i]) < 1e-12) ? 0 : (b[i] - sum) / A[i][i];
    }
    return x;
}

} // namespace slae