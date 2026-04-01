#include "generators.hpp"
#include <random>

namespace slae {

Matrix generateRandomMatrix(int n, unsigned int seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    Matrix A(n, Vector(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = dist(gen);
    return A;
}

Vector generateRandomVector(int n, unsigned int seed) {
    std::mt19937 gen(seed + 1000);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    Vector b(n);
    for (int i = 0; i < n; ++i) b[i] = dist(gen);
    return b;
}

Matrix generateHilbertMatrix(int n) {
    Matrix H(n, Vector(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            H[i][j] = 1.0 / static_cast<double>(i + j + 1);
    return H;
}

} // namespace slae