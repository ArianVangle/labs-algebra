#ifndef SLAE_GENERATORS_HPP
#define SLAE_GENERATORS_HPP

#include "types.hpp"

namespace slae {
    Matrix generateRandomMatrix(int n, unsigned int seed);
    Vector generateRandomVector(int n, unsigned int seed);
    Matrix generateHilbertMatrix(int n);
}

#endif