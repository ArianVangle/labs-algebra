#ifndef SLAE_GAUSS_HPP
#define SLAE_GAUSS_HPP

#include "types.hpp"

namespace slae {
    Vector gaussianEliminationNoPivot(Matrix A, Vector b);
    Vector gaussianEliminationPivot(Matrix A, Vector b);
}

#endif