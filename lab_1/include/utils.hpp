#ifndef SLAE_UTILS_HPP
#define SLAE_UTILS_HPP

#include "types.hpp"
#include <string>
#include <vector>

namespace slae {
    double norm(const Vector& v);
    Matrix copyMatrix(const Matrix& A);
    Vector copyVector(const Vector& v);
    void saveCSV(const std::string& filename, 
                 const std::string& header, 
                 const std::vector<std::string>& rows);
}

#endif