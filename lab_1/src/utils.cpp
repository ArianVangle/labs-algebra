#include "utils.hpp"
#include <cmath>
#include <fstream>
#include <iostream>

namespace slae {

double norm(const Vector& v) {
    double sum = 0.0;
    for (double val : v) sum += val * val;
    return std::sqrt(sum);
}

Matrix copyMatrix(const Matrix& A) { return A; }
Vector copyVector(const Vector& v) { return v; }

void saveCSV(const std::string& filename, 
             const std::string& header, 
             const std::vector<std::string>& rows) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << "\n";
        return;
    }
    file << header << "\n";
    for (const auto& row : rows) file << row << "\n";
    file.close();
    std::cout << "[OK] Saved: " << filename << "\n";
}

} // namespace slae