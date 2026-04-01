#include "experiments.hpp"
#include "generators.hpp"
#include "gauss.hpp"
#include "lu.hpp"
#include "matrix_ops.hpp"
#include "utils.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <string>
#include <vector>

using namespace std::chrono;

namespace slae {

void experiment4_1() {
    std::cout << "\n=== Experiment 4.1: Single system timing ===\n";
    std::vector<int> sizes = {100, 200, 500, 1000};
    std::vector<std::string> rows;
    std::string header = "n,Gauss_No_Pivot_ms,Gauss_Pivot_ms,LU_Total_ms,LU_Decomp_ms,LU_Solve_ms";

    for (int n : sizes) {
        Matrix A = generateRandomMatrix(n, 42);
        Vector b = generateRandomVector(n, 42);

        auto t1 = high_resolution_clock::now();
        gaussianEliminationNoPivot(copyMatrix(A), copyVector(b));
        auto t2 = high_resolution_clock::now();
        double t_gauss_no = duration<double, std::milli>(t2 - t1).count();

        auto t3 = high_resolution_clock::now();
        gaussianEliminationPivot(copyMatrix(A), copyVector(b));
        auto t4 = high_resolution_clock::now();
        double t_gauss_pivot = duration<double, std::milli>(t4 - t3).count();

        auto t5 = high_resolution_clock::now();
        auto lu = luDecomposition(A);
        auto t6 = high_resolution_clock::now();
        double t_lu_decomp = duration<double, std::milli>(t6 - t5).count();

        auto t7 = high_resolution_clock::now();
        solveLU(lu.first, lu.second, b);
        auto t8 = high_resolution_clock::now();
        double t_lu_solve = duration<double, std::milli>(t8 - t7).count();

        double t_lu_total = t_lu_decomp + t_lu_solve;

        std::cout << "n=" << n << ": Gauss(no)=" << std::fixed << std::setprecision(4) << t_gauss_no
                  << ", Gauss(pivot)=" << t_gauss_pivot << ", LU(total)=" << t_lu_total << "\n";

        rows.push_back(std::to_string(n) + "," + std::to_string(t_gauss_no) + "," +
                       std::to_string(t_gauss_pivot) + "," + std::to_string(t_lu_total) + "," +
                       std::to_string(t_lu_decomp) + "," + std::to_string(t_lu_solve));
    }
    saveCSV("exp4_1.csv", header, rows);
}

void experiment4_2() {
    std::cout << "\n=== Experiment 4.2: Multiple RHS (n=500) ===\n";
    int n = 500;
    std::vector<int> k_vals = {1, 10, 100};
    Matrix A = generateRandomMatrix(n, 42);
    auto lu = luDecomposition(A);

    std::vector<std::string> rows;
    std::string header = "k,Gauss_Total_ms,LU_Total_ms";

    for (int k : k_vals) {
        double t_gauss = 0;
        for (int i = 0; i < k; ++i) {
            Vector b = generateRandomVector(n, 42 + i);
            auto t1 = high_resolution_clock::now();
            gaussianEliminationPivot(copyMatrix(A), b);
            auto t2 = high_resolution_clock::now();
            t_gauss += duration<double, std::milli>(t2 - t1).count();
        }

        auto t_start = high_resolution_clock::now();
        auto lu_local = luDecomposition(A);
        auto t_mid = high_resolution_clock::now();
        double t_decomp = duration<double, std::milli>(t_mid - t_start).count();

        double t_solve_sum = 0;
        for (int i = 0; i < k; ++i) {
            Vector b = generateRandomVector(n, 42 + i);
            auto t1 = high_resolution_clock::now();
            solveLU(lu_local.first, lu_local.second, b);
            auto t2 = high_resolution_clock::now();
            t_solve_sum += duration<double, std::milli>(t2 - t1).count();
        }
        double t_lu_total = t_decomp + t_solve_sum;

        std::cout << "k=" << k << ": Gauss=" << std::fixed << std::setprecision(4) << t_gauss
                  << ", LU=" << t_lu_total << "\n";
        rows.push_back(std::to_string(k) + "," + std::to_string(t_gauss) + "," + std::to_string(t_lu_total));
    }
    saveCSV("exp4_2.csv", header, rows);
}

void experiment4_3() {
    std::cout << "\n=== Experiment 4.3: Hilbert matrix accuracy ===\n";
    std::vector<int> sizes = {5, 10, 15};
    std::vector<std::string> rows;
    std::string header = "n,Method,Relative_Error,Residual";

    for (int n : sizes) {
        Matrix H = generateHilbertMatrix(n);
        Vector x_exact(n, 1.0);
        Vector b = matVecMult(H, x_exact);

        Vector x1 = gaussianEliminationNoPivot(copyMatrix(H), copyVector(b));
        double e1 = computeRelativeError(x1, x_exact);
        double r1 = computeResidual(H, x1, b);
        std::cout << "n=" << n << " NoPivot: Err=" << std::scientific << e1 << ", Res=" << r1 << "\n";
        rows.push_back(std::to_string(n) + ",Gauss_No_Pivot," + std::to_string(e1) + "," + std::to_string(r1));

        Vector x2 = gaussianEliminationPivot(copyMatrix(H), copyVector(b));
        double e2 = computeRelativeError(x2, x_exact);
        double r2 = computeResidual(H, x2, b);
        std::cout << "n=" << n << " Pivot:   Err=" << std::scientific << e2 << ", Res=" << r2 << "\n";
        rows.push_back(std::to_string(n) + ",Gauss_Pivot," + std::to_string(e2) + "," + std::to_string(r2));

        auto lu = luDecomposition(H);
        Vector x3 = solveLU(lu.first, lu.second, b);
        double e3 = computeRelativeError(x3, x_exact);
        double r3 = computeResidual(H, x3, b);
        std::cout << "n=" << n << " LU:      Err=" << std::scientific << e3 << ", Res=" << r3 << "\n";
        rows.push_back(std::to_string(n) + ",LU_Decomposition," + std::to_string(e3) + "," + std::to_string(r3));
    }
    saveCSV("exp4_3.csv", header, rows);
}

} // namespace slae