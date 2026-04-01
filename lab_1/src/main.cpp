#include <iostream>
#include "experiments.hpp"

int main() {
    std::cout << "Lab 1: SLAE Methods Comparison (C++)\n";
    std::cout << "=====================================\n";
    
    slae::experiment4_1();
    slae::experiment4_2();
    slae::experiment4_3();

    std::cout << "\n✓ All experiments completed. Check CSV files.\n";
    return 0;
}