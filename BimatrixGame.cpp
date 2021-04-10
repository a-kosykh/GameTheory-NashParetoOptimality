#include "BimatrixGame.h"
#include <ctime>
#include <iostream>
#include <iomanip>

BimatrixGame createRandomGame(size_t n)
{
    BimatrixGame rv;
    
    std::srand(std::time(nullptr));

    for (size_t i = 0; i < n; i++) {
        std::vector<std::pair<double, double>> row;
        for (size_t j = 0; j < n; ++j) {
            std::pair<double, double> pair = 
                std::make_pair<double,double>(
                    (std::rand() % 100) - 50, 
                    (std::rand() % 100) - 50);
            row.push_back(pair);
        }
        rv.push_back(row);
    }

    return rv;
}

void printGame(BimatrixGame bg)
{
    for (const auto &row : bg) {
        for (const auto& pair : row) {
            std::cout << "(" << std::setw(3) << pair.first << "/" <<
                std::setw(3) << pair.second << ")  ";
        }
        std::cout << std::endl;
    }
}
