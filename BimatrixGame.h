#pragma once
#include <vector>  // std::vector
#include <utility> // std::pair

using BimatrixGame = std::vector<std::vector<std::pair<double, double>>>;

BimatrixGame createRandomGame(size_t size);
void printGame(BimatrixGame bg);
