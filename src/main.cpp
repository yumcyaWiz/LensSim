#include <iostream>

#include "parallel.h"

int main() {
  Prl2::Parallel parallel;
  parallel.parallelFor1D(
      [](unsigned int i) { std::cout << "asdf" << std::endl; }, 4, 4);
}