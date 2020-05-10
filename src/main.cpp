#include <iostream>

// ext
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// prl2
#include "parallel.h"

// LensSym
#include "lens-system.h"

int main() {
  LensSystem lsys("../data/dgauss.50mm.json");

  Ray ray_in(Vec3(0, 0, 0), Vec3(0, 0, 1));
  Ray ray_out;
  if (lsys.raytraceFromObject(ray_in, ray_out)) {
    std::cout << ray_out << std::endl;
  }

  return 0;
}