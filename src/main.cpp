#include <iostream>

// ext
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// prl2
#include "parallel.h"

// LensSym
#include "lens.h"

int main() {
  Ray ray(Vec3(0, 0, 0), Vec3(0, 0, 1));
  Lens lens(0, 1, 1, 1);
  lens.z = 10;
  Hit res;
  if (lens.intersect(ray, res)) {
    std::cout << res.hitPos << std::endl;
  }
  return 0;
}