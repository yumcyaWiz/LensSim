#include <iostream>

// ext
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// prl2
#include "parallel.h"

// LensSym
#include "lens-system.h"

int main() {
  // LensSystem lsys("../data/dgauss.50mm.json");
  // lsys.focus(-1);
  Film film(512, 512);
  for (int i = 0; i < film.width; ++i) {
    for (int j = 0; j < film.height; ++j) {
      film.setPixel(i, j,
                    Vec3(float(i) / film.width, float(j) / film.height, 1.0f));
    }
  }

  film.writePPM("output.ppm");

  return 0;
}