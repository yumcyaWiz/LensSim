#include <iostream>

// ext
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// prl2
#include "parallel.h"

// LensSym
#include "lens-system.h"

int main() {
  LensSystem lsys;
  lsys.loadJSON("../data/dgauss.50mm.json");

  return 0;
}