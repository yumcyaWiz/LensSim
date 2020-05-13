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
  lsys.focus(-1);

  return 0;
}