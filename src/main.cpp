#include <iostream>

// ext
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// prl2
#include "parallel.h"

// LensSym
#include "lens-system.h"

int main() {
  std::shared_ptr<Film> film = std::make_shared<Film>(512, 512);
  LensSystem lsys("../data/dgauss.50mm.json", film);

  lsys.computeExitPupilBounds();

  return 0;
}