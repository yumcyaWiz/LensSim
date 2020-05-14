#include <iostream>

// ext
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// prl2
#include "samplers/random.h"

// LensSym
#include "lens-system.h"

int main() {
  std::shared_ptr<Sampler> sampler = std::make_shared<RandomSampler>();

  std::shared_ptr<Film> film = std::make_shared<Film>(512, 512);
  LensSystem lsys("../data/dgauss.50mm.json", film);

  lsys.computeExitPupilBounds();
  Ray ray_out;
  for (int i = 0; i < 10; ++i) {
    Real u = (i + 0.5) / 10.0f;
    if (lsys.sampleRay(u, 0, *sampler, ray_out)) {
      std::cout << ray_out << std::endl;
    }
  }

  return 0;
}