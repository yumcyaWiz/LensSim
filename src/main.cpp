#include <iostream>

// ext
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// prl2
#include "random.h"

// LensSym
#include "lens-system.h"

int main() {
  std::shared_ptr<Sampler> sampler = std::make_shared<RandomSampler>();

  std::shared_ptr<Film> film = std::make_shared<Film>(512, 512);
  LensSystem lsys("../data/dgauss.50mm.json", film);

  // Ray ray_in(Vec3(0, 0, 0), normalize(Vec3(0, 0, -1)));
  // Ray ray_out;
  // if (lsys.raytrace(ray_in, ray_out)) {
  //   std::cout << ray_out << std::endl;
  // }

  lsys.computeExitPupilBounds();
  Ray ray_out;
  for (int i = 0; i < 10; ++i) {
    if (lsys.sampleRay(Vec2(0, 0), *sampler, ray_out)) {
      std::cout << ray_out << std::endl;
    }
  }

  return 0;
}