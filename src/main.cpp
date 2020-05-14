#include <iostream>
#include <string>

// ext
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "ibl.h"
#include "lens-system/lens-system.h"
#include "parallel/parallel.h"
#include "samplers/random.h"

int main() {
  constexpr int width = 128;
  constexpr int height = 128;
  constexpr int num_samples = 100;
  const std::string path_to_lens = "../data/dgauss.50mm.json";

  std::shared_ptr<Sampler> sampler = std::make_shared<RandomSampler>();

  std::shared_ptr<Film> film = std::make_shared<Film>(width, height);
  LensSystem lsys(path_to_lens, film);

  lsys.computeExitPupilBounds();

  for (int j = 0; j < film->height; ++j) {
    for (int i = 0; i < film->width; ++i) {
      const Real u =
          (2.0f * (i + sampler->getNext()) - film->width) / film->height;
      const Real v =
          (2.0f * (j + sampler->getNext()) - film->height) / film->height;

      Vec3 col;
      for (int k = 0; k < num_samples; ++k) {
        // sample ray
        Ray ray;
        if (!lsys.sampleRay(u, v, *sampler, ray)) continue;
        col += 0.5f * (ray.direction + 1.0f);
      }

      // divide by samples
      col /= num_samples;

      // set color on pixel
      film->setPixel(i, j, col);
    }
  }

  // output ppm
  film->writePPM("output.ppm");

  return 0;
}