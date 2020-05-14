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
  constexpr int width = 512;
  constexpr int height = 512;
  constexpr int num_samples = 100;
  const std::string path_to_lens = "../data/dgauss.50mm.json";

  std::shared_ptr<Sampler> sampler = std::make_shared<RandomSampler>();

  std::shared_ptr<Film> film = std::make_shared<Film>(width, height);
  LensSystem lsys(path_to_lens, film);

  IBL ibl("../data/PaperMill_E_3k.hdr");

  lsys.focus(-1);
  lsys.computeExitPupilBounds();

  Parallel parallel;
  parallel.parallelFor2D(
      [&](unsigned int i, unsigned int j) {
        const Real u =
            (2.0f * (i + sampler->getNext()) - film->width) / film->height;
        const Real v =
            (2.0f * (j + sampler->getNext()) - film->height) / film->height;

        Vec3 col;
        for (int k = 0; k < num_samples; ++k) {
          // sample ray
          Ray ray;
          if (!lsys.sampleRay(u, v, *sampler, ray)) continue;
          col += ibl.getRadiance(ray);
        }

        // divide by samples
        col /= num_samples;

        // set color on pixel
        film->setPixel(i, j, col);
      },
      32, 32, width, height);

  // output ppm
  film->gammaCorrection();
  film->writePPM("output.ppm");

  return 0;
}