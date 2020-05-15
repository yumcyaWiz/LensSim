#include <iostream>
#include <string>

// ext
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "core/spectrum.h"
#include "ibl.h"
#include "lens-system/lens-system.h"
#include "parallel/parallel.h"
#include "samplers/random.h"

int main() {
  constexpr int width = 512;
  constexpr int height = 512;
  constexpr int num_samples = 100;
  const std::string path_to_lens = "../data/wide.22mm.json";

  std::shared_ptr<Sampler> sampler = std::make_shared<RandomSampler>();

  std::shared_ptr<Film> film = std::make_shared<Film>(width, height);
  LensSystem lsys(path_to_lens, film);
  std::cout << lsys.image_focal_length << std::endl;

  IBL ibl("../data/PaperMill_E_3k.hdr");

  lsys.focus(-0.2);
  lsys.computeExitPupilBounds();

  Parallel parallel;
  parallel.parallelFor2D(
      [&](unsigned int i, unsigned int j) {
        const Real u =
            (2.0f * (i + sampler->getNext()) - film->width) / film->height;
        const Real v =
            (2.0f * (j + sampler->getNext()) - film->height) / film->height;

        for (int k = 0; k < num_samples; ++k) {
          // sample wavelength
          Real lambda_pdf = 1.0f / (SPD::LAMBDA_MAX - SPD::LAMBDA_MIN);
          Real lambda =
              sampler->getNext() * (SPD::LAMBDA_MAX - SPD::LAMBDA_MIN) +
              SPD::LAMBDA_MIN;
          lambda = std::clamp(lambda, SPD::LAMBDA_MIN, SPD::LAMBDA_MAX);

          // sample ray
          Ray ray;
          ray.lambda = lambda;
          if (!lsys.sampleRay(u, v, *sampler, ray)) continue;

          // IBL
          Real radiance = ibl.getRadiance(ray);

          film->addPixel(i, j, ray.lambda, radiance);
        }
      },
      32, 32, width, height);

  // output ppm
  film->divide(num_samples);
  film->gammaCorrection();
  film->writePPM("output.ppm");

  return 0;
}