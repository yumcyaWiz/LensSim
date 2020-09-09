#include <cmath>
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
  const bool reflection = false;
  const std::string path_to_lens = "../data/dgauss50mm.json";

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
          const Real lambda_pdf = 1.0f / (SPD::LAMBDA_MAX - SPD::LAMBDA_MIN);
          const Real lambda =
              sampler->getNext() * (SPD::LAMBDA_MAX - SPD::LAMBDA_MIN) +
              SPD::LAMBDA_MIN;

          // sample ray
          Ray ray;
          Real ray_pdf;
          if (!lsys.sampleRay(u, v, lambda, *sampler, ray, ray_pdf, reflection))
            continue;
          const Real cos = std::abs(dot(ray.direction, Vec3(0, 0, -1)));

          // IBL
          const Real radiance =
              ibl.getRadiance(ray) * cos / (ray_pdf * lambda_pdf);

          film->addPixel(i, j, ray.lambda, radiance);
        }
      },
      32, 32, width, height);

  // output ppm
  film->divide(num_samples);
  film->writeEXR("output.exr");

  return 0;
}