#include <cmath>
#include <iostream>
#include <string>

// ext
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "core/spectrum.h"
#include "ibl.h"
#include "lens-system/lens-system.h"
#include "linear.h"
#include "parallel/parallel.h"
#include "samplers/random.h"
#include "sphere.h"

int main() {
  constexpr int width = 512;
  constexpr int height = 512;
  constexpr int num_samples = 100;
  const std::string path_to_lens = "../data/wide.22mm.json";

  std::shared_ptr<Sampler> sampler = std::make_shared<RandomSampler>();

  std::shared_ptr<Film> film =
      std::make_shared<Film>(width, height, 0.024, 0.024);
  LensSystem lsys(path_to_lens, film);

  // Scene
  LinearIntersector intersector;
  intersector.add(std::make_shared<Sphere>(Vec3(0, 0, -1), 0.2));
  intersector.add(std::make_shared<Sphere>(Vec3(0.5, 0, -2), 0.2));
  intersector.add(std::make_shared<Sphere>(Vec3(1, 0, -3), 0.2));

  lsys.focus(-1);
  lsys.computeExitPupilBounds();

  Parallel parallel;
  parallel.parallelFor2D(
      [&](unsigned int i, unsigned int j) {
        const Real u = (2.0f * (i + 0.5f) - film->width) / film->width;
        const Real v = (2.0f * (j + 0.5f) - film->height) / film->height;

        for (int k = 0; k < num_samples; ++k) {
          // sample ray
          Ray ray;
          Real ray_pdf;
          Real lambda = 550.0f;
          if (!lsys.sampleRay(u, v, lambda, *sampler, ray, ray_pdf, false))
            continue;

          Hit res;
          if (intersector.intersect(ray, res)) {
            film->addPixel(i, j, 0.5f * (res.hitNormal + 1.0f));
          } else {
            film->addPixel(i, j, Vec3(0));
          }
        }
      },
      32, 32, width, height);

  // output ppm
  film->writePPM("output.ppm");
  film->writeEXR("output.exr");

  return 0;
}
