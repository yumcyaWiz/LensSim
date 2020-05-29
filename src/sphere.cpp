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
  constexpr int num_samples = 1;
  const std::string path_to_lens = "../data/wide.22mm.json";

  constexpr Real plane_z = -1.0f;
  constexpr Real line_width = 0.1;
  constexpr Real grid_interval = 0.2;

  std::shared_ptr<Sampler> sampler = std::make_shared<RandomSampler>();

  std::shared_ptr<Film> film =
      std::make_shared<Film>(width, height, 0.024, 0.024);
  LensSystem lsys(path_to_lens, film);
  std::cout << lsys.image_focal_length << std::endl;

  lsys.focus(plane_z);
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
          const Real cos = std::abs(dot(ray.direction, Vec3(0, 0, -1)));

          // Compute Intersection with Plane
          const Real t = -(ray.origin.z() - plane_z) / ray.direction.z();
          const Vec3 hitPos = ray(t);

          // Grid Color
          const Real sx = std::sin(2 * PI / grid_interval * hitPos.x());
          const Real sy = std::sin(2 * PI / grid_interval * hitPos.y());
          Real intensity = 0;
          if (std::abs(sx) < line_width || std::abs(sy) < line_width) {
            intensity = 1.0f;
          }

          film->addPixel(i, j, intensity);
        }
      },
      32, 32, width, height);

  // output ppm
  film->divide(num_samples);
  film->writePPM("output.ppm");
  film->writeCSV("output.csv");

  return 0;
}
