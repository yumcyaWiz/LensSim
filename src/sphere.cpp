#include <cmath>
#include <iostream>
#include <string>

// ext
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "core/spectrum.h"
#include "lens-system/lens-system.h"
#include "linear.h"
#include "parallel/parallel.h"
#include "samplers/random.h"
#include "sphere.h"

int main() {
  constexpr int width = 1024;
  constexpr int height = 1024;
  constexpr int num_samples = 300;
  const std::string path_to_lens = "../data/wide.22mm.json";

  std::shared_ptr<Sampler> sampler = std::make_shared<RandomSampler>();

  std::shared_ptr<Film> film =
      std::make_shared<Film>(width, height, 0.024, 0.024);
  LensSystem lsys(path_to_lens, film);

  std::shared_ptr<Film> film_depth =
      std::make_shared<Film>(width, height, 0.024, 0.024);

  // Scene
  LinearIntersector intersector;
  constexpr Real radius = 0.2f;
  intersector.add(std::make_shared<Sphere>(Vec3(0, 0, -1), radius));
  intersector.add(std::make_shared<Sphere>(Vec3(0.5, 0, -2), radius));
  intersector.add(std::make_shared<Sphere>(Vec3(1, 0, -3), radius));

  lsys.focus(-2);
  lsys.computeExitPupilBounds();

  Parallel parallel;
  parallel.parallelFor2D(
      [&](unsigned int i, unsigned int j) {
        for (int k = 0; k < num_samples; ++k) {
          const Real u =
              (2.0f * (i + sampler->getNext()) - film->width) / film->width;
          const Real v =
              (2.0f * (j + sampler->getNext()) - film->height) / film->height;

          // sample ray
          Ray ray;
          Real ray_pdf;
          Real lambda = 550.0f;
          if (!lsys.sampleRay(u, v, lambda, *sampler, ray, ray_pdf, false))
            continue;

          Hit res;
          if (intersector.intersect(ray, res)) {
            // add depth
            film_depth->addPixel(i, j, Vec3(res.t));

            // compute spherical coordinate
            Real phi = std::atan2(res.hitNormal.z(), res.hitNormal.x());
            if (phi < 0) phi += PI_MUL_2;
            const Real theta =
                std::acos(std::clamp(res.hitNormal.y(), -1.0f, 1.0f));

            // checkerboard
            Vec3 color(0.1, 0.1, 0.1);
            constexpr Real sl = 0.1;
            if (std::sin(phi / sl) * std::sin(theta / sl) > 0) {
              color = Vec3(0.8, 0.8, 0.8);
            }

            // simple shading
            const Vec3 sunDir = normalize(Vec3(1, 1, 1));
            color *= std::max(dot(res.hitNormal, sunDir), 0.0f);

            // add color
            film->addPixel(i, j, color);
          } else {
            film->addPixel(i, j, Vec3(0));
            film_depth->addPixel(i, j, Vec3(0));
          }
        }
      },
      32, 32, width, height);

  // output ppm
  // film->writePPM("output.ppm");
  // film->writeEXR("output.exr");
  film->writeTIFF("output.tiff");
  film_depth->writeTIFF("depth.tiff");

  return 0;
}
