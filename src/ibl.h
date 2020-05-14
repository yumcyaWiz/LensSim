#ifndef _IBL_H
#define _IBL_H

#include <string>

// stb
#include "stb_image.h"

// prl2
#include "core/ray.h"
#include "core/vec3.h"

using namespace Prl2;

class IBL {
 public:
  int width;
  int height;

  Real* pixels;

  IBL(const std::string& filename) {
    int c;
    pixels = stbi_loadf(filename.c_str(), &width, &height, &c, 3);
  }
  ~IBL() { stbi_image_free(pixels); }

  Vec3 getRadiance(const Ray& ray) const {
    Real theta, phi;
    cartesianToSpherical(ray.direction, theta, phi);

    const Real u = phi * INV_PI_MUL_2;
    const Real v = theta * INV_PI;

    const int i = u * width;
    const int j = v * height;

    const Real r = pixels[3 * i + 3 * width * j];
    const Real g = pixels[3 * i + 3 * width * j + 1];
    const Real b = pixels[3 * i + 3 * width * j + 2];

    return Vec3(r, g, b);
  }
};

#endif