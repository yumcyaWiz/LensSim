#ifndef _FILM_H
#define _FILM_H

#include <cmath>
#include <fstream>
#include <iostream>

#include "core/vec3.h"

using namespace Prl2;

class Film {
 public:
  unsigned int width;
  unsigned int height;
  Real width_length;
  Real height_length;
  Real diagonal_length;

  Vec3* pixels;

  Film(unsigned int _width, unsigned int _height, Real _width_length = 0.036f,
       Real _height_length = 0.024f)
      : width(_width),
        height(_height),
        width_length(_width_length),
        height_length(_height_length) {
    pixels = new Vec3[width * height];
    diagonal_length =
        std::sqrt(width_length * width_length + height_length * height_length);
  }

  ~Film() { delete[] pixels; }

  Vec3 getPixel(unsigned int i, unsigned int j) const {
    return pixels[i + width * j];
  }

  void setPixel(unsigned int i, unsigned int j, const Vec3& c) {
    pixels[i + width * j] = c;
  }

  Vec2 computePosition(Real u, Real v) const {
    return Vec2(0.5f * width_length * u, 0.5f * height_length * v);
  }

  void writePPM(const std::string& filename) const {
    std::ofstream file(filename);
    file << "P3" << std::endl;
    file << width << " " << height << std::endl;
    file << "255" << std::endl;

    for (int j = 0; j < height; ++j) {
      for (int i = 0; i < width; ++i) {
        const Vec3 c = getPixel(i, j);
        unsigned int r = static_cast<unsigned int>(255 * c.x());
        unsigned int g = static_cast<unsigned int>(255 * c.y());
        unsigned int b = static_cast<unsigned int>(255 * c.z());
        file << r << " " << g << " " << b << std::endl;
      }
    }

    file.close();
  }
};

#endif