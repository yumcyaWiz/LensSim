#ifndef _FILM_H
#define _FILM_H

#include <cmath>

#include "vec3.h"

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
};

#endif