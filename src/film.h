#ifndef _FILM_H
#define _FILM_H

// prl2
#include "core/vec2.h"
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
       Real _height_length = 0.024f);

  ~Film();

  Vec3 getPixel(unsigned int i, unsigned int j) const;
  void setPixel(unsigned int i, unsigned int j, const Vec3& c);

  void addPixel(unsigned int i, unsigned int j, Real lambda, Real radiance);

  void divide(unsigned int k);

  Vec2 computePosition(Real u, Real v) const;

  void writePPM(const std::string& filename) const;
  void writeEXR(const std::string& filename) const;
};

#endif