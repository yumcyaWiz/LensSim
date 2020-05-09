#ifndef _LENS_H
#define _LENS_H

#include "vec3.h"
using namespace Prl2;

class Lens {
 public:
  unsigned int index;
  Real curvature_radius;
  Real aperture_radius;
  Real thickness;
  Real ior;
  Real z;

  Lens(unsigned int _index, Real _curvature_radius, Real _thickness, Real _ior)
      : index(_index),
        aperture_radius(_curvature_radius),
        thickness(_thickness),
        ior(_ior),
        z(0){};
};

#endif