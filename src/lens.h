#ifndef _LENS_H
#define _LENS_H

#include <cmath>

#include "hit.h"
#include "ray.h"
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

  bool intersect(const Ray& ray, Hit& res) const {
    // aperture
    if (curvature_radius == 0) {
      Real t = -(ray.origin.z() - z) / ray.direction.z();
      Vec3 hitPos = ray(t);

      Real r = hitPos.x() * hitPos.x() + hitPos.y() * hitPos.y();
      if (r > aperture_radius * aperture_radius) return false;

      res.t = t;
      res.hitPos = ray(t);
      res.hitNormal = Vec3(0, 0, -1);
      return true;
    }
    // spherical lens
    else {
      Vec3 center(0, 0, z + curvature_radius);
      Real b = dot(ray.origin - center, ray.direction);
      Real c =
          length2(ray.origin - center) - curvature_radius * curvature_radius;
      Real D = b * b - c;
      if (D < 0) return false;

      Real t0 = -b - std::sqrt(D);
      Real t1 = -b + std::sqrt(D);
      Real t = curvature_radius * ray.direction.z() > 0 ? t0 : t1;
      Vec3 hitPos = ray(t);

      Real r = hitPos.x() * hitPos.x() + hitPos.y() * hitPos.y();
      if (r > aperture_radius * aperture_radius) return false;

      res.t = t;
      res.hitPos = ray(t);
      res.hitNormal = normalize(res.hitPos - center);
    }
  };
};

#endif