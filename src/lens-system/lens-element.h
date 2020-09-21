#ifndef _LENS_ELEMENT_H
#define _LENS_ELEMENT_H

#include <cmath>

#include "core/hit.h"
#include "core/ray.h"
#include "core/vec2.h"
#include "core/vec3.h"
#include "grid-data.h"
#include "lens-system/ior.h"

using namespace Prl2;

inline Vec3 alignNormal(const Vec3& v, const Vec3& n) {
  return dot(v, n) < 0 ? n : -n;
}

class LensElement {
 public:
  unsigned int index;
  Real curvature_radius;
  Real aperture_radius;
  Real thickness;

  Real z;

  std::shared_ptr<IOREquation> ior_equation;

  bool is_stop;

  LensElement(unsigned int _index, Real _aperture_radius, Real _thickness,
              Real _curvature_radius,
              const std::shared_ptr<IOREquation>& _ior_equation,
              bool _is_aperture)
      : index(_index),
        curvature_radius(_curvature_radius),
        aperture_radius(_aperture_radius),
        thickness(_thickness),
        z(0),
        ior_equation(_ior_equation),
        is_stop(_is_aperture) {}

  bool intersect(const Ray& ray, Hit& res) const {
    // if element is aperture or curvature radius is too big, treat element as
    // plane
    if (is_stop || curvature_radius > 10000) {
      Real t = -(ray.origin.z() - z) / ray.direction.z();
      Vec3 hitPos = ray(t);

      Real r = hitPos.x() * hitPos.x() + hitPos.y() * hitPos.y();
      if (r > aperture_radius * aperture_radius) return false;

      res.t = t;
      res.hitPos = ray(t);
      res.hitNormal = alignNormal(ray.direction, Vec3(0, 0, -1));
      return true;
    } else {
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
      res.hitPos = hitPos;
      res.hitNormal =
          alignNormal(ray.direction, normalize(res.hitPos - center));
      return true;
    }
  }

  Real ior(Real lambda) const { return ior_equation->ior(lambda); }

  GridData<Vec3> samplePoints(unsigned int N) const {
    GridData<Vec3> ret(N, N);

    // grid sampling
    for (int j = 0; j < N; ++j) {
      const Real v = (2.0 * (j + 0.5) - N) / N;
      const Real y = v * aperture_radius;
      for (int i = 0; i < N; ++i) {
        const Real u = (2.0 * (i + 0.5) - N) / N;
        const Real x = u * aperture_radius;
        ret.set(j, i, Vec3(x, y, z));
      }
    }

    return ret;
  }
};

#endif