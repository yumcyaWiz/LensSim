#ifndef _LENS_ELEMENT_H
#define _LENS_ELEMENT_H

#include <cmath>

#include "core/hit.h"
#include "core/ray.h"
#include "core/vec2.h"
#include "core/vec3.h"

using namespace Prl2;

inline Vec3 alignNormal(const Vec3& v, const Vec3& n) {
  return dot(v, n) < 0 ? n : -n;
}

class LensElement {
 public:
  unsigned int index;
  Real aperture_radius;
  Real thickness;
  Real z;
  LensElement(unsigned int _index, Real _aperture_radius, Real _thickness)
      : index(_index),
        aperture_radius(_aperture_radius),
        thickness(_thickness),
        z(0){};

  virtual bool intersect(const Ray& ray, Hit& res) const = 0;
};

class Aperture : public LensElement {
 public:
  Aperture(unsigned int _index, Real _aperture_radius, Real _thickness)
      : LensElement(_index, _aperture_radius, _thickness){};

  bool intersect(const Ray& ray, Hit& res) const override {
    Real t = -(ray.origin.z() - z) / ray.direction.z();
    Vec3 hitPos = ray(t);

    Real r = hitPos.x() * hitPos.x() + hitPos.y() * hitPos.y();
    if (r > aperture_radius * aperture_radius) return false;

    res.t = t;
    res.hitPos = ray(t);
    res.hitNormal = alignNormal(ray.direction, Vec3(0, 0, -1));
    return true;
  };
};

class Lens : public LensElement {
 public:
  Real curvature_radius;
  Real ior;

  Lens(unsigned int _index, Real _aperture_radius, Real _thickness,
       Real _curvature_radius, Real _ior)
      : LensElement(_index, _aperture_radius, _thickness),
        curvature_radius(_curvature_radius),
        ior(_ior){};

  bool intersect(const Ray& ray, Hit& res) const override {
    Vec3 center(0, 0, z + curvature_radius);
    Real b = dot(ray.origin - center, ray.direction);
    Real c = length2(ray.origin - center) - curvature_radius * curvature_radius;
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
    res.hitNormal = alignNormal(ray.direction, normalize(res.hitPos - center));
    return true;
  }
};

#endif