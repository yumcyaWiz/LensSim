#ifndef _SPHERE_H
#define _SPHERE_H

#include <cmath>

#include "core/hit.h"
#include "core/ray.h"
#include "core/type.h"
#include "core/vec3.h"

class Sphere {
 public:
  Vec3 center;
  Real radius;

  Sphere(const Vec3& _center, Real _radius)
      : center(_center), radius(_radius) {}

  bool intersect(const Ray& ray, Hit& res) const {
    const Real b = dot(ray.direction, ray.origin - center);
    const Real c = length2(ray.origin - center) - radius * radius;
    const Real D = b * b - c;
    if (D < 0) return false;

    const Real t0 = -b - std::sqrt(D);
    const Real t1 = -b + std::sqrt(D);
    Real t = t0;
    if (t < 0) {
      t = t1;
      if (t < 0) return false;
    }

    res.t = t;
    res.hitPos = ray(t);
    res.hitNormal = normalize(res.hitPos - center);
    return true;
  }
};

#endif