#ifndef _LENS_SYSTEM_H
#define _LENS_SYSTEM_H

#include <cmath>
#include <vector>

#include "lens-element.h"

using namespace Prl2;

// 反射ベクトルを返す
inline Vec3 reflect(const Vec3& v, const Vec3& n) {
  return -v + 2 * dot(v, n) * n;
}

// フレネル係数を計算する
inline Real fresnel(const Vec3& wo, const Vec3& n, const Real& n1,
                    const Real& n2) {
  const float f0 = std::pow((n1 - n2) / (n1 + n2), 2.0f);
  return f0 + (1.0f - f0) * std::pow(1.0f - dot(wo, n), 5.0f);
}

// 屈折ベクトルを返す
inline bool refract(const Vec3& wi, Vec3& wt, const Vec3& n, const Real& ior1,
                    const Real& ior2) {
  const Real eta = ior1 / ior2;
  const Real cosThetaI = dot(wi, n);
  const Real sin2ThetaI = std::max(0.0f, 1.0f - cosThetaI * cosThetaI);
  const Real sin2ThetaT = eta * eta * sin2ThetaI;
  if (sin2ThetaT >= 1.0f) return false;
  const Real cosThetaT = std::sqrt(1.0f - sin2ThetaT);
  wt = eta * (-wi) + (eta * cosThetaI - cosThetaT) * n;
  return true;
}

class LensSystem {
 public:
  std::vector<LensElement> elements;

  // TODO: load lens-system from json
  LensSystem(const std::string& filename){};

  bool raytrace_from_object(const Ray& ray_in, Ray& ray_out,
                            bool reflection = false) const {
    int element_index = -1;
    Ray ray = ray_in;
    Real ior = 1.0f;

    while (true) {
      element_index += ray.direction.z() > 0 ? 1 : -1;
      if (element_index < 0 || element_index >= elements.size()) break;
      const LensElement* element = &elements[element_index];

      if (const Aperture* aperture = dynamic_cast<const Aperture*>(element)) {
        Hit res;
        if (aperture->intersect(ray, res)) {
          ray = Ray(res.hitPos, ray.direction);
          ior = 1.0f;
        } else {
          return false;
        }
      } else if (const Lens* lens = dynamic_cast<const Lens*>(element)) {
        Hit res;
        if (lens->intersect(ray, res)) {
          if (reflection) {
            // TODO: implement this
          } else {
            Vec3 next_direction;
            if (!refract(-ray.direction, next_direction, res.hitNormal, ior,
                         next_ior))
              return false;
          }
        } else {
          return false;
        }
      }
    }
  };
};

#endif