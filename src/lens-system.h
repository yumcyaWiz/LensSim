#ifndef _LENS_SYSTEM_H
#define _LENS_SYSTEM_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

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

  bool loadJSON(const std::string& filename) {
    // open file
    std::ifstream stream(filename);
    if (!stream) {
      std::cerr << "failed to open " << filename << std::endl;
      return false;
    }

    // parse JSON
    json j;
    stream >> j;

    std::cout << "a" << std::endl;

    return true;
  };

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
        if (!aperture->intersect(ray, res)) return false;
        ray = Ray(res.hitPos, ray.direction);
        ior = 1.0f;
      } else if (const Lens* lens = dynamic_cast<const Lens*>(element)) {
        // Compute Next Element
        const int next_element_index =
            ray.direction.z() > 0 ? element_index : element_index - 1;
        const LensElement* next_element = &elements[element_index + 1];

        // Compute Next Element IOR
        Real next_ior = 0;
        if (const Aperture* next_aperture =
                dynamic_cast<const Aperture*>(next_element)) {
          next_ior = 1.0f;
        } else if (const Lens* next_lens =
                       dynamic_cast<const Lens*>(next_element)) {
          next_ior = next_lens->ior;
        }

        // Compute Intersection with Lens
        Hit res;
        if (!lens->intersect(ray, res)) return false;

        // Refract and Reflect
        if (reflection) {
          // TODO: implement this
        } else {
          Vec3 next_direction;
          if (!refract(-ray.direction, next_direction, res.hitNormal, ior,
                       next_ior))
            return false;

          // Set Next Ray
          ray = Ray(res.hitPos, next_direction);
        }
      }
    }

    ray_out = ray;

    return true;
  };
};

#endif