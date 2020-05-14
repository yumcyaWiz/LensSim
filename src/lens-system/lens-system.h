#ifndef _LENS_SYSTEM_H
#define _LENS_SYSTEM_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "nlohmann/json.hpp"
using JSON = nlohmann::json;

#include "core/bounds2.h"
#include "film.h"
#include "lens-system/lens-element.h"
#include "parallel/parallel.h"

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

// 点を平面で原点中心に回転する
inline Vec2 rotate2D(const Vec2& p, Real theta) {
  return Vec2(p.x() * std::cos(theta) - p.y() * std::sin(theta),
              p.x() * std::sin(theta) + p.y() * std::cos(theta));
}

class LensSystem {
 public:
  std::shared_ptr<Film> film;

  std::vector<std::shared_ptr<LensElement>> elements;

  Real object_focal_z;
  Real object_principal_z;
  Real object_focal_length;
  Real image_focal_z;
  Real image_principal_z;
  Real image_focal_length;

  static constexpr unsigned int num_exit_pupil_bounds = 64;
  static constexpr unsigned int num_exit_pupil_bounds_samples = 1024;
  std::vector<Bounds2> exit_pupil_bounds;

  LensSystem(const std::string& filename, const std::shared_ptr<Film> _film)
      : film(_film) {
    // load json
    if (!loadJSON(filename)) exit(EXIT_FAILURE);

    // compute system length and z
    Real length = 0;
    for (auto itr = elements.rbegin(); itr != elements.rend(); itr++) {
      length += (*itr)->thickness;
      (*itr)->z = -length;
    }

    // compute cardinal points
    if (!computeCardinalPoints()) exit(EXIT_FAILURE);

    // focus at z = -inf
    if (!focus(-10000)) {
      std::cerr << "failed to focus lens at z = -inf" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  bool loadJSON(const std::string& filename) {
    // open file
    std::ifstream stream(filename);
    if (!stream) {
      std::cerr << "failed to open " << filename << std::endl;
      return false;
    }

    // parse JSON
    JSON json;
    stream >> json;

    // push lens elements
    for (const auto& [key, value] : json.items()) {
      // unit conversion from [mm] to [m]
      const unsigned int index = value["index"].get<unsigned int>();
      const Real curvature_radius =
          value["curvature_radius"].get<Real>() * 1e-3f;
      const Real thickness = value["thickness"].get<Real>() * 1e-3f;
      const Real ior = value["eta"].get<Real>();
      const Real aperture_radius =
          0.5f * value["aperture_diameter"].get<Real>() * 1e-3f;

      // Aperture
      if (curvature_radius == 0.0f) {
        auto element =
            std::make_shared<Aperture>(index, aperture_radius, thickness);
        elements.push_back(element);
      }
      // Lens
      else {
        auto element = std::make_shared<Lens>(index, aperture_radius, thickness,
                                              curvature_radius, ior);
        elements.push_back(element);
      }
    }

    // sort lens elements by index
    std::sort(elements.begin(), elements.end(),
              [](const std::shared_ptr<LensElement>& x1,
                 const std::shared_ptr<LensElement>& x2) {
                return x1->index < x2->index;
              });

    return true;
  };

  bool raytrace(const Ray& ray_in, Ray& ray_out,
                bool reflection = false) const {
    int element_index = ray_in.direction.z() > 0 ? -1 : elements.size();
    Ray ray = ray_in;
    Real ior = 1.0f;

    while (true) {
      // update element index
      element_index += ray.direction.z() > 0 ? 1 : -1;
      if (element_index < 0 || element_index >= elements.size()) break;
      const auto element = elements[element_index];

      // Aperture
      if (const std::shared_ptr<Aperture> aperture =
              std::dynamic_pointer_cast<Aperture>(element)) {
        Hit res;
        if (!aperture->intersect(ray, res)) return false;

        // Update ray
        ray = Ray(res.hitPos, ray.direction);

        // Update ior
        ior = 1.0f;
      }
      // Lens
      else if (const std::shared_ptr<Lens> lens =
                   std::dynamic_pointer_cast<Lens>(element)) {
        // Compute Next Element IOR
        Real next_ior = 1.0f;

        // Compute Next Element
        const int next_element_index =
            ray.direction.z() > 0 ? element_index : element_index - 1;
        if (next_element_index >= 0) {
          const std::shared_ptr<LensElement> next_element =
              elements[next_element_index];
          if (const std::shared_ptr<Lens> next_lens =
                  std::dynamic_pointer_cast<Lens>(next_element)) {
            next_ior = next_lens->ior;
          }
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
          ray = Ray(res.hitPos, normalize(next_direction));

          // update ior
          ior = next_ior;
        }
      } else {
        std::cerr << "invalid lens element" << std::endl;
        return false;
      }
    }

    ray_out = ray;

    return true;
  }

  bool computeCardinalPoints() {
    const Real height = 0.001f;

    // raytrace from object plane
    Ray ray_in(Vec3(0, height, elements.front()->z - 1.0f), Vec3(0, 0, 1));
    Ray ray_out;
    if (!raytrace(ray_in, ray_out)) {
      std::cerr << "failed to compute cardinal points" << std::endl;
      return false;
    }

    // compute image focal point
    Real t = -ray_out.origin.y() / ray_out.direction.y();
    image_focal_z = ray_out(t).z();

    // compute image principal point
    t = -(ray_out.origin.y() - height) / ray_out.direction.y();
    image_principal_z = ray_out(t).z();

    // compute image focal length
    image_focal_length = image_focal_z - image_principal_z;

    // raytrace from image plane
    ray_in = Ray(Vec3(0, height, 0), Vec3(0, 0, -1));
    if (!raytrace(ray_in, ray_out)) {
      std::cerr << "failed to compute cardinal points" << std::endl;
      return false;
    }

    // compute object focal point
    t = -ray_out.origin.y() / ray_out.direction.y();
    object_focal_z = ray_out(t).z();

    // compute object principal point
    t = -(ray_out.origin.y() - height) / ray_out.direction.y();
    object_principal_z = ray_out(t).z();

    // compute object focal length
    object_focal_length = object_focal_z - object_principal_z;

    return true;
  }

  bool focus(Real focus_z) {
    const Real delta =
        0.5f * (object_principal_z - focus_z + image_principal_z -
                std::sqrt((object_principal_z - focus_z - image_principal_z) *
                          (object_principal_z - focus_z -
                           4 * image_focal_length - image_principal_z)));

    // move lens elements
    for (auto& element : elements) {
      element->z -= delta;
    }

    // recompute cardinal points
    if (!computeCardinalPoints()) return false;

    return true;
  }

  Bounds2 computeExitPupilBound(const Vec2& p) const {
    Bounds2 bounds;

    const auto lastElement = elements.back();
    Ray ray_out;
    for (int i = 0; i < num_exit_pupil_bounds_samples; ++i) {
      for (int j = 0; j < num_exit_pupil_bounds_samples; ++j) {
        // sample point on last element surface
        const Real u =
            2.0f * static_cast<Real>(i) / num_exit_pupil_bounds_samples - 1.0f;
        const Real v =
            2.0f * static_cast<Real>(j) / num_exit_pupil_bounds_samples - 1.0f;
        const Vec3 samplePoint =
            Vec3(lastElement->aperture_radius * u,
                 lastElement->aperture_radius * v, lastElement->z);

        // make ray
        const Vec3 origin(p.x(), p.y(), 0);
        const Ray ray_in(origin, normalize(samplePoint - origin));

        // raytrace
        if (!raytrace(ray_in, ray_out)) continue;

        // extend bounding box
        bounds = extendBounds(bounds, Vec2(samplePoint.x(), samplePoint.y()));
      }
    }

    return bounds;
  }

  bool computeExitPupilBounds() {
    exit_pupil_bounds.resize(num_exit_pupil_bounds);

    Parallel parallel;

    parallel.parallelFor1D(
        [&](unsigned int idx) {
          const Real r = static_cast<Real>(idx) / num_exit_pupil_bounds * 0.5f *
                         film->diagonal_length;
          exit_pupil_bounds[idx] = computeExitPupilBound(Vec2(r, 0));

          std::cout << "finished " << idx
                    << "th computation of exit pupil bounds" << std::endl;
          std::cout << exit_pupil_bounds[idx] << std::endl;
        },
        16, num_exit_pupil_bounds);

    return true;
  }

  bool sampleRay(Real u, Real v, Sampler& sampler, Ray& ray_out) const {
    // compute position on film
    const Vec2 p = film->computePosition(u, v);

    // choose exit pupil bound
    const Real r = length(p);
    const unsigned int exit_pupil_bounds_index =
        r / (0.5f * film->diagonal_length) * num_exit_pupil_bounds;
    const Bounds2& exit_pupil_bound =
        exit_pupil_bounds[exit_pupil_bounds_index];
    if (!exit_pupil_bound.isValid()) return false;

    // sample point on exit pupil bound
    Real pdf_area;
    Vec2 pBound = exit_pupil_bound.samplePoint(sampler, pdf_area);

    // rotate sampled point
    if (r > 0) {
      const Real theta = std::atan2(v, u);
      pBound = rotate2D(pBound, theta);
    }

    // make input ray
    const Vec3 origin = Vec3(p.x(), p.y(), 0);
    const Vec3 direction =
        normalize(Vec3(pBound.x(), pBound.y(), elements.back()->z) - origin);
    const Ray ray_in(origin, direction);

    // raytrace
    Ray ray_tmp;
    if (!raytrace(ray_in, ray_tmp)) return false;

    ray_out = ray_tmp;

    return true;
  }
};

#endif