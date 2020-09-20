#ifndef _LENS_SYSTEM_H
#define _LENS_SYSTEM_H

#include <array>
#include <tuple>
#include <vector>

#include "core/bounds2.h"
#include "film.h"
#include "grid-data.h"
#include "lens-system/lens-element.h"

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
  const Real sin2ThetaI = std::max(
      static_cast<Real>(0.0), static_cast<Real>(1.0) - cosThetaI * cosThetaI);
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

// 近軸光線追跡のレイ
struct ParaxialRay {
  Real u;  // z軸となす角度
  Real h;  // 光線の高さ

  ParaxialRay() {}
  ParaxialRay(Real _u, Real _h) : u(_u), h(_h) {}
};

class LensSystem {
 public:
  std::shared_ptr<Film> film;

  std::vector<LensElement> elements;
  unsigned int aperture_index;

  Real system_length;  // レンズ系の長さ

  Real horizontal_fov;  // 水平FOV [radian]
  Real vertical_fov;    // 垂直FOV [radian]
  Real diagonal_fov;    // 対角FOV [radian]

  Real object_focal_z;       // 物側焦点位置
  Real object_principal_z;   // 物側主点位置
  Real object_focal_length;  // 物側焦点距離
  Real image_focal_z;        // 像側焦点位置
  Real image_principal_z;    // 像側主点位置
  Real image_focal_length;   // 像側焦点距離

  static constexpr unsigned int num_exit_pupil_bounds = 64;
  static constexpr unsigned int num_exit_pupil_bounds_samples = 1024;
  std::vector<Bounds2> exit_pupil_bounds;

  LensSystem(const std::string& filename, const std::shared_ptr<Film>& _film);

  // load lens json
  bool loadJSON(const std::string& filename);

  // compute effective focal length
  Real effective_focal_length() const;
  // compute front focal length
  Real front_focal_length() const;
  // compute back focal length
  Real back_focal_length() const;

  // raytrace
  bool raytrace(const Ray& ray_in, Ray& ray_out, bool reflection = false,
                Sampler* sampler = nullptr) const;
  // raytrace and return raytraced path
  std::vector<Ray> raytracePath(const Ray& ray_in) const;
  // raytrace many rays
  GridData<std::pair<bool, Ray>> raytraceN(const GridData<Ray>& rays_in,
                                           bool reflection = false,
                                           Sampler* sampler = nullptr) const;

  // paraxial raytrace
  std::vector<ParaxialRay> raytraceParaxial(const ParaxialRay& ray_in,
                                            int start = 0, int end = -1,
                                            Real lambda = 550.0) const;

  // compute principal, focal points, focal length
  void computeCardinalPoints();

  // focus lens at z = focus_z
  bool focus(Real focus_z);

  // compute exit pupil bounds at given point on film
  Bounds2 computeExitPupilBound(const Vec2& p) const;
  // compute exit pupil bounds
  bool computeExitPupilBounds();

  // compute exit pupil
  std::pair<GridData<Real>, std::array<Real, 4>> computeExitPupil(
      const Vec2& pFilm, unsigned int n_grids = 512) const;

  // compute spot diagram
  std::vector<Vec3> computeSpotDiagram(const Vec3& origin,
                                       unsigned int n_grids) const;

  std::pair<GridData<Real>, std::array<Real, 4>> computeGeometricPSF(
      const Vec3& origin, unsigned int n_rays = 512,
      unsigned int n_grids = 512) const;

  // compute primary ray
  bool computePrimaryRay(const Vec3& origin, Ray& primary_ray,
                         unsigned int n_grids = 512) const;

  // sample ray going from image sensor to object space
  bool sampleRay(Real u, Real v, Real lambda, Sampler& sampler, Ray& ray_out,
                 Real& pdf, bool reflection = false) const;
};

#endif