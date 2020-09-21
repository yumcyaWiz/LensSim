#include "lens-system/lens-system.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include "nlohmann/json.hpp"
using JSON = nlohmann::json;

#include "parallel/parallel.h"

LensSystem::LensSystem(const std::string& filename,
                       const std::shared_ptr<Film>& _film)
    : film(_film) {
  // load json
  if (!loadJSON(filename)) exit(EXIT_FAILURE);

  // compute aperture index
  for (int i = 0; i < elements.size(); ++i) {
    if (elements[i].is_stop) {
      aperture_index = i;
      break;
    }
  }

  // compute system length and z
  Real _length = 0;
  for (auto itr = elements.rbegin(); itr != elements.rend(); itr++) {
    _length += (*itr).thickness;
    (*itr).z = -_length;
  }
  system_length = _length;

  // compute cardinal points
  computeCardinalPoints();

  // compute fov
  horizontal_fov =
      2.0f * std::atan2(film->width_length, 2.0f * image_focal_length);
  vertical_fov =
      2.0f * std::atan2(film->height_length, 2.0f * image_focal_length);
  diagonal_fov =
      2.0f * std::atan2(film->diagonal_length, 2.0f * image_focal_length);
}

bool LensSystem::loadJSON(const std::string& filename) {
  // open file
  std::ifstream stream(filename);
  if (!stream) {
    std::cerr << "failed to open " << filename << std::endl;
    return false;
  }

  // parse JSON
  JSON json;
  stream >> json;
  stream.close();

  // push lens elements
  for (const auto& [key, value] : json.items()) {
    // required
    const unsigned int index = value["index"].get<unsigned int>();
    const Real curvature_radius = value["curvature_radius"].get<Real>();
    const Real thickness = value["thickness"].get<Real>();
    const Real aperture_radius = 0.5f * value["aperture_diameter"].get<Real>();

    // optional
    bool is_stop = false;
    if (value.count("is_stop") > 0) {
      is_stop = value["is_stop"].get<bool>();
    }

    std::optional<Real> nd;
    if (value.count("nd") > 0) {
      nd = value["nd"].get<Real>();
    }

    std::optional<Real> nD;
    if (value.count("nD") > 0) {
      nD = value["nD"].get<Real>();
    }

    std::optional<Real> nF;
    if (value.count("nF") > 0) {
      nF = value["nF"].get<Real>();
    }

    // select ior equation
    std::shared_ptr<IOREquation> ior_equation;
    if (nd) {
      ior_equation = std::make_shared<ConstantIOR>(nd.value());
    } else if (nD && nF) {
      ior_equation =
          std::make_shared<CauthyEquation>(fitCauthy(nD.value(), nF.value()));
    } else {
      std::cerr << "failed to create ior equation for this lens element"
                << std::endl;
      return false;
    }

    // make lens element
    const auto element = LensElement(index, aperture_radius, thickness,
                                     curvature_radius, ior_equation, is_stop);
    elements.push_back(element);
  }

  // sort lens elements by index
  std::sort(elements.begin(), elements.end(),
            [](const LensElement& x1, const LensElement& x2) {
              return x1.index < x2.index;
            });

  return true;
}

Real LensSystem::effective_focal_length() const { return image_focal_length; }
Real LensSystem::front_focal_length() const {
  return std::abs(elements.front().z - object_focal_z);
}
Real LensSystem::back_focal_length() const {
  return std::abs(image_focal_z - elements.back().z);
}

bool LensSystem::raytrace(const Ray& ray_in, Ray& ray_out, bool reflection,
                          Sampler* sampler) const {
  int element_index = ray_in.direction.z() > 0 ? -1 : elements.size();
  const int initial_element_index = element_index;

  Ray ray = ray_in;
  Real ior = 1.0f;

  while (true) {
    // update element index
    element_index += ray.direction.z() > 0 ? 1 : -1;
    if (element_index < 0 || element_index >= elements.size()) break;
    const LensElement& element = elements[element_index];

    // Aperture
    if (element.is_stop) {
      Hit res;
      if (!element.intersect(ray, res)) return false;

      // Update ray
      ray.origin = res.hitPos;

      // Update ior
      ior = 1.0f;
    }
    // Lens
    else {
      // Compute Next Element IOR
      Real next_ior = 1.0f;

      // Compute Next Element
      const int next_element_index =
          ray.direction.z() > 0 ? element_index : element_index - 1;
      if (next_element_index >= 0) {
        const LensElement& next_element = elements[next_element_index];
        if (!next_element.is_stop) {
          next_ior = next_element.ior(ray.lambda);
        }
      }

      // Compute Intersection with Lens
      Hit res;
      if (!element.intersect(ray, res)) return false;

      // Refract and Reflect
      Vec3 next_direction;
      if (reflection) {
        const Real fr = fresnel(-ray.direction, res.hitNormal, ior, next_ior);
        if (sampler->getNext() < fr) {
          // reflection
          next_direction = reflect(-ray.direction, res.hitNormal);
        } else {
          // refract
          if (!refract(-ray.direction, next_direction, res.hitNormal, ior,
                       next_ior)) {
            // total reflection
            next_direction = reflect(-ray.direction, res.hitNormal);
          }
        }
      } else {
        if (!refract(-ray.direction, next_direction, res.hitNormal, ior,
                     next_ior))
          return false;
      }

      // Set Next Ray
      ray.origin = res.hitPos;
      ray.direction = normalize(next_direction);

      // update ior
      ior = next_ior;
    }
  }

  // if ray exits from the same side
  if (element_index == initial_element_index) return false;

  ray_out = ray;

  return true;
}

std::vector<Ray> LensSystem::raytracePath(const Ray& ray_in) const {
  std::vector<Ray> ret;

  int element_index = ray_in.direction.z() > 0 ? -1 : elements.size();
  const int initial_element_index = element_index;

  Ray ray = ray_in;
  Real ior = 1.0f;
  while (true) {
    // push ray origin
    ret.push_back(ray);

    // update element index
    element_index += ray.direction.z() > 0 ? 1 : -1;
    if (element_index < 0 || element_index >= elements.size()) break;
    const LensElement& element = elements[element_index];

    // Aperture
    if (element.is_stop) {
      Hit res;
      if (!element.intersect(ray, res)) break;

      // Update ray
      ray.origin = res.hitPos;

      // Update ior
      ior = 1.0f;
    }
    // Lens
    else {
      // Compute Next Element IOR
      Real next_ior = 1.0f;

      // Compute Next Element
      const int next_element_index =
          ray.direction.z() > 0 ? element_index : element_index - 1;
      if (next_element_index >= 0) {
        const LensElement& next_element = elements[next_element_index];
        if (!next_element.is_stop) {
          next_ior = next_element.ior(ray.lambda);
        }
      }

      // Compute Intersection with Lens
      Hit res;
      if (!element.intersect(ray, res)) break;

      // Refract and Reflect
      Vec3 next_direction;
      if (!refract(-ray.direction, next_direction, res.hitNormal, ior,
                   next_ior))
        break;

      // Set Next Ray
      ray.origin = res.hitPos;
      ray.direction = normalize(next_direction);

      // update ior
      ior = next_ior;
    }
  }

  return ret;
}

std::vector<ParaxialRay> LensSystem::raytraceParaxial(const ParaxialRay& ray_in,
                                                      int start, int end,
                                                      Real lambda) const {
  std::vector<ParaxialRay> ret;
  ret.push_back(ray_in);

  // compute start, end index
  int start_index = start;
  if (start_index == -1) {
    start_index = elements.size() - 1;
  }
  int end_index = end;
  if (end_index == -1) {
    end_index = elements.size() - 1;
  }

  // invalid index case
  if (start_index < 0 || start_index >= elements.size() || end_index < 0 ||
      end_index >= elements.size()) {
    std::cerr << "invalid start, end index" << std::endl;
    return ret;
  }

  // paraxial raytrace
  Real ior;
  Real ior_prev = (start_index == 0 || start_index == elements.size() - 1)
                      ? 1.0
                      : elements[start_index].ior(lambda);
  Real u_prev = ray_in.u;
  Real h_prev = ray_in.h;

  if (start_index < end_index) {
    for (int i = start_index; i <= end_index; ++i) {
      const auto& element = elements[i];

      // compute curvature radius
      Real r = element.curvature_radius;
      if (element.is_stop) {
        r = 1e9;
      }

      // compute ior of element
      ior = element.ior(lambda);

      // compute thickness
      Real thickness = 0;
      if (i != elements.size() - 1) {
        thickness = element.thickness;
      }

      // compute paraxial ray
      const Real u =
          ior_prev / ior * u_prev + (ior - ior_prev) / (ior * r) * h_prev;
      const Real h = h_prev - thickness * u;

      // save paraxial ray
      ret.push_back(ParaxialRay(u, h));

      // update
      ior_prev = ior;
      u_prev = u;
      h_prev = h;
    }
  }
  // reverse paraxial raytrace
  else {
    for (int i = start_index; i >= end_index; --i) {
      // compute curvature radius
      Real r = -elements[i].curvature_radius;
      if (elements[i].is_stop) {
        r = -1e9;
      }

      // compute ior of element
      if (i > 0) {
        ior = elements[i - 1].ior(lambda);
      } else {
        ior = 1.0;
      }

      // compute thickness
      Real thickness = 0;
      if (i > 0) {
        thickness = elements[i - 1].thickness;
      }

      // compute paraxial ray
      const Real u =
          ior_prev / ior * u_prev + (ior - ior_prev) / (ior * r) * h_prev;
      const Real h = h_prev - thickness * u;

      // save paraxial ray
      ret.push_back(ParaxialRay(u, h));

      // update
      ior_prev = ior;
      u_prev = u;
      h_prev = h;
    }
  }

  return ret;
}

void LensSystem::computeCardinalPoints() {
  // compute image focal point
  // paraxial raytrace with (u, h) = (0, 1)
  auto result = raytraceParaxial(ParaxialRay(0, 1));
  image_focal_z = elements.back().z + result.back().h / result.back().u;

  // compute principal point
  image_principal_z = elements.back().z +
                      (result.back().h - result.front().h) / result.back().u;

  // compute image focal length
  image_focal_length = image_focal_z - image_principal_z;

  // compute object focal point
  // paraxial reverse raytrace with (u, h) = (0, 1)
  result = raytraceParaxial(ParaxialRay(0, 1), -1, 0);
  object_focal_z = elements.front().z - result.back().h / result.back().u;

  // compute object principal point
  object_principal_z = elements.front().z -
                       (result.back().h - result.front().h) / result.back().u;

  // compute object focal length
  object_focal_length = -(object_focal_z - object_principal_z);
}

bool LensSystem::focus(Real focus_z) {
  const Real delta =
      0.5f * (object_principal_z - focus_z + image_principal_z -
              std::sqrt((object_principal_z - focus_z - image_principal_z) *
                        (object_principal_z - focus_z - 4 * image_focal_length -
                         image_principal_z)));

  // move lens elements
  for (auto& element : elements) {
    element.z -= delta;
  }

  // recompute cardinal points
  computeCardinalPoints();

  return true;
}

Bounds2 LensSystem::computeExitPupilBound(const Vec2& p) const {
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
          Vec3(lastElement.aperture_radius * u, lastElement.aperture_radius * v,
               lastElement.z);

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

bool LensSystem::computeExitPupilBounds() {
  exit_pupil_bounds.resize(num_exit_pupil_bounds);

  Parallel parallel;

  parallel.parallelFor1D(
      [&](unsigned int idx) {
        const Real r = static_cast<Real>(idx) / num_exit_pupil_bounds * 0.5f *
                       film->diagonal_length;
        exit_pupil_bounds[idx] = computeExitPupilBound(Vec2(r, 0));

        std::cout << "finished " << idx << "th computation of exit pupil bounds"
                  << std::endl;
        std::cout << exit_pupil_bounds[idx] << std::endl;
      },
      16, num_exit_pupil_bounds);

  return true;
}

bool LensSystem::sampleRay(Real u, Real v, Real lambda, Sampler& sampler,
                           Ray& ray_out, Real& pdf, bool reflection) const {
  // compute position on film
  const Vec2 p = film->computePosition(u, v);

  // choose exit pupil bound
  const Real r = length(p);
  const unsigned int exit_pupil_bounds_index =
      r / (0.5f * film->diagonal_length) * num_exit_pupil_bounds;
  const Bounds2& exit_pupil_bound = exit_pupil_bounds[exit_pupil_bounds_index];
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
  const Vec3 pBound3 = Vec3(pBound.x(), pBound.y(), elements.back().z);
  const Vec3 direction = normalize(pBound3 - origin);
  const Ray ray_in(origin, direction, lambda);

  // convert area pdf to solid angle pdf
  const Real l = length(pBound3 - origin);
  pdf = l * l / std::abs(dot(direction, Vec3(0, 0, -1))) * pdf_area;

  // raytrace
  if (!raytrace(ray_in, ray_out, reflection, &sampler)) return false;

  return true;
}

GridData<std::pair<bool, Ray>> LensSystem::raytraceN(
    const GridData<Ray>& rays_in, bool reflection, Sampler* sampler) const {
  GridData<std::pair<bool, Ray>> ret(rays_in.nrows, rays_in.ncols);

  Parallel parallel;
  parallel.parallelFor2D(
      [&](unsigned int i, unsigned int j) {
        bool traced;
        Ray ray_out;
        traced = raytrace(rays_in.get(i, j), ray_out, reflection, sampler);

        ret.set(i, j, {traced, ray_out});
      },
      16, 16, rays_in.nrows, rays_in.ncols);

  return ret;
}

std::pair<GridData<Real>, std::array<Real, 4>> LensSystem::computeExitPupil(
    const Vec2& pFilm, unsigned int n_grids) const {
  const auto& lastElement = elements.back();

  // compute extends
  std::array<Real, 4> extends;
  extends[0] = -lastElement.aperture_radius;
  extends[1] = lastElement.aperture_radius;
  extends[2] = -lastElement.aperture_radius;
  extends[3] = lastElement.aperture_radius;

  // compute grids
  const GridData<Vec3> grids = elements.back().samplePoints(n_grids);

  // make rays
  GridData<Ray> rays_in(n_grids, n_grids);
  for (int i = 0; i < n_grids; ++i) {
    for (int j = 0; j < n_grids; ++j) {
      const Vec3 pFilm3 = Vec3(pFilm.x(), pFilm.y(), 0);
      rays_in.set(i, j, Ray(pFilm3, normalize(grids.get(i, j) - pFilm3)));
    }
  }

  // raytrace
  const auto result = raytraceN(rays_in);

  // represent exit pupil as 0, 1
  GridData<Real> exit_pupil(result.nrows, result.ncols);
  for (int i = 0; i < result.nrows; ++i) {
    for (int j = 0; j < result.ncols; ++j) {
      exit_pupil.set(i, j, result.get(i, j).first);
    }
  }

  return {exit_pupil, extends};
}

bool LensSystem::computePrimaryRay(const Vec3& origin, Ray& primary_ray,
                                   unsigned int n_grids) const {
  // compute grids
  const GridData<Vec3> grids = elements.front().samplePoints(n_grids);

  // make rays
  GridData<Ray> rays_in(n_grids, n_grids);
  for (int i = 0; i < n_grids; ++i) {
    for (int j = 0; j < n_grids; ++j) {
      rays_in.set(i, j, Ray(origin, normalize(grids.get(i, j) - origin)));
    }
  }

  // raytrace
  const auto result = raytraceN(rays_in);

  // compute entrance pupil center
  unsigned int n_average = 0;
  Vec3 entrance_pupil_center;
  for (int i = 0; i < n_grids; ++i) {
    for (int j = 0; j < n_grids; ++j) {
      if (result.get(i, j).first) {
        entrance_pupil_center += grids.get(i, j);
        n_average++;
      }
    }
  }
  if (n_average == 0) {
    std::cerr << "failed to compute primary ray" << std::endl;
    return false;
  }
  entrance_pupil_center /= n_average;

  primary_ray = Ray(origin, normalize(entrance_pupil_center - origin));

  return true;
}

std::vector<Vec3> LensSystem::computeSpotDiagram(const Vec3& origin,
                                                 unsigned int n_grids) const {
  std::vector<Vec3> ret;

  // compute grids
  const GridData<Vec3> grids = elements.front().samplePoints(n_grids);

  // make rays
  GridData<Ray> rays_in(n_grids, n_grids);
  for (int i = 0; i < n_grids; ++i) {
    for (int j = 0; j < n_grids; ++j) {
      rays_in.set(i, j, Ray(origin, normalize(grids.get(i, j) - origin)));
    }
  }

  // raytrace
  const auto result = raytraceN(rays_in);

  // compute intersect position at gaussian plane
  for (int i = 0; i < n_grids; ++i) {
    for (int j = 0; j < n_grids; ++j) {
      if (result.get(i, j).first) {
        const Ray& ray = result.get(i, j).second;
        const Real t = -(ray.origin.z() - image_focal_z) / ray.direction.z();
        const Vec3 pFilm = ray(t);
        ret.push_back(pFilm);
      }
    }
  }

  return ret;
}

std::pair<GridData<Real>, std::array<Real, 4>> LensSystem::computeGeometricPSF(
    const Vec3& origin, unsigned int n_rays, unsigned int n_grids) const {
  // compute grids
  const GridData<Vec3> grids = elements.front().samplePoints(n_rays);

  // make rays
  GridData<Ray> rays_in(n_rays, n_rays);
  for (int i = 0; i < n_rays; ++i) {
    for (int j = 0; j < n_rays; ++j) {
      rays_in.set(i, j, Ray(origin, normalize(grids.get(i, j) - origin)));
    }
  }

  // raytrace
  const auto result = raytraceN(rays_in);

  // compute intersect position at gaussian plane
  std::vector<Vec3> points;
  for (int i = 0; i < n_rays; ++i) {
    for (int j = 0; j < n_rays; ++j) {
      if (result.get(i, j).first) {
        const Ray& ray = result.get(i, j).second;
        const Real t = -(ray.origin.z() - image_focal_z) / ray.direction.z();
        const Vec3 pFilm = ray(t);
        points.push_back(pFilm);
      }
    }
  }

  // compute primary ray
  Ray primary_ray;
  if (!computePrimaryRay(origin, primary_ray)) {
    std::cerr << "failed to compute primary ray" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // raytrace
  Ray primary_ray_out;
  if (!raytrace(primary_ray, primary_ray_out)) {
    std::cerr << "failed to raytrace primary ray" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // compute intersect position at gaussian plane
  const Vec3 pFilm_primary =
      primary_ray_out(-(primary_ray_out.origin.z() - image_focal_z) /
                      primary_ray_out.direction.z());

  // compute mean, variance
  const Vec3 p_mean =
      std::accumulate(points.begin(), points.end(), Vec3()) / points.size();
  const Vec3 p_var =
      std::accumulate(points.begin(), points.end(), Vec3(),
                      [&](const Vec3& v1, const Vec3& v2) {
                        return v1 + (v2 - p_mean) * (v2 - p_mean);
                      }) /
      points.size();

  // compute extent
  const Real xstd = std::sqrt(p_var.x());
  const Real ystd = std::sqrt(p_var.y());
  const Real grid_width = xstd > ystd ? 3 * xstd : 3 * ystd;
  const Real xmin = pFilm_primary.x() - grid_width;
  const Real xmax = pFilm_primary.x() + grid_width;
  const Real ymin = pFilm_primary.y() - grid_width;
  const Real ymax = pFilm_primary.y() + grid_width;
  const std::array<Real, 4> extent = {xmin, xmax, ymin, ymax};

  // init
  GridData<Real> psf(n_grids, n_grids);
  for (int i = 0; i < n_grids; ++i) {
    for (int j = 0; j < n_grids; ++j) {
      psf.set(i, j, 0);
    }
  }

  // compute psf by grid
  for (const auto& p : points) {
    // compute grid index
    const unsigned int i = (p.x() - xmin) / (xmax - xmin) * n_grids;
    const unsigned int j = (p.y() - ymin) / (ymax - ymin) * n_grids;

    // accumulate
    if (i >= 0 && i < n_grids && j >= 0 && j < n_grids)
      psf.set(j, i, psf.get(j, i) + 1);
  }

  // normalize
  const Real max_value = *std::max_element(psf.data.begin(), psf.data.end());
  for (int i = 0; i < psf.data.size(); ++i) {
    psf.data[i] /= max_value;
  }

  return {psf, extent};
}