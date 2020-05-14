#ifndef PRL2_RAY_H
#define PRL2_RAY_H

#include <iostream>

#include "core/type.h"
#include "core/vec3.h"

namespace Prl2 {

struct Ray {
  Vec3 origin;     // 始点
  Vec3 direction;  // 方向
  Real lambda;     // 波長[nm]

  static constexpr Real tmin = 1e-6f;  //最小衝突距離
  static constexpr Real tmax = 1e6f;   //最大衝突距離

  Ray() {}
  Ray(const Vec3& _origin, const Vec3& _direction, const Real& _lambda = 550)
      : origin(_origin), direction(_direction), lambda(_lambda) {}

  Vec3 operator()(const Real& t) const { return origin + t * direction; }
};

inline std::ostream& operator<<(std::ostream& stream, const Ray& ray) {
  stream << "origin: " << ray.origin << std::endl;
  stream << "direction: " << ray.direction << std::endl;
  return stream;
}

}  // namespace Prl2

#endif