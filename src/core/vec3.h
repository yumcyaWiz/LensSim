#ifndef PRL2_VEC3_H
#define PRL2_VEC3_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "core/constant.h"
#include "core/type.h"

namespace Prl2 {

class alignas(16) Vec3 {
 public:
  Real v[3];

  explicit Vec3() { v[0] = v[1] = v[2] = 0; }
  explicit Vec3(const Real& _x) {
    assert(!std::isnan(_x));
    v[0] = v[1] = v[2] = _x;
  }
  explicit Vec3(const Real& _x, const Real& _y, const Real& _z) {
    assert(!std::isnan(_x) && !std::isnan(_y) && !std::isnan(_z));
    v[0] = _x;
    v[1] = _y;
    v[2] = _z;
  }

  Real operator[](int i) const {
    assert(i >= 0 && i < 3);
    return v[i];
  }
  Real& operator[](int i) {
    assert(i >= 0 && i < 3);
    return v[i];
  }

  Real x() const { return v[0]; }
  Real y() const { return v[1]; }
  Real z() const { return v[2]; }

  Vec3 operator-() const { return Vec3(-x(), -y(), -z()); }

  Vec3& operator+=(const Vec3& v2) {
    v[0] += v2.x();
    v[1] += v2.y();
    v[2] += v2.z();
    return *this;
  }
  Vec3& operator+=(const Real& k) {
    v[0] += k;
    v[1] += k;
    v[2] += k;
    return *this;
  }
  Vec3& operator-=(const Vec3& v2) {
    v[0] -= v2.x();
    v[1] -= v2.y();
    v[2] -= v2.z();
    return *this;
  }
  Vec3& operator-=(const Real& k) {
    v[0] -= k;
    v[1] -= k;
    v[2] -= k;
    return *this;
  }
  Vec3& operator*=(const Vec3& v2) {
    v[0] *= v2.x();
    v[1] *= v2.y();
    v[2] *= v2.z();
    return *this;
  }
  Vec3& operator*=(const Real& k) {
    v[0] *= k;
    v[1] *= k;
    v[2] *= k;
    return *this;
  }
  Vec3& operator/=(const Vec3& v2) {
    v[0] /= v2.x();
    v[1] /= v2.y();
    v[2] /= v2.z();
    return *this;
  }
  Vec3& operator/=(const Real& k) {
    v[0] /= k;
    v[1] /= k;
    v[2] /= k;
    return *this;
  }
};

inline Vec3 operator+(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}
inline Vec3 operator+(const Vec3& v, const Real& k) {
  return Vec3(v.x() + k, v.y() + k, v.z() + k);
}
inline Vec3 operator+(const Real& k, const Vec3& v) {
  return Vec3(k + v.x(), k + v.y(), k + v.z());
}

inline Vec3 operator-(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}
inline Vec3 operator-(const Vec3& v, const Real& k) {
  return Vec3(v.x() - k, v.y() - k, v.z() - k);
}
inline Vec3 operator-(const Real& k, const Vec3& v) {
  return Vec3(k - v.x(), k - v.y(), k - v.z());
}

inline Vec3 operator*(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x() * v2.x(), v1.y() * v2.y(), v1.z() * v2.z());
}
inline Vec3 operator*(const Vec3& v, const Real& k) {
  return Vec3(v.x() * k, v.y() * k, v.z() * k);
}
inline Vec3 operator*(const Real& k, const Vec3& v) {
  return Vec3(k * v.x(), k * v.y(), k * v.z());
}

inline Vec3 operator/(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x() / v2.x(), v1.y() / v2.y(), v1.z() / v2.z());
}
inline Vec3 operator/(const Vec3& v, const Real& k) {
  return Vec3(v.x() / k, v.y() / k, v.z() / k);
}
inline Vec3 operator/(const Real& k, const Vec3& v) {
  return Vec3(k / v.x(), k / v.y(), k / v.z());
}

inline Real length(const Vec3& v) {
  return std::sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
}
inline Real length2(const Vec3& v) {
  return v.x() * v.x() + v.y() * v.y() + v.z() * v.z();
}

inline Real dot(const Vec3& v1, const Vec3& v2) {
  return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}
inline Vec3 cross(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.y() * v2.z() - v1.z() * v2.y(),
              v1.z() * v2.x() - v1.x() * v2.z(),
              v1.x() * v2.y() - v1.y() * v2.x());
}

inline Vec3 normalize(const Vec3& v) { return v / length(v); }

inline std::ostream& operator<<(std::ostream& stream, const Vec3& v) {
  stream << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
  return stream;
}

// 方向ベクトルから球面座標系を計算
inline void cartesianToSpherical(const Vec3& v, Real& theta, Real& phi) {
  phi = std::atan2(v.z(), v.x());
  if (phi < 0) phi += PI_MUL_2;
  theta = std::acos(std::clamp(v.y(), -1.0f, 1.0f));
}

// 球面座標から方向ベクトルを生成
inline Vec3 sphericalToCartesian(const Real& theta, const Real& phi) {
  const Real sinTheta = std::sin(theta);
  return Vec3(std::cos(phi) * sinTheta, std::cos(theta),
              std::sin(phi) * sinTheta);
}

// ２つの方向ベクトルの間の角度を計算
inline Real radianBetween(const Vec3& v1, const Vec3& v2) {
  const Real cos = dot(v1, v2);
  return std::acos(cos);
}

inline Vec3 clamp(const Vec3& v, const Vec3& vmin, const Vec3& vmax) {
  Vec3 ret;
  ret[0] = std::clamp(v.x(), vmin.x(), vmax.x());
  ret[1] = std::clamp(v.y(), vmin.y(), vmax.y());
  ret[2] = std::clamp(v.z(), vmin.z(), vmax.z());
  return ret;
}

// 正規直交基底を作る
// Duff et al.
// https://shikihuiku.wordpress.com/2018/07/09/%E6%AD%A3%E8%A6%8F%E7%9B%B4%E4%BA%A4%E5%9F%BA%E5%BA%95%E3%81%AE%E4%BD%9C%E3%82%8A%E6%96%B9%E3%81%AB%E3%81%A4%E3%81%84%E3%81%A6%E3%80%81%E6%94%B9%E3%82%81%E3%81%A6%E5%8B%89%E5%BC%B7%E3%81%97%E3%81%BE/
inline void orthonormalBasis(const Vec3& n, Vec3& s, Vec3& t) {
  const Real sign = std::copysign(1.0f, n.z());
  const Real a = -1.0f / (sign + n.z());
  const Real b = n.x() * n.y() * a;
  s = Vec3(1.0f + sign * n.x() * n.x() * a, sign * b, -sign * n.x());
  t = Vec3(b, sign + n.y() * n.y() * a, -n.y());
}

inline Vec3 lerp3(const Real& u, const Real& v, const Vec3& p0, const Vec3& p1,
                  const Vec3& p2) {
  return (1.0f - u - v) * p0 + u * p1 + v * p2;
}

}  // namespace Prl2
#endif
