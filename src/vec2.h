#ifndef PRL2_VEC2_H
#define PRL2_VEC2_H

#include <cassert>
#include <cmath>
#include <iostream>

#include "type.h"

namespace Prl2 {

class alignas(8) Vec2 {
 public:
  explicit Vec2() { v[0] = v[1] = 0; }
  explicit Vec2(const Real& _x) {
    assert(!std::isnan(_x));
    v[0] = v[1] = _x;
  }
  explicit Vec2(const Real& _x, const Real& _y) {
    assert(!std::isnan(_x) && !std::isnan(_y));
    v[0] = _x;
    v[1] = _y;
  }

  Real& operator[](int i) {
    assert(i >= 0 && i < 2);
    return v[i];
  }

  Real x() const { return v[0]; }
  Real y() const { return v[1]; }

 private:
  Real v[2];
};

inline Vec2 operator+(const Vec2& v1, const Vec2& v2) {
  return Vec2(v1.x() + v2.x(), v1.y() + v2.y());
}
inline Vec2 operator+(const Vec2& v, const Real& k) {
  return Vec2(v.x() + k, v.y() + k);
}
inline Vec2 operator+(const Real& k, const Vec2& v) {
  return Vec2(k + v.x(), k + v.y());
}

inline Vec2 operator-(const Vec2& v1, const Vec2& v2) {
  return Vec2(v1.x() - v2.x(), v1.y() - v2.y());
}
inline Vec2 operator-(const Vec2& v, const Real& k) {
  return Vec2(v.x() - k, v.y() - k);
}
inline Vec2 operator-(const Real& k, const Vec2& v) {
  return Vec2(k - v.x(), k - v.y());
}

inline Vec2 operator*(const Vec2& v1, const Vec2& v2) {
  return Vec2(v1.x() * v2.x(), v1.y() * v2.y());
}
inline Vec2 operator*(const Vec2& v, const Real& k) {
  return Vec2(v.x() * k, v.y() * k);
}
inline Vec2 operator*(const Real& k, const Vec2& v) {
  return Vec2(k * v.x(), k * v.y());
}

inline Vec2 operator/(const Vec2& v1, const Vec2& v2) {
  return Vec2(v1.x() / v2.x(), v1.y() / v2.y());
}
inline Vec2 operator/(const Vec2& v, const Real& k) {
  return Vec2(v.x() / k, v.y() / k);
}
inline Vec2 operator/(const Real& k, const Vec2& v) {
  return Vec2(k / v.x(), k / v.y());
}

inline Real length(const Vec2& v) {
  return std::sqrt(v.x() * v.x() + v.y() * v.y());
}
inline Real length2(const Vec2& v) { return v.x() * v.x() + v.y() * v.y(); }

inline Real dot(const Vec2& v1, const Vec2& v2) {
  return v1.x() * v2.x() + v1.y() * v2.y();
}

inline Vec2 normalize(const Vec2& v) { return v / length(v); }

inline Vec2 lerp3(const Real& u, const Real& v, const Vec2& p0, const Vec2& p1,
                  const Vec2& p2) {
  return (1.0f - u - v) * p0 + u * p1 + v * p2;
}

inline std::ostream& operator<<(std::ostream& stream, const Vec2& v) {
  stream << "(" << v.x() << ", " << v.y() << ")";
  return stream;
}

}  // namespace Prl2
#endif