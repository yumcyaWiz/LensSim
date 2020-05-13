#ifndef _BOUNDS2_H
#define _BOUNDS2_H
#include <iostream>
#include <limits>

#include "vec2.h"

using namespace Prl2;

class Bounds2 {
 public:
  Vec2 p0;
  Vec2 p1;

  Bounds2()
      : p0(Vec2(std::numeric_limits<Real>::max())),
        p1(Vec2(std::numeric_limits<Real>::lowest())) {}
  Bounds2(const Vec2& _p0, const Vec2& _p1) : p0(_p0), p1(_p1) {}

  Real area() const { return (p1.x() - p0.x()) * (p1.y() - p0.y()); }
};

inline Bounds2 extendBounds(const Bounds2& b, const Vec2& p) {
  const Real p0x = std::min(b.p0.x(), p.x());
  const Real p0y = std::min(b.p0.y(), p.y());
  const Real p1x = std::max(b.p1.x(), p.x());
  const Real p1y = std::max(b.p1.y(), p.y());
  return Bounds2(Vec2(p0x, p0y), Vec2(p1x, p1y));
}

inline std::ostream& operator<<(std::ostream& stream, const Bounds2& bounds) {
  stream << bounds.p0 << ", " << bounds.p1;
  return stream;
}

#endif