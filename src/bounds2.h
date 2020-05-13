#ifndef _BOUNDS2_H
#define _BOUNDS2_H
#include <iostream>

#include "vec2.h"

using namespace Prl2;

class Bounds2 {
 public:
  Vec2 p0;
  Vec2 p1;

  Bounds2() {}
  Bounds2(const Vec2& _p0, const Vec2& _p1) : p0(_p0), p1(_p1) {}
};

inline Bounds2 extendBounds(const Bounds2& b, const Vec2& p) {
  const Real p0x = std::min(b.p0.x(), p.x());
  const Real p0y = std::min(b.p0.y(), p.y());
  const Real p1x = std::max(b.p1.x(), p.x());
  const Real p1y = std::max(b.p1.y(), p.y());
  return Bounds2(Vec2(p0x, p0y), Vec2(p1x, p1y));
}

std::ofstream& operator<<(std::ofstream& stream, const Bounds2& bounds) {
  stream << bounds.p0 << ", " << bounds.p1 << std::endl;
  return stream;
}

#endif