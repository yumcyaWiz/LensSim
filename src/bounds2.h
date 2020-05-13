#ifndef _BOUNDS2_H
#define _BOUNDS2_H
#include "vec2.h"

using namespace Prl2;

class Bounds2 {
 public:
  Vec2 p0;
  Vec2 p1;

  Bounds2(const Vec2& _p0, const Vec2& _p1) : p0(_p0), p1(_p1) {}
}

#endif