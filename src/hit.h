#ifndef _HIT_H
#define _HIT_H
#include "type.h"
#include "vec3.h"

using namespace Prl2;

struct Hit {
  Real t;
  Vec3 hitPos;
  Vec3 hitNormal;

  Hit(){};
};

#endif