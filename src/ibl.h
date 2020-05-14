#ifndef _IBL_H
#define _IBL_H

#include <string>

#include "core/vec3.h"

using namespace Prl2;

class IBL {
 public:
  unsigned int width;
  unsigned int height;

  Vec3* pixels;

  IBL(const std::string& filename) {}
};

#endif