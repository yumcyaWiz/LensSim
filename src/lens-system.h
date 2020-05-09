#ifndef _LENS_SYSTEM_H
#define _LENS_SYSTEM_H

#include <vector>

#include "lens-element.h"

using namespace Prl2;

class LensSystem {
 public:
  std::vector<LensElement> elements;

  LensSystem(const std::string& filename){};
};

#endif