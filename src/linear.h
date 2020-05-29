#ifndef _LINEAR_H
#define _LINEAR_H

#include <memory>
#include <vector>

#include "core/hit.h"
#include "core/ray.h"
#include "sphere.h"

class LinearIntersector {
 public:
  std::vector<std::shared_ptr<Sphere>> prims;

  LinearIntersector() {}

  void add(const std::shared_ptr<Sphere>& prim) { prims.push_back(prim); }

  bool intersect(const Ray& ray, Hit& res) const {
    bool hit = false;
    res.t = 1e9;
    for (const auto& prim : prims) {
      Hit tmp;
      if (prim->intersect(ray, tmp)) {
        if (tmp.t < res.t) {
          hit = true;
          res = tmp;
        }
      }
    }
    return hit;
  }
};

#endif