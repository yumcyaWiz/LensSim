#ifndef _SELLMEIER_H
#define _SELLMEIER_H

#include <cmath>

#include "core/type.h"

struct SellmeierCofficient {
  Real alpha;
  Real B[3];
  Real C[3];

  SellmeierCofficient(Real ior550, Real b0, Real b1, Real b2, Real c0, Real c1,
                      Real c2) {
    B[0] = b0;
    B[1] = b1;
    B[2] = b2;

    C[0] = c0;
    C[1] = c1;
    C[2] = c2;

    const Real lambda2 = 0.550 * 0.550;
    alpha = ior550 - std::sqrt(1 + B[0] * lambda2 / (lambda2 - C[0]) +
                               B[1] * lambda2 / (lambda2 - C[1]) +
                               B[2] * lambda2 / (lambda2 - C[2]));
  }

  Real ior(Real lambda) const {
    const Real l = lambda * 1e-3f;
    const Real l2 = l * l;
    return std::sqrt(1 + B[0] * l2 / (l2 - C[0]) + B[1] * l2 / (l2 - C[1]) +
                     B[2] * l2 / (l2 - C[2])) +
           alpha;
  }
};

#endif