#ifndef _IOR_H
#define _IOR_H

#include <cmath>

#include "core/type.h"

class IOREquation {
 public:
  virtual Real ior(Real lambda) const = 0;
};

class ConstantIOR : public IOREquation {
 public:
  Real _ior;

  ConstantIOR(Real _ior) : _ior(_ior) {}

  Real ior(Real lambda) const override { return _ior; }
};

class CauthyEquation : public IOREquation {
 public:
  Real A;
  Real B;

  CauthyEquation() {}
  CauthyEquation(Real _A, Real _B) : A(_A), B(_B) {}

  Real ior(Real lambda) const override {
    const Real mu = 1e-3 * lambda;
    return A + B / (mu * mu);
  }
};

inline CauthyEquation fitCauthy(Real nD, Real nF) {
  constexpr Real lambdaD = 0.5893;
  constexpr Real lambdaF = 0.4861;

  const Real alpha1 = 1 / (lambdaD * lambdaD);
  const Real alpha2 = 1 / (lambdaF * lambdaF);
  const Real det = 1 / (alpha2 - alpha1);

  return CauthyEquation(det * (alpha2 * nD - alpha1 * nF), det * (nF - nD));
}

struct SellmeierCofficient : public IOREquation {
 public:
  Real alpha;
  Real B[3];
  Real C[3];

  SellmeierCofficient() {}
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

  Real ior(Real lambda) const override {
    const Real l = lambda * 1e-3f;
    const Real l2 = l * l;
    return std::sqrt(1 + B[0] * l2 / (l2 - C[0]) + B[1] * l2 / (l2 - C[1]) +
                     B[2] * l2 / (l2 - C[2])) +
           alpha;
  }
};

#endif