#ifndef PRL2_CONSTANT_H
#define PRL2_CONSTANT_H

#include <cmath>
#include <limits>

#include "core/type.h"

namespace Prl2 {

static constexpr Real PI = static_cast<Real>(3.141592653589793238462643383279);
static constexpr Real PI_MUL_2 = PI * 2;
static constexpr Real PI_MUL_4 = PI * 4;
static constexpr Real PI_DIV_2 = PI / 2;

static constexpr Real INV_PI = 1 / PI;
static constexpr Real INV_PI_MUL_2 = 1 / PI_MUL_2;
static constexpr Real INV_PI_MUL_4 = 1 / PI_MUL_4;

static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
static constexpr Real ONE_MINUS_EPS =
    static_cast<Real>(1) - std::numeric_limits<Real>::epsilon();

static constexpr Real INF = std::numeric_limits<Real>::infinity();

}  // namespace Prl2

#endif