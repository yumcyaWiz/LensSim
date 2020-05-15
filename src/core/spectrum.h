#ifndef PRL2_SPECTRUM_H
#define PRL2_SPECTRUM_H

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <vector>

#include "core/type.h"
#include "core/vec3.h"

namespace Prl2 {

using XYZ = Vec3;
using RGB = Vec3;

// XYZをsRGB色空間に変換する
// XYZ to sRGB(D65)
// http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
inline RGB XYZ2RGB(const XYZ& xyz) {
  return RGB(
      3.2404542f * xyz.x() - 1.5371385f * xyz.y() - 0.4985314f * xyz.z(),
      -0.9692660f * xyz.x() + 1.8760108f * xyz.y() + 0.0415560f * xyz.z(),
      0.0556434f * xyz.x() - 0.2040259f * xyz.y() + 1.0572252f * xyz.z());
}

//等間隔にサンプリングされたSPDを表現する
//波長と放射束のサンプリング列を保持する
//波長は[nm]で保持する
//波長の分割幅より狭いピークを持つSPDは適切に表現されない可能性がある
//本当はサンプルをそのまま保持して非等間隔のSPDを表現できるようにしたかったが、データサイズが大きすぎるので諦めた
class SPD {
 public:
  // SPDに格納する波長の範囲
  static constexpr Real LAMBDA_MIN = 380;
  static constexpr Real LAMBDA_MAX = 780;

  //波長の分割数
  static constexpr size_t LAMBDA_SAMPLES = 80;

  //分割された波長幅
  static constexpr Real LAMBDA_INTERVAL =
      (LAMBDA_MAX - LAMBDA_MIN) / LAMBDA_SAMPLES;

  std::array<Real, LAMBDA_SAMPLES> phi;  //放射束

  // 0で初期化
  SPD() {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] = 0;
    }
  }

  // ある値で初期化
  SPD(const Real& v) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] = v;
    }
  }

  //任意の波長と放射束のサンプリング列から等間隔のSPDを構築
  //波長と対応する放射束は昇順で並んでいると仮定している
  SPD(const std::vector<Real>& _lambda, const std::vector<Real>& _phi);

  // i番目の放射束を返す
  Real operator[](size_t i) const {
    assert(i < SPD::LAMBDA_SAMPLES);
    return phi[i];
  }

  // クリアする
  void clear() {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] = 0;
    }
  }

  //分光放射束を加算する
  void addPhi(const Real& _lambda, const Real& _phi);

  //黒色か返す
  bool isBlack() const {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      if (phi[i] != 0.0f) {
        return false;
      }
    }
    return true;
  }

  //指定した波長の放射束を線形補間して返す
  // l : 波長[nm]
  Real sample(const Real& l) const;

  // XYZ色空間に変換する
  XYZ toXYZ() const;

  // sRGB色空間に変換する
  // XYZ to sRGB(D65)
  // http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
  RGB toRGB() const {
    XYZ xyz = this->toXYZ();
    RGB rgb = XYZ2RGB(xyz);
    return clamp(rgb, Vec3(0), Vec3(INF));
  }

  //演算
  SPD& operator+=(const SPD& spd) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] += spd.phi[i];
    }
    return *this;
  }
  SPD& operator+=(const Real& k) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] += k;
    }
    return *this;
  }
  SPD& operator-=(const SPD& spd) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] -= spd.phi[i];
    }
    return *this;
  }
  SPD& operator-=(const Real& k) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] -= k;
    }
    return *this;
  }
  SPD& operator*=(const SPD& spd) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] *= spd.phi[i];
    }
    return *this;
  }
  SPD& operator*=(const Real& k) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] *= k;
    }
    return *this;
  }
  SPD& operator/=(const SPD& spd) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] /= spd.phi[i];
    }
    return *this;
  }
  SPD& operator/=(const Real& k) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] /= k;
    }
    return *this;
  }

  //等色関数(CIE1931)
  // http://cvrl.ucl.ac.uk/cmfs.htm
  static constexpr int color_matching_func_samples = 85;
  static constexpr Real color_matching_func_x[color_matching_func_samples] = {
      0.001368000000f, 0.002236000000f, 0.004243000000f, 0.007650000000f,
      0.014310000000f, 0.023190000000f, 0.043510000000f, 0.077630000000f,
      0.134380000000f, 0.214770000000f, 0.283900000000f, 0.328500000000f,
      0.348280000000f, 0.348060000000f, 0.336200000000f, 0.318700000000f,
      0.290800000000f, 0.251100000000f, 0.195360000000f, 0.142100000000f,
      0.095640000000f, 0.057950010000f, 0.032010000000f, 0.014700000000f,
      0.004900000000f, 0.002400000000f, 0.009300000000f, 0.029100000000f,
      0.063270000000f, 0.109600000000f, 0.165500000000f, 0.225749900000f,
      0.290400000000f, 0.359700000000f, 0.433449900000f, 0.512050100000f,
      0.594500000000f, 0.678400000000f, 0.762100000000f, 0.842500000000f,
      0.916300000000f, 0.978600000000f, 1.026300000000f, 1.056700000000f,
      1.062200000000f, 1.045600000000f, 1.002600000000f, 0.938400000000f,
      0.854449900000f, 0.751400000000f, 0.642400000000f, 0.541900000000f,
      0.447900000000f, 0.360800000000f, 0.283500000000f, 0.218700000000f,
      0.164900000000f, 0.121200000000f, 0.087400000000f, 0.063600000000f,
      0.046770000000f, 0.032900000000f, 0.022700000000f, 0.015840000000f,
      0.011359160000f, 0.008110916000f, 0.005790346000f, 0.004109457000f,
      0.002899327000f, 0.002049190000f, 0.001439971000f, 0.000999949300f,
      0.000690078600f, 0.000476021300f, 0.000332301100f, 0.000234826100f,
      0.000166150500f, 0.000117413000f, 0.000083075270f, 0.000058706520f,
      0.000041509940f};
  static constexpr Real color_matching_func_y[color_matching_func_samples] = {
      0.000039000000f, 0.000064000000f, 0.000120000000f, 0.000217000000f,
      0.000396000000f, 0.000640000000f, 0.001210000000f, 0.002180000000f,
      0.004000000000f, 0.007300000000f, 0.011600000000f, 0.016840000000f,
      0.023000000000f, 0.029800000000f, 0.038000000000f, 0.048000000000f,
      0.060000000000f, 0.073900000000f, 0.090980000000f, 0.112600000000f,
      0.139020000000f, 0.169300000000f, 0.208020000000f, 0.258600000000f,
      0.323000000000f, 0.407300000000f, 0.503000000000f, 0.608200000000f,
      0.710000000000f, 0.793200000000f, 0.862000000000f, 0.914850100000f,
      0.954000000000f, 0.980300000000f, 0.994950100000f, 1.000000000000f,
      0.995000000000f, 0.978600000000f, 0.952000000000f, 0.915400000000f,
      0.870000000000f, 0.816300000000f, 0.757000000000f, 0.694900000000f,
      0.631000000000f, 0.566800000000f, 0.503000000000f, 0.441200000000f,
      0.381000000000f, 0.321000000000f, 0.265000000000f, 0.217000000000f,
      0.175000000000f, 0.138200000000f, 0.107000000000f, 0.081600000000f,
      0.061000000000f, 0.044580000000f, 0.032000000000f, 0.023200000000f,
      0.017000000000f, 0.011920000000f, 0.008210000000f, 0.005723000000f,
      0.004102000000f, 0.002929000000f, 0.002091000000f, 0.001484000000f,
      0.001047000000f, 0.000740000000f, 0.000520000000f, 0.000361100000f,
      0.000249200000f, 0.000171900000f, 0.000120000000f, 0.000084800000f,
      0.000060000000f, 0.000042400000f, 0.000030000000f, 0.000021200000f,
      0.000014990000f};
  static constexpr Real color_matching_func_z[color_matching_func_samples] = {
      0.006450001000f, 0.010549990000f, 0.020050010000f, 0.036210000000f,
      0.067850010000f, 0.110200000000f, 0.207400000000f, 0.371300000000f,
      0.645600000000f, 1.039050100000f, 1.385600000000f, 1.622960000000f,
      1.747060000000f, 1.782600000000f, 1.772110000000f, 1.744100000000f,
      1.669200000000f, 1.528100000000f, 1.287640000000f, 1.041900000000f,
      0.812950100000f, 0.616200000000f, 0.465180000000f, 0.353300000000f,
      0.272000000000f, 0.212300000000f, 0.158200000000f, 0.111700000000f,
      0.078249990000f, 0.057250010000f, 0.042160000000f, 0.029840000000f,
      0.020300000000f, 0.013400000000f, 0.008749999000f, 0.005749999000f,
      0.003900000000f, 0.002749999000f, 0.002100000000f, 0.001800000000f,
      0.001650001000f, 0.001400000000f, 0.001100000000f, 0.001000000000f,
      0.000800000000f, 0.000600000000f, 0.000340000000f, 0.000240000000f,
      0.000190000000f, 0.000100000000f, 0.000049999990f, 0.000030000000f,
      0.000020000000f, 0.000010000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f, 0.000000000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f, 0.000000000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f, 0.000000000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f, 0.000000000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f, 0.000000000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f, 0.000000000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f};
};

// SPDどうしの演算
//要素ごとに演算を行う
inline SPD operator+(const SPD& spd1, const SPD& spd2) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd1.phi[i] + spd2.phi[i];
  }
  return ret;
}
inline SPD operator-(const SPD& spd1, const SPD& spd2) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd1.phi[i] - spd2.phi[i];
  }
  return ret;
}
inline SPD operator*(const SPD& spd1, const SPD& spd2) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd1.phi[i] * spd2.phi[i];
  }
  return ret;
}
inline SPD operator/(const SPD& spd1, const SPD& spd2) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd1.phi[i] / spd2.phi[i];
  }
  return ret;
}

// SPDとRealの演算
inline SPD operator+(const SPD& spd, const Real& k) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] + k;
  }
  return ret;
}
inline SPD operator+(const Real& k, const SPD& spd) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] + k;
  }
  return ret;
}
inline SPD operator-(const SPD& spd, const Real& k) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] - k;
  }
  return ret;
}
inline SPD operator-(const Real& k, const SPD& spd) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = k - spd.phi[i];
  }
  return ret;
}
inline SPD operator*(const SPD& spd, const Real& k) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] * k;
  }
  return ret;
}
inline SPD operator*(const Real& k, const SPD& spd) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] * k;
  }
  return ret;
}
inline SPD operator/(const SPD& spd, const Real& k) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] / k;
  }
  return ret;
}
inline SPD operator/(const Real& k, const SPD& spd) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = k / spd.phi[i];
  }
  return ret;
}

// 正規化
inline SPD normalize(const SPD& spd) {
  const auto m = std::max_element(spd.phi.begin(), spd.phi.end());
  return spd / *m;
}

// SPDの出力
inline std::ostream& operator<<(std::ostream& stream, const SPD& spd) {
  stream << std::setw(12) << "lambda" << std::setw(12) << "phi" << std::endl;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    const Real lambda = SPD::LAMBDA_MIN + i * SPD::LAMBDA_INTERVAL;
    stream << std::setw(12) << lambda << std::setw(12) << spd.phi[i]
           << std::endl;
  }
  return stream;
}

// An RGB to Spectrum Conversion for Reflectances, Smits(2001)
SPD RGB2Spectrum(const RGB& rgb);

}  // namespace Prl2

#endif