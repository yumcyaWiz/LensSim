#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>

#include "core/spectrum.h"

namespace Prl2 {

//非等間隔なSPDから等間隔なSPDを構成する
//非等間隔なSPDを線形補間して計算する
SPD::SPD(const std::vector<Real>& _lambda, const std::vector<Real>& _phi) {
  assert(_lambda.size() == _phi.size());

  //サンプル数が1つの場合
  if (_lambda.size() == 1) {
    //他の波長を0で初期化
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] = 0;
    }
    addPhi(_lambda[0], _phi[0]);
  }
  //複数のサンプルがある場合
  else {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      //等間隔側の波長
      const Real lambda_value = LAMBDA_MIN + LAMBDA_INTERVAL * i;

      //指定した波長が非等間隔のSPDの範囲に含まれていない場合、放射束を0とする
      if (lambda_value < _lambda.front() || lambda_value > _lambda.back()) {
        phi[i] = 0;
      } else if (lambda_value == _lambda.front()) {
        phi[i] = _phi.front();
      } else if (lambda_value == _lambda.back()) {
        phi[i] = _phi.back();
      } else {
        //非等間隔のSPDを線形補間
        const size_t lambda1_index = static_cast<size_t>(
            std::lower_bound(_lambda.begin(), _lambda.end(), lambda_value) -
            _lambda.begin());
        const size_t lambda0_index = lambda1_index - 1;
        assert(lambda0_index != lambda1_index);

        const Real t = (lambda_value - _lambda[lambda0_index]) /
                       (_lambda[lambda1_index] - _lambda[lambda0_index]);
        assert(t >= 0 && t <= 1);

        const Real interpolated_phi =
            (1.0f - t) * _phi[lambda0_index] + t * _phi[lambda1_index];
        phi[i] = interpolated_phi;
      }
    }
  }
}

void SPD::addPhi(const Real& _lambda, const Real& _phi) {
  //範囲外の寄与は加算しない
  if (_lambda < LAMBDA_MIN || _lambda >= LAMBDA_MAX) {
    return;
  }

  //対応する波長のインデックスを計算
  const size_t lambda_index = (_lambda - LAMBDA_MIN) / LAMBDA_INTERVAL;

  //両側に寄与を分配する
  const Real lambda0 =
      LAMBDA_MIN + lambda_index * LAMBDA_INTERVAL;  //左側の波長
  const Real t =
      (_lambda - lambda0) / LAMBDA_INTERVAL;  //与えられた波長の区間中の位置
  phi[lambda_index] += (1.0f - t) * _phi;
  phi[lambda_index + 1] += t * _phi;
}

Real SPD::sample(const Real& l) const {
  assert(l >= LAMBDA_MIN && l < LAMBDA_MAX);

  //対応する波長のインデックスを計算
  const size_t lambda_index = (l - LAMBDA_MIN) / LAMBDA_INTERVAL;

  //放射束を線形補間して計算
  if (lambda_index == 0 || lambda_index == LAMBDA_SAMPLES - 1) {
    return phi[lambda_index];
  } else {
    const Real lambda_nearest = LAMBDA_MIN + lambda_index * LAMBDA_INTERVAL;
    const Real t = (l - lambda_nearest) / LAMBDA_INTERVAL;
    assert(t >= 0 && t <= 1);
    return (1.0f - t) * phi[lambda_index] + t * phi[lambda_index + 1];
  }
}

//等色関数は線形補間して使用する
XYZ SPD::toXYZ() const {
  XYZ xyz;

  for (std::size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
    const Real lambda_value = LAMBDA_MIN + LAMBDA_INTERVAL * i;
    const Real phi_value = phi[i];

    //放射束が0の場合は計算をスキップ
    if (phi_value == 0) continue;

    //対応する等色関数のインデックスを計算
    const int index = (lambda_value - 380) / 5;
    assert(index >= 0 && index < color_matching_func_samples);

    //等色関数を線形補間
    Real cmf_x, cmf_y, cmf_z;
    if (index != color_matching_func_samples - 1) {
      const Real cmf_lambda = 5 * index + 380;
      const Real t = (lambda_value - cmf_lambda) / 5;
      assert(t >= 0 && t <= 1);

      cmf_x = (1.0f - t) * color_matching_func_x[index] +
              t * color_matching_func_x[index + 1];
      cmf_y = (1.0f - t) * color_matching_func_y[index] +
              t * color_matching_func_y[index + 1];
      cmf_z = (1.0f - t) * color_matching_func_z[index] +
              t * color_matching_func_z[index + 1];
    } else {
      cmf_x = color_matching_func_x[index];
      cmf_y = color_matching_func_y[index];
      cmf_z = color_matching_func_z[index];
    }

    // XYZを計算(短冊近似)
    xyz[0] += cmf_x * phi_value;
    xyz[1] += cmf_y * phi_value;
    xyz[2] += cmf_z * phi_value;
  }

  return xyz;
}

SPD RGB2Spectrum(const RGB& rgb) {
  static const std::vector<Real> sampled_lambda = {
      380, 417.7, 455.55, 493.33, 531.11, 568.88, 606.66, 644.44, 682.22, 720};

  static const SPD white_spectrum =
      SPD(sampled_lambda, {1, 1, 0.9999, 0.9993, 0.9992, 0.9998, 1, 1, 1, 1});

  static const SPD cyan_spectrum =
      SPD(sampled_lambda,
          {0.9710, 0.9426, 1.0007, 1.0007, 1.0007, 1.0007, 0.1564, 0, 0, 0});

  static const SPD magenta_spectrum = SPD(
      sampled_lambda, {1, 1, 0.9685, 0.2229, 0, 0.0458, 0.8369, 1, 1, 0.9959});

  static const SPD yellow_spectrum =
      SPD(sampled_lambda,
          {0.0001, 0, 0.1088, 0.6651, 1, 1, 0.9996, 0.9586, 0.9685, 0.9840});

  static const SPD red_spectrum =
      SPD(sampled_lambda,
          {0.1012, 0.0515, 0, 0, 0, 0, 0.8325, 1.0149, 1.0149, 1.0149});

  static const SPD green_spectrum = SPD(
      sampled_lambda, {0, 0, 0.0273, 0.7937, 1, 0.9418, 0.1719, 0, 0, 0.0025});

  static const SPD blue_spectrum =
      SPD(sampled_lambda,
          {1, 1, 0.8916, 0.3323, 0, 0, 0.0003, 0.0369, 0.0483, 0.0496});

  SPD ret;
  const Real red = rgb.x();
  const Real green = rgb.y();
  const Real blue = rgb.z();
  if (red <= green && red <= blue) {
    ret += red * white_spectrum;
    if (green <= blue) {
      ret += (green - red) * cyan_spectrum;
      ret += (blue - green) * blue_spectrum;
    } else {
      ret += (blue - red) * cyan_spectrum;
      ret += (green - blue) * green_spectrum;
    }
  } else if (green <= red && green <= blue) {
    ret += green * white_spectrum;
    if (red <= blue) {
      ret += (red - green) * magenta_spectrum;
      ret += (blue - red) * blue_spectrum;
    } else {
      ret += (blue - green) * magenta_spectrum;
      ret += (red - blue) * red_spectrum;
    }
  } else {
    ret += blue * white_spectrum;
    if (red <= green) {
      ret += (red - blue) * yellow_spectrum;
      ret += (green - red) * green_spectrum;
    } else {
      ret += (green - blue) * yellow_spectrum;
      ret += (red - green) * red_spectrum;
    }
  }

  return ret;
}

}  // namespace Prl2