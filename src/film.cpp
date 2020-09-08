#include "film.h"

#include <cmath>
#include <fstream>
#include <iostream>

#include "core/spectrum.h"

// ext
#include "tinyexr.h"

Film::Film(unsigned int _width, unsigned int _height, Real _width_length,
           Real _height_length)
    : width(_width),
      height(_height),
      width_length(_width_length),
      height_length(_height_length) {
  pixels = new Vec3[width * height];
  diagonal_length =
      std::sqrt(width_length * width_length + height_length * height_length);
}

Film::~Film() { delete[] pixels; }

Vec3 Film::getPixel(unsigned int i, unsigned int j) const {
  return pixels[i + width * j];
}

void Film::setPixel(unsigned int i, unsigned int j, const Vec3& c) {
  pixels[i + width * j] = c;
}

void Film::addPixel(unsigned int i, unsigned int j, Real lambda,
                    Real radiance) {
  //対応する等色関数のインデックスを計算
  const int index = (lambda - 380) / 5;

  // CMFを線形補間して計算
  Vec3 XYZ;
  if (index >= 0 && index <= SPD::color_matching_func_samples - 1) {
    XYZ[0] = radiance * SPD::color_matching_func_x[index];
    XYZ[1] = radiance * SPD::color_matching_func_y[index];
    XYZ[2] = radiance * SPD::color_matching_func_z[index];
  }

  // (i, j)に加算
  pixels[i + width * j] += XYZ;
}

void Film::divide(unsigned int k) {
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      pixels[i + width * j] /= k;
    }
  }
}

Vec2 Film::computePosition(Real u, Real v) const {
  return Vec2(0.5f * width_length * u, 0.5f * height_length * v);
}

void Film::writePPM(const std::string& filename) const {
  std::ofstream file(filename);
  file << "P3" << std::endl;
  file << width << " " << height << std::endl;
  file << "255" << std::endl;

  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      const Vec3 rgb = XYZ2RGB(getPixel(i, j));
      unsigned int r = std::min(
          static_cast<unsigned int>(255 * std::pow(rgb.x(), 1 / 2.2)), 255U);
      unsigned int g = std::min(
          static_cast<unsigned int>(255 * std::pow(rgb.y(), 1 / 2.2)), 255U);
      unsigned int b = std::min(
          static_cast<unsigned int>(255 * std::pow(rgb.z(), 1 / 2.2)), 255U);

      file << r << " " << g << " " << b << std::endl;
    }
  }

  file.close();
}

void Film::writeEXR(const std::string& filename) const {
  EXRHeader header;
  InitEXRHeader(&header);

  EXRImage image;
  InitEXRImage(&image);

  image.num_channels = 3;

  std::vector<float> images[3];
  images[0].resize(width * height);
  images[1].resize(width * height);
  images[2].resize(width * height);

  for (size_t i = 0; i < width * height; ++i) {
    const Vec3 rgb = XYZ2RGB(pixels[i]);
    images[0][i] = rgb.x();
    images[1][i] = rgb.y();
    images[2][i] = rgb.z();
  }

  float* image_ptr[3];
  image_ptr[0] = &(images[2].at(0));
  image_ptr[1] = &(images[1].at(0));
  image_ptr[2] = &(images[0].at(0));

  image.images = reinterpret_cast<unsigned char**>(image_ptr);
  image.width = width;
  image.height = height;

  header.num_channels = 3;
  header.channels = new EXRChannelInfo[header.num_channels];
  strncpy(header.channels[0].name, "B", 255);
  header.channels[0].name[strlen("B")] = '\0';
  strncpy(header.channels[1].name, "G", 255);
  header.channels[1].name[strlen("G")] = '\0';
  strncpy(header.channels[2].name, "R", 255);
  header.channels[2].name[strlen("R")] = '\0';

  header.pixel_types = new int[header.num_channels];
  header.requested_pixel_types = new int[header.num_channels];
  for (int i = 0; i < header.num_channels; ++i) {
    header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;
    header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF;
  }

  const char* err = nullptr;
  int ret = SaveEXRImageToFile(&image, &header, filename.c_str(), &err);
  if (ret != TINYEXR_SUCCESS) {
    fprintf(stderr, "Save EXR error: %s\n", err);
    FreeEXRErrorMessage(err);
    return;
  }
  printf("Saved EXR file. [%s] \n", filename.c_str());

  delete[] header.channels;
  delete[] header.pixel_types;
  delete[] header.requested_pixel_types;
}