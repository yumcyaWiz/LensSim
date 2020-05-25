#ifndef _FILM_H
#define _FILM_H

#include <cmath>
#include <fstream>
#include <iostream>

// ext
#include "tinyexr.h"

// prl2
#include "core/vec3.h"

using namespace Prl2;

class Film {
 public:
  unsigned int width;
  unsigned int height;
  Real width_length;
  Real height_length;
  Real diagonal_length;

  Real* pixels;

  Film(unsigned int _width, unsigned int _height, Real _width_length = 0.036f,
       Real _height_length = 0.024f)
      : width(_width),
        height(_height),
        width_length(_width_length),
        height_length(_height_length) {
    pixels = new Real[width * height];
    diagonal_length =
        std::sqrt(width_length * width_length + height_length * height_length);

    // Init Pixels
    for(unsigned int j = 0; j < height; ++j) {
      for(unsigned int i = 0; i < width; ++i) {
        pixels[i + width*j] = 0.0f;
      }
    }
  }

  ~Film() { delete[] pixels; }

  Real getPixel(unsigned int i, unsigned int j) const {
    return pixels[i + width * j];
  }

  void setPixel(unsigned int i, unsigned int j, const Real& c) {
    pixels[i + width * j] = c;
  }

  void addPixel(unsigned int i, unsigned int j, Real radiance) {
    pixels[i + width * j] += radiance;
  }

  void divide(unsigned int k) {
    for (int j = 0; j < height; ++j) {
      for (int i = 0; i < width; ++i) {
        pixels[i + width * j] /= k;
      }
    }
  }

  Vec2 computePosition(Real u, Real v) const {
    return Vec2(0.5f * width_length * u, 0.5f * height_length * v);
  }

  void writePPM(const std::string& filename) const {
    std::ofstream file(filename);
    file << "P3" << std::endl;
    file << width << " " << height << std::endl;
    file << "255" << std::endl;

    for (int j = 0; j < height; ++j) {
      for (int i = 0; i < width; ++i) {
        const Vec3 rgb(getPixel(i, j));
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

  void writeEXR(const std::string& filename) const {
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
      const Vec3 rgb(pixels[i]);
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

  void writeCSV(const std::string& filename) const {
    std::ofstream file(filename);

    for (int j = 0; j < height; ++j) {
      for (int i = 0; i < width; ++i) {
        const Real r = getPixel(i, j);
        if (i < width - 1) {
          file << r << ",";
        } else {
          file << r;
        }
      }
      file << std::endl;
    }

    file.close();
  }
};

#endif
