#ifndef _FILM_H
#define _FILM_H

#include <cmath>
#include <fstream>
#include <iostream>

// ext
#include "tinyexr.h"
//
#include "tiny_dng_loader.h"
//
#include "examples/dngwriter/tiny_dng_writer.h"

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

  Vec3* pixels;
  unsigned int* samples;

  Film(unsigned int _width, unsigned int _height, Real _width_length = 0.036f,
       Real _height_length = 0.024f)
      : width(_width),
        height(_height),
        width_length(_width_length),
        height_length(_height_length) {
    pixels = new Vec3[width * height];
    samples = new unsigned int[width * height];
    diagonal_length =
        std::sqrt(width_length * width_length + height_length * height_length);

    for (int j = 0; j < height; ++j) {
      for (int i = 0; i < width; ++i) {
        samples[i + width * j] = 0;
      }
    }
  }

  ~Film() { delete[] pixels; }

  Vec3 getPixel(unsigned int i, unsigned int j) const {
    unsigned int s = samples[i + width * j];
    if (s > 0) {
      return pixels[i + width * j] / s;
    } else {
      return Vec3(0, 0, 0);
    }
  }

  void setPixel(unsigned int i, unsigned int j, const Vec3& c) {
    pixels[i + width * j] = c;
  }

  void addPixel(unsigned int i, unsigned int j, const Vec3& radiance) {
    pixels[i + width * j] += radiance;
    samples[i + width * j] += 1;
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
        const Vec3 rgb = getPixel(i, j);
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
      unsigned int s = samples[i];
      Vec3 rgb;
      if (s > 0) {
        rgb = pixels[i] / s;
      }
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
        const Vec3 r = getPixel(i, j);
        if (i < width - 1) {
          file << r.x() << ",";
        } else {
          file << r.x();
        }
      }
      file << std::endl;
    }

    file.close();
  }

  // void writeTIFF(const std::string& filename) const {
  //   std::vector<float> image(3 * width * height);
  //   for (int j = 0; j < height; ++j) {
  //     for (int i = 0; i < width; ++i) {
  //       unsigned int s = samples[i];
  //       Vec3 rgb;
  //       if (s > 0) {
  //         rgb = pixels[i] / s;
  //       }

  //       image[3 * i + 3 * width * j + 0] = rgb.x();
  //       image[3 * i + 3 * width * j + 1] = rgb.y();
  //       image[3 * i + 3 * width * j + 2] = rgb.z();
  //     }
  //   }

  //   tinydngwriter::DNGImage tiff_image;
  //   tiff_image.SetSubfileType(tinydngwriter::FILETYPE_REDUCEDIMAGE);
  //   tiff_image.SetImageWidth(width);
  //   tiff_image.SetImageLength(height);
  //   tiff_image.SetRowsPerStrip(height);
  //   uint16_t bps = 32;
  //   tiff_image.SetBitsPerSample(1, &bps);
  //   tiff_image.SetPlanarConfig(tinydngwriter::PLANARCONFIG_CONTIG);
  //   tiff_image.SetCompression(tinydngwriter::COMPRESSION_NONE);
  //   tiff_image.SetPhotometric(tinydngwriter::PHOTOMETRIC_RGB);
  //   tiff_image.SetSamplesPerPixel(1);

  //   uint16_t format = tinydngwriter::SAMPLEFORMAT_IEEEFP;
  //   tiff_image.SetSampleFormat(1, &format);

  //   tiff_image.SetImageData(reinterpret_cast<unsigned char*>(image.data()),
  //                           image.size() * sizeof(float));

  //   tinydngwriter::DNGWriter tiff_writer(false);
  //   tiff_writer.AddImage(&tiff_image);

  //   std::string err;
  //   bool ret = tiff_writer.WriteToFile(filename.c_str(), &err);

  //   if (!err.empty()) {
  //     std::cerr << err << std::endl;
  //   }
  // }
};

#endif
