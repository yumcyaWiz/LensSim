# LensSim

LensSim is a raytracer for Photographic Lens System.

LensSim can simulate image of Lens System with IBL(Image Based Lighting) by Montecarlo Light Transportation.

## Features

- [x] Ray Tracing
- [x] Focusing
- [x] Sampling Ray from Exit Pupil
- [x] Vignetting
- [x] Lens Flare
- [x] Chromatic Aberration
- [x] Paraxial Ray Tracing
- [x] Spectral Rendering with IBL
- [x] Python Binding
- [ ] Geometric Point Spread Function(PSF)
- [ ] Geometric Optical Transfer Function(OTF), Modulation Transfer Function(MTF)
- [ ] Wavefront Aberration
- [ ] Zernike Polynomial
- [ ] Diffraction Point Spread Function(PSF)
- [ ] Diffraction Optical Transfer Function(OTF), Modulation Transfer Function(MTF)

### Supported Lens Element

- [x] Aperture
- [x] Spherical Lens
- [ ] Aspheric Lens

## Requirements

* C++17
* CMake 3.12 or Higher

## Setup

```
git submodule update --init
```

## Build

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
## Python Binding

See `python/example.py`


## Externals

* [progschj/ThreadPool](https://github.com/progschj/ThreadPool) - zlib License
* [nlohmann/json](https://github.com/nlohmann/json) - MIT License
* [nothings/stb](https://github.com/nothings/stb) - Public Domain or MIT License 
* [syoyo/tinyexr](https://github.com/syoyo/tinyexr) - 3-clause BSD License
* [pybind/pybind11](https://github.com/pybind/pybind11) - BSD-style license 
