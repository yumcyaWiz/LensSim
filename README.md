# LensSim

LensSim is a raytracer for Photographic Lens System.

LensSim can simulate image of Lens System with IBL(Image Based Lighting) by Montecarlo Light Transportation.

## Features

- [x] Ray Tracing
- [x] Focusing
- [x] Exit Pupil Sampling
- [x] Spectral Rendering with IBL
- [x] Vignetting
- [x] Lens Flare
- [x] Chromatic Aberration

### Lens Elements

* Aperture
* Spherical Lens

## Build

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Externals

* [progschj/ThreadPool](https://github.com/progschj/ThreadPool) - zlib License
* [nlohmann/json](https://github.com/nlohmann/json) - MIT License
* [nothings/stb](https://github.com/nothings/stb) - Public Domain or MIT License 
* [syoyo/tinyexr](https://github.com/syoyo/tinyexr) - 3-clause BSD License
