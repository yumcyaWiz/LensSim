#include <memory>
#include <string>

#include "lens-system/lens-system.h"
#include "pybind11/pybind11.h"

namespace py = pybind11;

std::unique_ptr<LensSystem> lsysFactory(const std::string& filename,
                                        unsigned int width, unsigned int height,
                                        Real width_length = 0.025,
                                        Real height_length = 0.025) {
  const std::shared_ptr<Film> film =
      std::make_shared<Film>(width, height, width_length, height_length);
  return std::make_unique<LensSystem>(filename, film);
}

PYBIND11_MODULE(LensSim, m) {
  py::class_<LensSystem>(m, "LensSystem")
      .def(py::init(&lsysFactory), py::arg("filename"), py::arg("width"),
           py::arg("height"), py::arg("width_length") = 0.025,
           py::arg("height_length") = 0.025)
      .def_readonly("object_focal_length", &LensSystem::object_focal_length)
      .def_readonly("image_focal_length", &LensSystem::image_focal_length);
}