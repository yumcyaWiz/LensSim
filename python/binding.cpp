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

Vec3 vec3test() { return Vec3(1, 2, 3); }
Vec3 vec3test2(const Vec3& v1, const Vec3& v2) { return v1 + v2; }

PYBIND11_MODULE(LensSim, m) {
  py::class_<LensSystem>(m, "LensSystem")
      .def(py::init(&lsysFactory), py::arg("filename"), py::arg("width"),
           py::arg("height"), py::arg("width_length") = 0.025,
           py::arg("height_length") = 0.025)
      .def_readonly("object_focal_length", &LensSystem::object_focal_length)
      .def_readonly("image_focal_length", &LensSystem::image_focal_length);

  py::class_<Vec3>(m, "Vec3", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Real>())
      .def(py::init<Real, Real, Real>())
      .def_buffer([](Vec3& v) -> py::buffer_info {
        return py::buffer_info(v.v, sizeof(Real),
                               py::format_descriptor<Real>::format(), 1, {3},
                               {sizeof(Real)});
      });

  m.def("vec3test", &vec3test, "");
  m.def("vec3test2", &vec3test2, "");
}