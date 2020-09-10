#include <memory>
#include <string>

#include "lens-system/lens-system.h"

// pybind11
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;

std::unique_ptr<LensSystem> lsysFactory(const std::string& filename,
                                        unsigned int width, unsigned int height,
                                        Real width_length = 0.025,
                                        Real height_length = 0.025) {
  const std::shared_ptr<Film> film =
      std::make_shared<Film>(width, height, width_length, height_length);
  return std::make_unique<LensSystem>(filename, film);
}

std::string vec2string(const Vec3& v) {
  return "(" + std::to_string(v.x()) + ", " + std::to_string(v.y()) + ", " +
         std::to_string(v.z()) + ")";
}
std::string ray2string(const Ray& ray) {
  return "origin: " + vec2string(ray.origin) + " " +
         "direction: " + vec2string(ray.direction);
}

PYBIND11_MODULE(LensSim, m) {
  py::class_<Vec3>(m, "Vec3", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Real>())
      .def(py::init<Real, Real, Real>())

      .def_buffer([](Vec3& v) -> py::buffer_info {
        return py::buffer_info(v.v, sizeof(Real),
                               py::format_descriptor<Real>::format(), 1, {3},
                               {sizeof(Real)});
      })

      .def("__repr__", &vec2string);

  py::class_<Ray>(m, "Ray")
      .def(py::init<>())
      .def(py::init<Vec3, Vec3, Real>(), py::arg("origin"),
           py::arg("direction"), py::arg("lambda") = 550.0)

      .def_readonly("origin", &Ray::origin)
      .def_readonly("direction", &Ray::direction)

      .def("__repr__", &ray2string);

  py::class_<Sampler>(m, "Sampler");

  py::class_<LensElement>(m, "LensElement")
      .def_readonly("curvature_radius", &LensElement::curvature_radius)
      .def_readonly("aperture_radius", &LensElement::aperture_radius)
      .def_readonly("thickness", &LensElement::thickness)
      .def_readonly("eta", &LensElement::eta)
      .def_readonly("z", &LensElement::z)
      .def_readonly("is_aperture", &LensElement::is_aperture);

  py::class_<LensSystem>(m, "LensSystem")
      .def(py::init(&lsysFactory), py::arg("filename"), py::arg("width"),
           py::arg("height"), py::arg("width_length") = 0.025,
           py::arg("height_length") = 0.025)

      .def_readonly("object_focal_length", &LensSystem::object_focal_length)
      .def_readonly("image_focal_length", &LensSystem::image_focal_length)
      .def_readonly("elements", &LensSystem::elements)

      .def("raytrace", &LensSystem::raytrace, py::arg("ray_in"),
           py::arg("ray_out"), py::arg("reflection") = false,
           py::arg("sampler") = nullptr);
}