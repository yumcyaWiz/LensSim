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
  py::class_<Vec2>(m, "Vec2", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Real>())
      .def(py::init<Real, Real>())

      .def_buffer([](Vec2& v) -> py::buffer_info {
        return py::buffer_info(v.v, sizeof(Real),
                               py::format_descriptor<Real>::format(), 1, {2},
                               {sizeof(Real)});
      })

      .def("__getitem__", [](const Vec2& v, int i) { return v.v[i]; })
      .def("__repr__", [](const Vec2& v) {
        return "(" + std::to_string(v.x()) + ", " + std::to_string(v.y()) + ")";
      });

  py::class_<Vec3>(m, "Vec3", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Real>())
      .def(py::init<Real, Real, Real>())

      .def_buffer([](Vec3& v) -> py::buffer_info {
        return py::buffer_info(v.v, sizeof(Real),
                               py::format_descriptor<Real>::format(), 1, {3},
                               {sizeof(Real)});
      })

      .def("__getitem__", [](const Vec3& v, int i) { return v.v[i]; })
      .def("__repr__", &vec2string)
      .def("__add__", [](const Vec3& v1, const Vec3& v2) { return v1 + v2; })
      .def("__sub__", [](const Vec3& v1, const Vec3& v2) { return v1 - v2; });

  py::class_<Ray>(m, "Ray")
      .def(py::init<>())
      .def(py::init<Vec3, Vec3, Real>(), py::arg("origin"),
           py::arg("direction"), py::arg("lambda") = 550.0)

      .def_readonly("origin", &Ray::origin)
      .def_readonly("direction", &Ray::direction)

      .def("__repr__", &ray2string);

  py::class_<Sampler>(m, "Sampler");

  py::class_<GridData<Real>>(m, "GridData", py::buffer_protocol())
      .def_buffer([](GridData<Real>& grid) -> py::buffer_info {
        return py::buffer_info(grid.data.data(), sizeof(Real),
                               py::format_descriptor<Real>::format(), 2,
                               {grid.nrows, grid.ncols},
                               {sizeof(Real) * grid.ncols, sizeof(Real)});
      });

  py::class_<ParaxialRay>(m, "ParaxialRay")
      .def(py::init<>())
      .def(py::init<Real, Real>())

      .def_readonly("u", &ParaxialRay::u)
      .def_readonly("h", &ParaxialRay::h)

      .def("__repr__", [](const ParaxialRay& ray) {
        return "(" + std::to_string(ray.u) + ", " + std::to_string(ray.h) + ")";
      });

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
      .def_readonly("object_focal_z", &LensSystem::object_focal_z)
      .def_readonly("object_principal_z", &LensSystem::object_principal_z)
      .def_readonly("image_focal_length", &LensSystem::image_focal_length)
      .def_readonly("image_focal_z", &LensSystem::image_focal_z)
      .def_readonly("image_principal_z", &LensSystem::image_principal_z)
      .def_readonly("elements", &LensSystem::elements)
      .def_readonly("system_length", &LensSystem::system_length)
      .def_readonly("vertical_fov", &LensSystem::vertical_fov)
      .def_readonly("horizontal_fov", &LensSystem::horizontal_fov)
      .def_readonly("diagonal_fov", &LensSystem::diagonal_fov)

      .def("raytrace", &LensSystem::raytrace, py::arg("ray_in"),
           py::arg("ray_out"), py::arg("reflection") = false,
           py::arg("sampler") = nullptr)
      .def("raytracePath", &LensSystem::raytracePath, py::arg("ray_in"))
      .def("raytraceParaxial", &LensSystem::raytraceParaxial, py::arg("ray_in"),
           py::arg("start") = 0, py ::arg("end") = -1,
           py::arg("lambda") = 550.0)
      .def("computeExitPupil", &LensSystem::computeExitPupil, py::arg("pFilm"),
           py::arg("n_grids") = 512)
      .def("computePrimaryRay", &LensSystem::computePrimaryRay,
           py::arg("origin"), py::arg("primary_ray"), py::arg("n_grids") = 512)
      .def("computeSpotDiagram", &LensSystem::computeSpotDiagram,
           py::arg("origin"), py::arg("n_grids"));
}