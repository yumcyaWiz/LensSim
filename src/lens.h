#ifndef _LENS_H
#define _LENS_H

class Lens {
 public:
  unsigned int index;
  double curvature_radius;
  double aperture_radius;
  double thickness;
  double ior;
  double z;

  Lens(unsigned int _index, double _curvature_radius, double _thickness,
       double _ior)
      : index(_index),
        aperture_radius(_curvature_radius),
        thickness(_thickness),
        ior(_ior),
        z(0){};
};

#endif