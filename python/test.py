import numpy as np
import matplotlib.pyplot as plt

# from LensSim import LensSystem, Vec3, Ray

# lsys = LensSystem(
#     filename="../data/dgauss50mm.json", width=512, height=512)
# print(lsys.object_focal_length)
# print(lsys.image_focal_length)
# print(lsys.elements[0].aperture_radius)

from LensSim import Vec3, Ray
from LensSystem import LensSystem

lsys = LensSystem("../data/dgauss50mm.json", width=512,
                  height=512, width_length=0.025, height_length=0.025)

ray_in = Ray(Vec3(0, 0, -1), Vec3(0, 0, 1))
ray_out = Ray()
lsys.lsys.raytrace(ray_in=ray_in,
                   ray_out=ray_out)
print(ray_out)
