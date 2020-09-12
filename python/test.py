import numpy as np
import matplotlib.pyplot as plt

from LensSim import Vec2, Vec3, Ray
from LensSystem import LensSystem

lsys = LensSystem("../data/dgauss50mm.json", width=512,
                  height=512, width_length=0.025, height_length=0.025)
print(lsys.focal_length())
print(np.degrees(lsys.vertical_fov()))
print(lsys.lsys.entrance_pupil_z)
print(lsys.lsys.exit_pupil_z)
lsys.optical_path_diagram()
plt.show()

# primary_ray = Ray()
# lsys.lsys.computePrimaryRay(Vec3(0, 0, -1), primary_ray)
# print(primary_ray)
