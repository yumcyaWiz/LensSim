import numpy as np
import matplotlib.pyplot as plt

# from LensSim import LensSystem, Vec3, Ray

# lsys = LensSystem(
#     filename="../data/dgauss50mm.json", width=512, height=512)
# print(lsys.object_focal_length)
# print(lsys.image_focal_length)
# print(lsys.elements[0].aperture_radius)

from LensSystem import LensSystem

lsys = LensSystem("../data/dgauss50mm.json", width=512,
                  height=512, width_length=0.025, height_length=0.025)

print(lsys.focal_length())

lsys.plot()
plt.show()
