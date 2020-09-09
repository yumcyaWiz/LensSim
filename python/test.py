import numpy as np
import matplotlib.pyplot as plt

import LensSim
from LensSim import Vec3, Ray

# lsys = LensSim.LensSystem(filename="../data/dgauss50mm.json", width=512, height=512)
# print(lsys.object_focal_length)
# print(lsys.image_focal_length)

print(Ray(Vec3(0, 0, 0), Vec3(0, 0, 1)))
