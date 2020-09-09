import numpy as np
import matplotlib.pyplot as plt

from LensSim import LensSystem, Vec3, Ray

lsys = LensSystem(
    filename="../data/dgauss50mm.json", width=512, height=512)
print(lsys.object_focal_length)
print(lsys.image_focal_length)
print(lsys.elements[0].aperture_radius)
