import numpy as np
import matplotlib.pyplot as plt

from LensSim import Vec2, Vec3, Ray
from LensSystem import LensSystem

lsys = LensSystem("../data/dgauss50mm.json", width=512,
                  height=512, width_length=0.025, height_length=0.025)
print(lsys.focal_length())
print(np.degrees(lsys.vertical_fov()))

# lsys.optical_path_diagram()
# plt.show()

grid, extends = lsys.lsys.computeExitPupil(Vec2(0, 0), 100)
plt.imshow(np.array(grid), extent=extends, cmap='gray')
plt.show()
