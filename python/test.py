import matplotlib.pyplot as plt

from LensSim import LensSystem

lsys = LensSystem(filename="../data/dgauss50mm.json", width=512, height=512)
print(lsys.object_focal_length)
print(lsys.image_focal_length)
