import matplotlib.pyplot as plt

import LensSim
from LensSystem import LensSystem

if __name__ == "__main__":
    # load lens
    lsys = LensSystem("../data/dgauss50mm.json", width=512,
                      height=512, width_length=0.025, height_length=0.025)

    # print lens parameters
    print("Focal Length: {0}".format(lsys.focal_length()))

    # plot optical path diagram
    lsys.optical_path_diagram()
    plt.show()
