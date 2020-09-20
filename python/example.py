import matplotlib.pyplot as plt

import LensSim
from LensSystem import LensSystem

if __name__ == "__main__":
    # load lens
    lsys = LensSystem("../data/dgauss50mm.json", width=512,
                      height=512, width_length=0.025, height_length=0.025)

    # print lens parameters
    print("EFL: {0}".format(lsys.effective_focal_length()))
    print("FFL: {0}".format(lsys.front_focal_length()))
    print("BFL: {0}".format(lsys.back_focal_length()))

    # plot optical path diagram
    lsys.optical_path_diagram()
    plt.show()
