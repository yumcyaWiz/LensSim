import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import LensSim


class LensSystem:
    def __init__(self, filepath: str, width: int, height: int, width_length: float, height_length: float):
        self.lsys = LensSim.LensSystem(
            filepath, width, height, width_length=width_length, height_length=height_length)

    def focal_length(self):
        return self.lsys.image_focal_length

    def plot(self):
        fig, ax = plt.subplots()
        for i in range(len(self.lsys.elements)):
            element = self.lsys.elements[i]
            # Draw Lens Element
            # Aperture
            if element.is_aperture:
                line_x = [element.z, element.z]
                line_y = [element.aperture_radius,
                          1.2*element.aperture_radius]
                ax.plot(line_x, line_y, c="blue")
                line_y = [-element.aperture_radius, -
                          1.2*element.aperture_radius]
                ax.plot(line_x, line_y, c="blue")
            # Spherical Lens
            else:
                z = element.z
                r = element.curvature_radius
                h = element.aperture_radius
                theta = abs(np.degrees(np.arcsin(h / r)))
                angle = 180 if r > 0 else 0
                arc = patches.Arc((z + r, 0), 2*abs(r), 2*abs(r),
                                  angle=angle, theta1=-theta, theta2=theta)
                ax.add_patch(arc)

            # Draw Lens Box
            if i > 0:
                element_prev = self.lsys.elements[i - 1]

                # current or previous element is aperture radius
                if element.is_aperture or element_prev.is_aperture:
                    continue

                # previous element is air
                if element_prev.eta == 1:
                    continue

                z = element.z
                r = element.curvature_radius
                h = element.aperture_radius
                l = r - \
                    np.sqrt(r**2 - h**2) if r > 0 else -(np.abs(r) -
                                                         np.sqrt(r**2 - h**2))
                zp = element_prev.z
                rp = element_prev.curvature_radius
                hp = element_prev.aperture_radius
                lp = rp - \
                    np.sqrt(rp**2 - hp**2) if rp > 0 else -(np.abs(rp) -
                                                            np.sqrt(rp**2 - hp**2))

                if h > hp:
                    ax.plot([zp + lp, z + l], [h, h], c="black")
                    ax.plot([zp + lp, z + l], [-h, -h], c="black")
                    ax.plot([zp + lp, zp + lp], [hp, h], c="black")
                    ax.plot([zp + lp, zp + lp], [-hp, -h], c="black")
                else:
                    ax.plot([zp + lp, z + l], [hp, hp], c="black")
                    ax.plot([zp + lp, z + l], [-hp, -hp], c="black")
                    ax.plot([z + l, z + l], [h, hp], c="black")
                    ax.plot([z + l, z + l], [-h, -hp], c="black")

        # figure width, height
        z_list = [element.z for element in self.lsys.elements]
        length = max(z_list) - min(z_list)
        max_aperture_radius = max([
            element.aperture_radius for element in self.lsys.elements])
        ax.set_xlim([min(z_list) - 0.3*length, max(z_list) + 0.3*length])
        ax.set_ylim([-1.1 * max_aperture_radius, 1.1*max_aperture_radius])
        ax.set_aspect('equal')
        ax.grid('on')
        plt.xlabel('$z \mathrm{[mm]}$')
        plt.ylabel('$y \mathrm{[mm]}$')
