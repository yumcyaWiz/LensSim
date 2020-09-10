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

    def vertical_fov(self):
        return self.lsys.vertical_fov

    def horizontal_fov(self):
        return self.lsys.horizontal_fov

    def diagonal_fov(self):
        return self.lsys.diagonal_fov

    def length(self):
        return self.lsys.system_length

    def plot(self):
        """
        plot LensSystem with matplotlib
        """
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
        ax.set_xlim([min(z_list) - 0.3*length, 0.3*length])
        ax.set_ylim([-1.1 * max_aperture_radius, 1.1*max_aperture_radius])
        ax.set_aspect('equal')
        ax.grid('on')
        plt.xlabel('$z \mathrm{[m]}$')
        plt.ylabel('$y \mathrm{[m]}$')

        return ax

    def optical_path_diagram(self, n_rays=10):
        # Plot Lenses
        ax = self.plot()

        # Plot Optical Path
        for i in range(n_rays):
            u = 2*(i + 0.5)/n_rays - 1
            h = self.lsys.elements[0].aperture_radius
            rays = self.lsys.raytrace_path(
                LensSim.Ray(
                    LensSim.Vec3(0, u*h, self.lsys.elements[0].z - 1),
                    LensSim.Vec3(0, 0, 1)
                )
            )

            # Optical Path
            line_x = list(map(lambda x: x.origin[2], rays))
            line_y = list(map(lambda x: x.origin[1], rays))

            # add last path
            if len(line_x) == len(self.lsys.elements) + 1:
                # compute intersection with y = 0
                if rays[-1].direction[1] != 0:
                    t = -rays[-1].origin[1]/rays[-1].direction[1]
                    line_x.append(rays[-1].origin[2] + t*rays[-1].direction[2])
                    line_y.append(rays[-1].origin[1] + t*rays[-1].direction[1])

            # Plot
            # ax.scatter(line_x, line_y, c="lime", s=10)
            ax.plot(line_x, line_y, c="lime")

        return ax
