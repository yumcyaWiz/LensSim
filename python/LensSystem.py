from typing import Tuple

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

    def object_focal_z(self):
        return self.lsys.object_focal_z

    def object_principal_z(self):
        return self.lsys.object_principal_z

    def image_focal_z(self):
        return self.lsys.image_focal_z

    def image_principal_z(self):
        return self.lsys.image_principal_z

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
        ax.set_xlim([min(z_list) - 0.3*length,
                     self.image_focal_z() + 0.3*length])
        ax.set_ylim([-1.1 * max_aperture_radius, 1.1*max_aperture_radius])
        ax.set_aspect('equal')
        ax.grid('on')
        plt.xlabel('$z \mathrm{[m]}$')
        plt.ylabel('$y \mathrm{[m]}$')

        return ax

    def optical_path_diagram(self, n_rays=10, theta=0, origin=None):
        # Plot Lenses
        ax = self.plot()

        # Plot Optical Path
        for i in range(n_rays):
            u = 2*(i + 0.5)/n_rays - 1
            h = self.lsys.elements[0].aperture_radius

            # make ray
            ray_direction = LensSim.Vec3(0, np.sin(theta), np.cos(theta))
            ray_origin = LensSim.Vec3(
                0, u*h, self.lsys.elements[0].z) - ray_direction
            if origin is not None:
                ray_origin = LensSim.Vec3(origin[0], origin[1], origin[2])
                ray_direction = LensSim.normalize(LensSim.Vec3(
                    0, u*h, self.lsys.elements[0].z) - ray_origin)

            # raytrace
            rays = self.lsys.raytracePath(
                LensSim.Ray(ray_origin, ray_direction))

            # Optical Path
            line_x = list(map(lambda x: x.origin[2], rays))
            line_y = list(map(lambda x: x.origin[1], rays))

            # add last path
            if len(line_x) == len(self.lsys.elements) + 1:
                # compute intersection with gaussian plane
                if rays[-1].direction[2] != 0:
                    t = -(rays[-1].origin[2] - self.image_focal_z()) / \
                        rays[-1].direction[2]
                    line_x.append(rays[-1].origin[2] + t*rays[-1].direction[2])
                    line_y.append(rays[-1].origin[1] + t*rays[-1].direction[1])

            # Plot
            # ax.scatter(line_x, line_y, c="lime", s=10)
            ax.plot(line_x, line_y, c="lime")

        # Plot Cardinal Points
        ax.scatter(self.object_focal_z(), 0, c="red")
        ax.scatter(self.object_principal_z(), 0, c="blue")
        ax.scatter(self.image_focal_z(), 0, c="red")
        ax.scatter(self.image_principal_z(), 0, c="blue")

        return ax

    def plot_exit_pupil(self, pFilm: Tuple[float, float] = (0, 0), n_grids: int = 100):
        fig, ax = plt.subplots()

        # compute exit pupil
        grid, extent = self.lsys.computeExitPupil(
            LensSim.Vec2(pFilm[0], pFilm[1]), n_grids)

        # Plot
        ax.imshow(grid, extent=extent, cmap="gray")
        ax.set_xlabel("$x \mathrm{[m]}$")
        ax.set_ylabel("$y \mathrm{[m]}$")

        return ax

    def spot_diagram(self, origin: np.array, n_grids: int):
        fig, ax = plt.subplots()

        # compute spot diagram
        spots = self.lsys.computeSpotDiagram(
            LensSim.Vec3(origin[0], origin[1], origin[2]),
            n_grids
        )

        # retrieve x, y
        x = list(map(lambda x: x[0], spots))
        y = list(map(lambda x: x[1], spots))

        # Plot
        ax.scatter(x, y, s=10)
        ax.set_xlabel("$x \mathrm{[m]}$")
        ax.set_ylabel("$y \mathrm{[m]}$")

        return ax
