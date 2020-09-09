import LensSim


class LensSystem:
    def __init__(self, filepath: str, width: int, height: int, width_length: float, height_length: float):
        self.lsys = LensSim.LensSystem(
            filepath, width, height, width_length=width_length, height_length=height_length)

    def focal_length(self):
        return self.lsys.image_focal_length
