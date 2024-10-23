import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import matplotlib
import pyFAI
from pathlib import Path
import pyFAI.azimuthalIntegrator
import pyFAI.detectors
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


class Nudger:

    def __init__(self, ai: pyFAI.AzimuthalIntegrator, data, flat_field=None, incident_angle=0, radii=[12, 5], IM_SIZE=(6.3, 3)):
        self.ai = ai
        self.incident_angle = np.radians(float(incident_angle))
        self.tilt_angle = 0
        self.nudge1 = 0
        self.nudge2 = 0
        self.radii = radii
        if flat_field is None:
            self.data = data
        else:
            self.data = data / flat_field
        self.fig = plt.figure(figsize=IM_SIZE)
        self.ax = None
        self.shape = data.shape
        self.show()

    def set_nudge(self, nudge1, nudge2):
        self.nudge1 = nudge1
        self.nudge2 = nudge2

    def nudge(self, nudge1, nudge2):
        self.nudge1 += nudge1
        self.nudge2 += nudge2

    def reset_nudge(self):
        self.nudge1 = 0
        self.nudge2 = 0

    def get_poni(self):
        poni1 = self.ai.poni1 + self.ai.pixel1 * self.nudge1
        poni2 = self.ai.poni2 + self.ai.pixel1 * self.nudge2
        print(f"poni1: {poni1}")
        print(f"poni2: {poni2}")
        return poni1, poni2

    def set_incident(self, incident_angle):
        self.incident_angle = np.radians(incident_angle)

    def set_tilt(self, tilt_angle):
        self.tilt_angle = np.radians(tilt_angle)
    
    def get_incident(self):
        return np.rad2deg(self.incident_angle)
    
    def get_tilt(self):
        return np.rad2deg(self.tilt)
    
    def show(self, display_max=None):
        self.fig.clf()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor("k")
        if display_max is None:
            display_max = self.data.max()
        pos = self.ax.imshow(self.data, norm=LogNorm(1, display_max))
        self.ax.set_title("Find PONI")
        self.ax.set_xlabel("column (pixels)")
        self.ax.set_ylabel("row (pixels)")
        self.ax.set_xlim(0, self.data.shape[1])
        self.ax.set_ylim(self.data.shape[0], 0)
        self.fig.tight_layout()

        poni1, poni2 = self.get_poni()

        x_pos = poni2 / self.ai.pixel2 - 0.5
        y_pos = self.data.shape[0] - poni1 / self.ai.pixel1 + 0.5
        radii = self.ai.dist / self.ai.pixel1 * np.tan(2. * np.arcsin(0.5 * self.ai.wavelength / np.array(self.radii)))

        self.ax.scatter(x_pos, y_pos, s=30, color="r")
        for r in radii:
            x = np.linspace(x_pos - r, x_pos + r, 1000)
            y = np.sqrt(r * r - (x - x_pos) ** 2)
            self.ax.plot(x, y + y_pos, "r", linewidth=0.5)
            self.ax.plot(x, -y + y_pos, "r", linewidth=0.5)
        y = np.linspace(0, self.data.shape[0], 100)
        vertical = np.tan(self.tilt_angle) * (y - y_pos) + x_pos
        self.ax.plot(vertical, y, "r", linewidth=0.5)

        x = np.linspace(0, self.data.shape[1], 100)
        horizontal = -np.tan(self.tilt_angle) * (x - x_pos) + y_pos
        horizon = horizontal - self.ai.dist / self.ai.pixel1 * np.tan(self.incident_angle)
        self.ax.plot(x, horizontal, "r", linewidth=0.5)
        self.ax.plot(x, horizon, "r", linewidth=0.5)

    def save(self, poni_name: Path, orientation: int = 2):
        self.ai.poni1 = self.ai.poni1 + self.nudge1 * self.ai.pixel1
        self.ai.poni2 = self.ai.poni2 + self.nudge2 * self.ai.pixel2
        self.ai.detector = pyFAI.detectors.Detector(pixel1=self.ai.pixel1, pixel2=self.ai.pixel2,
                                               max_shape=self.ai.detector.shape, orientation=orientation)
        print("Saving geometry:")
        print(self.ai)
        self.ai.save(poni_name)


def new(dist: float, poni1: float, poni2: float, shape: tuple, pixel1: float = 75e-6, pixel2: float = 75e-6,
        wavelength: float = 1.54185e-10, rot1: float = 0, rot2: float = 0, rot3: float = 0, orientation: int = 2):
    detector = pyFAI.detectors.Detector(pixel1=float(pixel1), pixel2=float(pixel2), max_shape=shape, orientation=orientation)
    ai = pyFAI.azimuthalIntegrator.AzimuthalIntegrator(float(dist), float(poni1), float(poni2), rot1, rot2, rot3, pixel1, pixel2,
                                                       detector=detector, wavelength=wavelength, orientation=orientation)
    return ai


def save_new(poni_name: Path, dist: float, poni1: float, poni2: float,
             shape: tuple, pixel1: float = 75e-6, pixel2: float = 75e-6,
             wavelength: float = 1.54185, rot1: float = 0, rot2: float = 0, rot3: float = 0, orientation: int = 2):
    ai = new(dist, poni1, poni2, shape, pixel1, pixel2, wavelength, rot1, rot2, rot3, orientation)
    ai.save(poni_name)
