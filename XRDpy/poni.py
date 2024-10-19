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

def show(data, poni1, poni2, det_dist, incident_angle=0, tilt_angle=0, nudge1=0, nudge2=0,
         flat_field=None, radii=[12, 5], pixel_size=75e-3, wavelength=1.54185, IM_SIZE=(6.3, 3)):
    """
    Display an image and overlay PONI (Point of Normal Incidence) markers for manual adjustment.
    Parameters:
    data (ndarray): The 2D array of image data.
    poni1 (float): The first PONI coordinate.
    poni2 (float): The second PONI coordinate.
    det_dist (float): The detector distance.
    incident_angle (float, optional): The incident angle in degrees. Default is 0.
    tilt_angle (float, optional): The tilt angle in degrees. Default is 0.
    nudge1 (float, optional): Adjustment for the first PONI coordinate in pixels. Default is 0.
    nudge2 (float, optional): Adjustment for the second PONI coordinate in pixels. Default is 0.
    flat_field (ndarray, optional): Flat field correction array. Default is None.
    radii (list, optional): List of circle radii in A for lattice spacings.
    pixel_size (float, optional): The size of a pixel in mm. Default is 75e-3.
    wavelength (float, optional): Wavelength in A.
    IM_SIZE (tuple, optional): The size of the figure. Default is (6.3, 3).
    Returns:
    fig, ax
    """
    fig, ax = plt.subplots(1, 1, figsize=IM_SIZE, facecolor="w")
    ax.set_facecolor("k")
    if flat_field is not None:
        data /= flat_field
    pos = ax.imshow(data, norm=LogNorm(1, np.max(data)))
    ax.set_title("Find PONI")
    ax.set_xlabel("column (pixels)")
    ax.set_ylabel("row (pixels)")
    ax.set_xlim(0, data.shape[1])
    ax.set_ylim(0, data.shape[0])
    # fig.colorbar(pos, ax=ax, shrink=0.7, label="counts")
    fig.tight_layout()

    det_dist /= pixel_size
    poni1 += pixel_size * nudge1
    poni2 += pixel_size * nudge2

    x_pos = poni2 / pixel_size - 0.5
    y_pos = data.shape[0] - poni1 / pixel_size + 0.5
    radii = det_dist * np.tan(2. * np.arcsin(0.5 * wavelength / np.radii(radii)))

    tilt_angle = np.radians(tilt_angle)
    incident_angle = np.radians(incident_angle)

    ax.scatter(x_pos, y_pos, s=30, color="r")
    for r in radii:
        x = np.linspace(x_pos - r, x_pos + r, 1000)
        y = np.sqrt(r * r - (x - x_pos) ** 2)
        ax.plot(x, y + y_pos, "r", linewidth=0.5)
        ax.plot(x, -y + y_pos, "r", linewidth=0.5)
    y = np.linspace(0, data.shape[0], 100)
    vertical = np.tan(tilt_angle) * (y - y_pos) + x_pos
    x = np.linspace(0, data.shape[1], 100)
    horizontal = -np.tan(tilt_angle) * (x - x_pos) + y_pos
    horizon = horizontal - det_dist * np.tan(incident_angle)

    ax.plot(x, horizontal, "r", linewidth=0.5)
    ax.plot(x, horizon, "r", linewidth=0.5)

    print(f"poni1: {poni1}")
    print(f"poni2: {poni2}")
    return fig, ax


def save_nudge(poni_name: Path, ai: pyFAI.azimuthalIntegrator, nudge1: float = 0, nudge2: float = 0, orientation: int = 2):
    ai.poni1 = ai.poni1 + nudge1 * ai.pixel1
    ai.poni2 = ai.poni2 + nudge2 * ai.pixel2
    ai.detector = pyFAI.detectors.Detector(pixel1=ai.pixel1, pixel2=ai.pixel2,
                                           max_shape=ai.detector.shape, orientation=orientation)
    ai.save(poni_name)


def new(dist: float, poni1: float, poni2: float, shape: tuple, pixel1: float = 75e-6, pixel2: float = 75e-6,
        wavelength: float = 1.54185, rot1: float = 0, rot2: float = 0, rot3: float = 0, orientation: int = 2):
    detector = pyFAI.detectors.Detector(pixel1=pixel1, pixel2=pixel2, max_shape=shape, orientation=orientation)
    ai = pyFAI.azimuthalIntegrator.AzimuthalIntegrator(dist, poni1, poni2, rot1, rot2, rot3, pixel1, pixel2,
                                                       detector=detector, wavelength=wavelength, orientation=orientation)
    return ai


def save_new(poni_name: Path, dist: float, poni1: float, poni2: float,
             shape: tuple, pixel1: float = 75e-6, pixel2: float = 75e-6,
             wavelength: float = 1.54185, rot1: float = 0, rot2: float = 0, rot3: float = 0, orientation: int = 2):
    ai = new(dist, poni1, poni2, shape, pixel1, pixel2, wavelength, rot1, rot2, rot3, orientation)
    ai.save(poni_name)
