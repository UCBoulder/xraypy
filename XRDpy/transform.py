import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from pathlib import Path
import fabio
from pyFAI.geometry import Geometry
import gixpy as gp
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


def open_tiff(file_name: Path) -> np.ndarray:
    return fabio.open(file_name).data


def show_tiff(tiff: np.ndarray, title=None):
    plt.figure()
    plt.imshow(np.log10(tiff + 1))
    if title is not None:
        plt.title(title)

def save_data(directory: Path, filename: str, numpy_array: np.ndarray):
    if filename[-4:] != ".csv":
        filename.replace('.', '_')
        filename += ".csv"
    numpy_array.tofile(directory / filename, sep=',')


class TransformGIX:

    def __init__(self, incident_angle_degrees: float, tilt_angle_degrees: float=0.0, waveguiding=False):
        self.calibration_poni = Geometry()
        self.incident_angle = np.radians(incident_angle_degrees)
        self.yaw_angle = 0.0
        if waveguiding:
            self.pitch_angle = self.incident_angle
        else:
            self.pitch_angle = 0.0
        self.tilt_angle = np.radians(tilt_angle_degrees)

    def load(self, poni_file_name: Path):
        self.calibration_poni.load(poni_file_name)
        self.shape = self.calibration_poni.get_shape()

    def transform_image(self, data, flat_field):
        data_t, _ =  gp.transform(
            data.astype(np.float64),
            self.incident_angle,
            self.calibration_poni.get_pixel1(),
            self.calibration_poni.get_pixel2(),
            self.calibration_poni.get_poni1(),
            self.calibration_poni.get_poni2(),
            self.calibration_poni.get_dist(),
            self.yaw_angle,
            self.pitch_angle,
            self.tilt_angle, 2
        )
        flat_t, beam_center =  gp.transform(
            flat_field.astype(np.float64),
            self.incident_angle,
            self.calibration_poni.get_pixel1(),
            self.calibration_poni.get_pixel2(),
            self.calibration_poni.get_poni1(),
            self.calibration_poni.get_poni2(),
            self.calibration_poni.get_dist(),
            self.yaw_angle,
            self.pitch_angle,
            self.tilt_angle, 2
        )
            
        print(f"Tranformed beam center: ({beam_center[0]}, {beam_center[1]})")
        return (data_t, flat_t), beam_center
    
    def prep_transform(self, rot2=0.0):
        det_dist = self.calibration_poni.get_dist()
        # x is to the right from the PONI
        x = (np.arange(self.shape[1], dtype=np.float64) + 0.5) * self.calibration_poni.pixel2 - self.calibration_poni.poni2
        # z is up from the PONI
        z = ((self.shape[0] - 0.5) - np.arange(self.shape[0], dtype=np.float64).reshape(-1, 1)) * self.calibration_poni.pixel1 - self.calibration_poni.poni1
        r = np.sqrt(x * x + z * z)
        self.solid_angle = np.cos(np.arctan2(r, det_dist)) ** 3

        if self.tilt_angle:
            cos_tilt = np.cos(self.tilt_angle)
            sin_tilt = np.sin(self.tilt_angle)
            x_rot = x * cos_tilt - z * sin_tilt
            z = z * cos_tilt + x * sin_tilt
            x = x_rot

        phi = np.arctan2(x, det_dist)
        cos_phi = np.cos(phi)
        alpha = np.arctan2(z, det_dist) / cos_phi - self.incident_angle - rot2
        cos_alpha = np.cos(alpha)
        cos_incid = np.cos(self.incident_angle)

        q_xy_sq = cos_alpha * cos_alpha + (cos_incid - 2. * cos_alpha * cos_phi) * cos_incid # 2pi / lambda
        q_xy = np.sign(x) * np.sqrt(q_xy_sq)
        q_z = np.sin(alpha) + np.sin(self.incident_angle) # 2pi / lambda
        q = np.sqrt(q_xy_sq + q_z * q_z)
        r = det_dist * np.tan(2. * np.arcsin(0.5 * q))
        x = r * q_xy / (q * self.calibration_poni.pixel2)
        z = r * q_z / (q * self.calibration_poni.pixel1)
        
        self.col = np.floor(self.x)
        self.row = np.floor(self.z)
        remainder_col = self.x - self.col
        remainder_row = self.z - self.row
        remainder_col_comp = 1 - self.remainder_col
        remainder_row_comp = 1 - self.remainder_row
        self.weight_current = remainder_col_comp * remainder_row_comp
        self.weight_col_neighbor = remainder_col * remainder_row_comp
        self.weight_row_neighbor = remainder_col_comp * remainder_row
        self.weight_dia_neighbor = remainder_col * remainder_row
        
    def transform_image_python(self, data, flat_field, rot2=0.0):

        for rr in range(self.shape[0]):
            for cc in range(self.shape[1]):
                

        