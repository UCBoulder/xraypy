import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from pathlib import Path
import fabio
from pyFAI.geometry import Geometry
import copy
import pyFAI
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


class GIXS:
    def __init__(self, incident_angle_degrees: float, tilt_angle_degrees: float = 0.0, poni_file: Path = None):
        self.incident_angle = np.radians(incident_angle_degrees)
        self.tilt_angle = np.radians(tilt_angle_degrees)
        if poni_file is None:
            self.ai_original = None
            self.dir = None
        else:
            self.load(poni_file)
        self.ai = None
        self.data = None
        self.flat_field = None
        self.mask
        self.header = {}

    def load(self, poni_file: Path):
        self.ai_original = pyFAI.load(poni_file)
        try:
            self.dir = poni_file.parent
        except AttributeError:
            if "/" in poni_file:
                self.dir = Path("/".join(poni_file.split("/")[:-1]))
            elif "\\" in poni_file:
                self.dir = Path("\\".join(poni_file.split("\\")[:-1]))
            else:
                self.dir = Path.cwd()

    def transform(self, image: np.ndarray, flat_field: np.ndarray = None, waveguiding: bool = False, refraction_angle_degrees: float = None):
        if self.ai_original is None:
            raise AttributeError("Must load a poni file first using .load(poni_file: Path)")
        if refraction_angle_degrees is None:
            if waveguiding:
                refraction_angle = self.incident_angle
            else:
                refraction_angle = 0.0
        else:
            refraction_angle = np.radians(refraction_angle_degrees)
        if flat_field is None:
            flat_field = np.ones_like(image)
        self.prep_transform(refraction_angle)
        self.data, self.flat_field, new_poni = self.transform_python(image, flat_field, refraction_angle)
        self.ai = copy.deepcopy(self.ai_original)
        self.ai.poni1 = new_poni[0]
        self.ai.poni2 = new_poni[1]
        self.ai.detector = pyFAI.detectors.Detector(pixel1=self.ai.pixel1, pixel2=self.ai.pixel2, max_shape=self.data_t.shape, orientation=2)
        self.ai.save(self.dir / "GIXS.poni")
        self.header = {
            "IncidentAngle": self.incident_angle,
            "TiltAngle": self.tilt_angle,
            "WaveGuiding": waveguiding,
            "RefractionAngle": refraction_angle,
        }
        return self.data, self.flat_field

    def save_edf(self, filename: str, directory_override: Path = None):
        if self.data is None:
            raise AttributeError("Must transform an image first using .transform(image: np.ndarray)")
        if directory_override is None:
            directory = self.dir
        else:
            directory = directory_override
        filename = filename.rstrip(".edf")
        edf_data_obj = fabio.edfimage.EdfImage(data=self.data, header=self.header)
        edf_flat_obj = fabio.edfimage.EdfImage(data=self.flat_field, header=self.header)
        edf_data_obj.write(directory / (filename + "_transformed.edf"))
        edf_flat_obj.write(directory / (filename + "_flat_field_t.edf"))

    def load_mask(self, mask_file: Path = None):
        if self.data is None:
            raise AttributeError("Must transform an image first using .transform(image: np.ndarray)")
        if mask_file is None:
            self.mask = np.logical_not(self.flat_field)
        else:
            mask = fabio.open(mask_file).data.astype(bool)
            self.mask = np.logical_or(mask, np.logical_not(np.logical_not(self.flat_field)))

    def integrate2d(self, q_bins: int = 500, azimuthal_bins: int = 180, radial_range: tuple = None,
                    azimuth_range: tuple = None, unit: str = "q_A^-1"):
        if self.data is None:
            raise AttributeError("Must transform an image first using .transform(image: np.ndarray)")
        if self.mask is None:
            self.load_mask()
        cake = self.ai.integrate2d_ng(
            self.data, q_bins, azimuthal_bins,
            radial_range=None,   # In units specified below
            azimuth_range=None,  # Start from 180 degrees to start from the axis to the right
            mask=self.mask, flat=self.flat,
            error_model="poisson",  unit=unit,
            polarization_factor=None, correctSolidAngle=False,
        )
        return cake

    def integrate1d(self, file_to_save: str = "reduction", q_range: tuple = None, azimuth_range: tuple = None,
                    exposure_time: float = None, q_bins: int = 1000, unit="q_A^-1"):
        if self.data is None:
            raise AttributeError("Must transform an image first using .transform(image: np.ndarray)")
        if self.mask is None:
            self.load_mask()
        if azimuth_range is None:
            azimuth_range = (-180, 0)
        file_to_save += ".edf"
        if (self.dir / file_to_save).is_file():
            (self.dir / file_to_save).unlink()
        if exposure_time is None:
            normalization_factor = 1
        else:
            normalization_factor = exposure_time / 60
        redu = self.ai.integrate1d_ng(
            self.data, q_bins, 
            radial_range=q_range,   # In units specified below
            azimuth_range=azimuth_range,  # Start from 180 degrees to start from the axis to the right
            mask=self.mask, flat=self.flat, error_model="poisson",
            correctSolidAngle=False,
            unit=unit, filename=self.dir / file_to_save, normalization_factor=normalization_factor
        )
        return redu

    def sector(self, file_to_save: str = "sector", q_range: tuple = None,
               azimuth_range: tuple = None, center: float = None, size: float= None,
               exposure_time: float = None,
               q_bins: int = 1000, unit="q_A^-1"):
        if azimuth_range is None and center is None and size is None:
            raise ValueError("Must provide either azimuth_range or center and size")
        if azimuth_range is not None:
            azimuth_range = (-azimuth_range[1], -azimuth_range[0])
        else:
            azimuth_range = (-center - 0.5 * size, -center + 0.5 * size)
        file_to_save += "_({},{})".format(*azimuth_range)
        return self.integrate1d(file_to_save, q_range, azimuth_range, exposure_time, q_bins, unit)
    
    def transform_python(self, data: np.ndarray, flat_field: np.ndarray, refraction_angle: float = 0.0):
        det_dist = self.ai_original.get_dist()
        # x is to the right from the PONI
        x = (np.arange(self.shape[1], dtype=np.float64) + 0.5) * self.ai_original.pixel2 - self.ai_original.poni2
        # z is up from the PONI
        z = ((self.shape[0] - 0.5) - np.arange(self.shape[0], dtype=np.float64).reshape(-1, 1)) * self.ai_original.pixel1 - self.ai_original.poni1
        sec_2theta = np.sqrt(x * x + z * z + det_dist * det_dist) / det_dist
        self.solid_angle = sec_2theta * sec_2theta * sec_2theta

        if self.tilt_angle:
            cos_tilt = np.cos(self.tilt_angle)
            sin_tilt = np.sin(self.tilt_angle)
            x_rot = x * cos_tilt - z * sin_tilt
            z = z * cos_tilt + x * sin_tilt
            x = x_rot
        
        alpha = np.arctan2(z, det_dist) # alpha_s + alpha_i        
        phi = np.arctan2(x * np.cos(alpha), det_dist)
        cos_phi = np.cos(phi)
        alpha -= self.incident_angle
        cos_alpha = np.cos(alpha)
        internal_angle = self.incident_angle - refraction_angle
        cos_internal = np.cos(internal_angle)
        
        q_xy_sq = cos_alpha * cos_alpha + cos_internal * cos_internal - 2.0 * cos_internal * cos_alpha * cos_phi
        q_xy = np.sqrt(q_xy_sq) * np.sign(x)
        q_z = np.sin(alpha) + np.sin(internal_angle)
        q_sq = q_xy_sq + q_z * q_z
        
        q_scaler = det_dist * np.sqrt(4 - q_sq) / (2 - q_sq)
        conv_px_x = 1. / self.ai_original.pixel2
        conv_px_z = 1. / self.ai_original.pixel1
        r_xy = q_xy * q_scaler * conv_px_x
        r_z = q_z * q_scaler * conv_px_z

        x_deto = r_xy - r_xy.min()
        z_deto = r_z.max() - r_z
        x_floor = np.floor(x_deto)
        z_floor = np.floor(z_deto)
        x_r = x_deto - x_floor
        x_rc = 1. - x_r
        z_r = z_deto - z_floor
        z_rc = 1. - z_r
        self.row = z_floor.astype(int)
        self.col = x_floor.astype(int)
        self.weight_current = x_rc * z_rc
        self.weight_col_neighbor = x_r * z_rc
        self.weight_row_neighbor = x_rc * z_r
        self.weight_dia_neighbor = x_r * z_r

        self.shape_transform = (int(self.row.max() + 2), int(self.col.max() + 2))
        
        self.poni2_transform = (-r_xy.min() + 0.5) * self.ai_original.pixel2
        self.poni1_transform = (self.shape_transform[0] - 0.5 - r_z.max()) * self.ai_original.pixel1
        print("New PONI:")
        print("poni1: {}".format(self.poni1_transform))
        print("poni2: {}".format(self.poni2_transform))
        
        print(self.shape_transform)
        data_transform = np.zeros(self.shape_transform)
        flat_transform = np.zeros(self.shape_transform)
        data *= self.solid_angle
        for rr in range(data.shape[0]):
            for cc in range(data.shape[1]):
                row = self.row[rr, cc]
                col = self.col[rr, cc]
                data_transform[row, col] += data[rr, cc] * self.weight_current[rr, cc]
                flat_transform[row, col] += flat_field[rr, cc] * self.weight_current[rr, cc]

                data_transform[row + 1, col] += data[rr, cc] * self.weight_row_neighbor[rr, cc]
                flat_transform[row + 1, col] += flat_field[rr, cc] * self.weight_row_neighbor[rr, cc]

                data_transform[row, col + 1] += data[rr, cc] * self.weight_col_neighbor[rr, cc]
                flat_transform[row, col + 1] += flat_field[rr, cc] * self.weight_col_neighbor[rr, cc]

                data_transform[row + 1, col + 1] += data[rr, cc] * self.weight_dia_neighbor[rr, cc]
                flat_transform[row + 1, col + 1] += flat_field[rr, cc] * self.weight_dia_neighbor[rr, cc]
                
        return data_transform, flat_transform

