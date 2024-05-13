import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from pathlib import Path
import fabio
from pyFAI.geometry import Geometry
try:
    import gixpy as gp
except ImportError:
    pass
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

    def __init__(self, incident_angle_degrees, tilt_angle_degrees=None):
        if tilt_angle_degrees is None:
            tilt_angle_degrees = 0.
        self.calibration_poni = Geometry()
        self.incident_angle_d = incident_angle_degrees
        self.tilt_angle_d = tilt_angle_degrees
        self.alpha_scattered = np.empty(0, dtype=np.float64)  # will be array size of data
        self.phi_scattered = np.empty(0, dtype=np.float64)  # will be array size of data
        self.q_unit_vec = np.empty(0, dtype=np.float64)  # will be 3 x size of data
        self.shape = (0, 0)  # will be shape of date
        self.beam_y_px = 0.0  # will be poni-y in pixels
        self.beam_x_px = 0.0  # will be poni-x in pixels
        self.det_dist_px = 0.0  # will be sample-detector distance in units of pixels

    def load(self, poni_file_name: Path):
        self.calibration_poni.load(poni_file_name)
        self.shape = self.calibration_poni.get_shape()

    def transform_image(self, data, flat_field):
        try:
            data_t =  gp.transform(data.astype(np.float64),
                                   self.incident_angle_d,
                                   self.calibration_poni.get_pixel1(),
                                   self.calibration_poni.get_poni1(),
                                   self.calibration_poni.get_poni2(),
                                   self.calibration_poni.get_dist(),
                                   self.tilt_angle_d)[0]
            flat_t, beam_center =  gp.transform(flat_field.astype(np.float64),
                                                self.incident_angle_d,
                                                self.calibration_poni.get_pixel1(),
                                                self.calibration_poni.get_poni1(),
                                                self.calibration_poni.get_poni2(),
                                                self.calibration_poni.get_dist(),
                                                self.tilt_angle_d)
        except NameError:
            print("gixpy not installed, it is recommended to install it to increase speed")
            data_t, flat_t, beam_center = self.transform_image_slow(data, weight)
            
        print(f"Tranformed beam center: ({beam_center[0]}, {beam_center[1]})")
        return (data_t, flat_t), beam_center
    
    def transform_image_slow(self, data, weight):
        to_pixel = 1. / self.calibration_poni.get_pixel1()
        beam_center = np.array((data.shape[0] - self.calibration_poni.get_poni1(), self.calibration_poni.get_poni2())) * to_pixel
        det_dist = self.calibration_poni.get_dist() * to_pixel
        det_dist_sq = det_dist * det_dist
        incident_angle = np.radians(self.incident_angle_d)
        cos_incident = np.cos(incident_angle)
        tilt_angle = np.radians(self.tilt_angle_d)
        x = beam_center[1] - np.arange(data.shape[1])
        y = (beam_center[0] - np.arange(data.shape[0])).reshape((data.shape[0], 1))
        if tilt_angle:
            costilt = np.cos(tilt_angle)
            sintilt = np.sin(tilt_angle)
            x_rot = x * costilt - y * sintilt
            y = x * sintilt + y * costilt
            x = x_rot
        hori_travel_sq = det_dist_sq + x * x
        sin_phi = x / np.sqrt(hori_travel_sq)
        cos_phi = np.sqrt(1 - sin_phi * sin_phi)
        alpha_scattered= np.arcsin(y / np.sqrt(hori_travel_sq + y * y)) - incident_angle
        cos_alpha = np.cos(alpha_scattered)
        q_xq_sq = cos_alpha * cos_alpha + cos_incident * cos_incident - 2.0 * cos_alpha * cos_incident * cos_phi
        q_xy = np.sqrt(q_xq_sq) * np.sign(x)
        q_z = np.sin(alpha_scattered) + np.sin(incident_angle)
        q_z_sq = q_z * q_z
        q_sq = q_xq_sq + q_z_sq
        q_scaler = det_dist * np.sqrt(0.5 + 1.0 / (2.0 - q_sq))
        r_xy = q_xy * q_scaler
        r_z = q_z * q_scaler
        new_beam_center = np.array((np.max(r_z), np.max(r_xy)))
        min_r = np.array(np.min(r_z), np.min(r_xy))
        new_shape = (int(np.ceil(new_beam_center[0] - min_r[0])) + 1, int(np.ceil(new_beam_center[1] - min_r[1])) + 1)
        x_deto = beam_center[1] - r_xy
        y_deto = beam_center[0] - r_z
        x_floor = np.floor(x_deto)
        y_floor = np.floor(y_deto)
        x_remainder = x_deto - x_floor
        y_remainder = y_deto - y_floor
        x_remainder_compliment = 1 - x_remainder
        y_remainder_compliment = 1 - y_remainder
        current_pixel_weight = x_remainder_compliment * y_remainder_compliment
        x_neighbor_weight = x_remainder * y_remainder_compliment
        y_neighbor_weight = x_remainder_compliment * y_remainder
        diag_neighbor_weight = x_remainder * y_remainder
        transformed_data = np.zeros(new_shape, np.float64)
        transformed_weights = np.zeros(new_shape, np.float64)
        for rr in range(data.shape[0]):
            for cc in range(data.shape[1]):
                row = int(y_floor[rr, cc])
                col = int(x_floor[rr, cc])
                transformed_data[row, col] += data[rr, cc] * current_pixel_weight[rr, cc]
                transformed_weights[row, col] += weight[rr, cc] * current_pixel_weight[rr, cc]

                transformed_data[row, col + 1] += data[rr, cc] * x_neighbor_weight[rr, cc]
                transformed_weights[row, col + 1] += weight[rr, cc] * x_neighbor_weight[rr, cc]
                
                transformed_data[row + 1, col] += data[rr, cc] * y_neighbor_weight[rr, cc]
                transformed_weights[row + 1, col] += weight[rr, cc] * y_neighbor_weight[rr, cc]
                
                transformed_data[row + 1, col + 1] += data[rr, cc] * diag_neighbor_weight[rr, cc]
                transformed_weights[row + 1, col + 1] += weight[rr, cc] * diag_neighbor_weight[rr, cc]
        return transformed_data, transformed_weights, new_beam_center


    def transform_reciprocal(self, data, exposure_time, weights=None, scale=1, plot=True):
        image_shape = (int(self.shape[0] * scale), int(self.shape[1] * scale))
        q_xy = np.sum(self.q_unit_vec[:1] * self.q_unit_vec[:1], axis=1) * (2 * np.pi / (self.calibration_poni.get_wavelength() * 1e10)) * ((self.q_unit_vec[1] > 0) * 2 - 1)
        q_z = self.q_unit_vec[2] * (2 * np.pi / (self.calibration_poni.get_wavelength() * 1e10))
        # beam center from top left corner of detector in q
        extent = (-np.max(q_xy), -np.min(q_xy), np.min(q_z), np.max(q_z))
        # qxy_extent = np.linspace(extent[0], extent[1], image_shape[1])
        # qz_extent = np.linspace(extent[2], extent[3], image_shape[0])
        beam_center_qxy_tranfromed = np.max(q_xy)
        beam_center_qz_transformed = np.max(q_z)
        # move into detector from (top left corner is origin)
        qxy_det_origin = beam_center_qxy_tranfromed - q_xy
        qz_det_origin = beam_center_qz_transformed - q_z
        # rescale so largest q reaches edge of detector in pixels
        qxy_det_origin *= (image_shape[1] - 2) / np.max(qxy_det_origin)
        qz_det_origin *= (image_shape[0] - 2) / np.max(qz_det_origin)

        # floor the locations to get pixel locations
        qxy_det_origin_floor = np.floor(qxy_det_origin)
        qz_det_origin_floor = np.floor(qz_det_origin)
        x_det_px_loc = qxy_det_origin_floor.astype(int)
        y_det_px_loc = qz_det_origin_floor.astype(int)
        x_remainder = qxy_det_origin - qxy_det_origin_floor
        y_remainder = qz_det_origin - qz_det_origin_floor
        x_remainder_compliment = 1 - x_remainder
        y_remainder_compliment = 1 - y_remainder
        current_pixel_weight = x_remainder_compliment * y_remainder_compliment
        x_neighbor_weight = x_remainder * y_remainder_compliment
        y_neighbor_weight = x_remainder_compliment * y_remainder
        diag_neighbor_weight = x_remainder * y_remainder
        
        data = data.astype(np.float64)
        transformed_data = np.zeros(image_shape, np.float64)
        transformed_weights = np.zeros(image_shape, np.float64)
        if weights is None:
            weights = np.ones(data.shape, dtype=np.float64) * exposure_time
        else:
            weights = weights.astype(np.float64)
        for rr in range(self.shape[0]):
            for cc in range(self.shape[1]):
                row = y_det_px_loc[rr, cc]
                col = x_det_px_loc[rr, cc]
                transformed_data[row, col] += data[rr, cc] * current_pixel_weight[rr, cc]
                transformed_weights[row, col] += weights[rr, cc] * current_pixel_weight[rr, cc]

                transformed_data[row, col + 1] += data[rr, cc] * x_neighbor_weight[rr, cc]
                transformed_weights[row, col + 1] += weights[rr, cc] * x_neighbor_weight[rr, cc]
                
                transformed_data[row + 1, col] += data[rr, cc] * y_neighbor_weight[rr, cc]
                transformed_weights[row + 1, col] += weights[rr, cc] * y_neighbor_weight[rr, cc]
                
                transformed_data[row + 1, col + 1] += data[rr, cc] * diag_neighbor_weight[rr, cc]
                transformed_weights[row + 1, col + 1] += weights[rr, cc] * diag_neighbor_weight[rr, cc]
        if plot:
            adjuster = exposure_time / transformed_weights
            adjuster[np.where(adjuster == np.infty)] = 0
            trasnformed_data_adj = transformed_data * adjuster
            fig = plt.figure(figsize=(10, 5), facecolor="w")
            ax1 = plt.subplot()
            for ax in fig.get_axes():
                ax.tick_params(which='both', color='k', direction = 'in')
                ax.set_facecolor("b")
            ax1.set_xlabel(r"q$_\mathregular{xy}\ (\mathregular{\AA}^{-1})$")
            ax1.set_ylabel(r"q$_\mathregular{z}\ (\mathregular{\AA}^{-1})$")
            ax1.yaxis.set_ticks_position('both')
            ax1.xaxis.set_ticks_position('both')
            pos = ax1.imshow(trasnformed_data_adj+1, extent=extent, norm=LogNorm(1, np.max(trasnformed_data_adj)))
            fig.colorbar(pos, ax=ax1, shrink=0.7)
        return transformed_data, transformed_weights, extent
        
    

if __name__ == "__main__":
    from tiff_loader import load_from
    # directory = Path("C:\\Users\\Teddy\\OneDrive - UCB-O365\\Rogerslab3\\Teddy\\TPP Films\\BTB-TPP\\2024 Film Growth\\Film 1\\GIWAXS TT5-06")
    # transformer = TransformGIWAXS(0.23, .15)
    
    directory = Path("C:\\Users\\Teddy\\OneDrive - UCB-O365\\Rogerslab3\\Teddy\\TPP Films\\BTB-TPP\\2024 Film Growth\\Film 2\\XRD\\TT5-09\\Thick sio2\\rotation 1")
    # directory = Path("C:\\Users\\Teddy\\OneDrive - UCB-O365\\Rogerslab3\\Teddy\\TPP Films\\BTB-TPP\\2024 Film Growth\\Film 2\\XRD\\calibration")
    # directory = Path("C:\\Users\\Teddy\\OneDrive - UCB-O365\\Rogerslab3\\Teddy\\TPP Films\\BTB-TPP\\2024 Film Growth\\Film 2\\XRD\\TT5-09\\Thick sio2\\rotation 2\\non-grazing")
    transformer = TransformGIX(8.81)

    transformer.load(directory / "cal.poni")

    # show_tiff(transformer.q_unit[0])
    # show_tiff(transformer.q_unit[1])
    # plt.plot(transformer.q_unit[1, :, int(transformer.shape[1] / 2)])
    # plt.plot(transformer.q_unit[1, int(transformer.shape[0] / 2), :])
    
    # data = open_tiff(directory / "GIWAXS-30min-BTBaTPP-teddy-20240130.tif")
    
    data, weight = load_from(directory)

    data_t, weights_t = transformer.transform_image(data, np.max(data), weight, .9)
    # transformer.transform_cake(data, 10, weight)
    plt.show()
    
    import pyFAI
    pyFAI.detectors.Detector()