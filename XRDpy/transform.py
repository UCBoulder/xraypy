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

    def __init__(self, incident_angle_degrees, tilt_angle_degrees=None):
        if tilt_angle_degrees is None:
            tilt_angle_degrees = 0.
        self.calibration_poni = Geometry()
        self.incident_angle_d = incident_angle_degrees
        self.tilt_angle_d = tilt_angle_degrees
        self.alpha_scattered = np.empty(0, dtype=np.float64)  # will be array size of data
        self.phi_scattered = np.empty(0, dtype=np.float64)  # will be array size of data
        self.q_unit_vec = np.empty(0, dtype=np.float64)  # will be 3 x size of data
        # self.q_xy = np.empty(0, dtype=np.float64)   # will be array size of data
        self.shape = (0, 0)  # will be shape of date
        self.beam_y_px = 0.0  # will be poni-y in pixels
        self.beam_x_px = 0.0  # will be poni-x in pixels
        self.det_dist_px = 0.0  # will be sample-detector distance in units of pixels

    def load(self, poni_file_name: Path):
        self.calibration_poni.load(poni_file_name)
        self.shape = self.calibration_poni.get_shape()
        # self.beam_x_px = self.calibration_poni.get_poni2() / self.calibration_poni.get_pixel2()
        # self.beam_y_px = float(self.shape[0]) - self.calibration_poni.get_poni1() / self.calibration_poni.get_pixel1()
        # self.det_dist_px = self.calibration_poni.get_dist() / self.calibration_poni.get_pixel1()
        # self.calculate_scattering_angles()
        # self.calculate_q_vector()

    def transform_image(self, data, weight):
        data_t =  gp.transform(data.astype(np.float64),
                               self.incident_angle_d,
                               self.calibration_poni.get_pixel1(),
                               self.calibration_poni.get_poni1(),
                               self.calibration_poni.get_poni2(),
                               self.calibration_poni.get_dist(),
                               self.tilt_angle_d)[0]
        weights_t, beam_center =  gp.transform(weight.astype(np.float64),
                                               self.incident_angle_d,
                                               self.calibration_poni.get_pixel1(),
                                               self.calibration_poni.get_poni1(),
                                               self.calibration_poni.get_poni2(),
                                               self.calibration_poni.get_dist(),
                                               self.tilt_angle_d)
        print(f"Tranformed beam center: ({beam_center[0]}, {beam_center[1]})")
        return (data_t, weights_t), beam_center


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