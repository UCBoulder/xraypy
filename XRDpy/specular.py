import matplotlib.offsetbox
import numpy as np
import fabio
import yaml
import pyFAI.detectors as detectors
from XRDpy.tiff_loader import Detector
from scipy.optimize import curve_fit, root_scalar
from scipy.special import erf
from pathlib import Path
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from matplotlib.animation import FuncAnimation
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

class SpecularScan:

    MAX_INT = (1 << 32) - 1

    def __init__(self, data_directory: Path, det_dist: float=150, anglular_range: float=1.5, beam_width: float=1,
                 standard_deviations: float=4, plot_name: str = "", plot_dpi: int = None, plot_title: str = None):
        self.plot_name = plot_name
        if plot_dpi is None:
            self.save_plots = False
        else:
            self.save_plots = True
        self.plot_title = plot_title
        self.plot_dpi = plot_dpi
        self.type = None
        self.z0 = 0
        self.det_dist=det_dist
        self.beam_width = beam_width
        self.data_directory = data_directory
        self.instrument = None
        self.detector = None
        self.file_type = None
        self.beam_center = None
        self.bc_sigma = None
        self.bc_amp = None
        self.where_max_angle = None
        self.perr = None
        self.max_angle = None
        self.standard_deviations = standard_deviations
        self.z_specular = None
        self.counts_total = None
        self.success_total = False
        self.success_dbeam = False
        self.success_specz = False
        self.determine_instrument()

        self.file_list = list(data_directory.glob("*" + self.file_type))

        self.pixel_size = self.detector.get_pixel1() * 1e3  # mm/pixel
        self.pixels_above = int(det_dist * np.radians(anglular_range) / self.pixel_size)   # number of pixels above the beam to keep (for crop)
        self.pixels_below = 8   # number of pixels below the beam to keep (for crop)
        self.pixel_horiontal_offset = 0
        self.z = np.arange(self.detector.shape[0])[::-1] * self.pixel_size
        
        self.base_mask = np.logical_not(self.detector.calc_mask())  # 1 keep, 0 mask (opposite of pyFAI standard, so it can be simply multiplied)

        self.angles = None
        self.z_motor = None
        self.intensity_data = None
        self.intensity_specular = None
        self.run()

    def save(self, directory: Path):
        if self.type == "om":
            np.save(directory / "angles.npy", self.angles)
            np.save(directory / "z_pos.npy", self.z)
            np.save(directory / "specular_data.npy", self.intensity_specular)
        elif self.type == "z":
            np.save(directory / "z_motor.npy", self.z_motor)
            np.save(directory / "z_pos.npy", self.z)
            np.save(directory / "specular_data.npy", self.intensity_specular)
        else:
            raise AttributeError("Has not determined type of specular scan.")
    
    def save_fig(self, fig):
        if self.save_plots:
            fig.savefig(self.plot_name, dpi=self.plot_dpi, bbox_inches="tight")

    def run(self):
        direct_beam_file_index = self.find_direct_beam_file_index()

        if direct_beam_file_index is None:
            with open(self.data_directory.parent / "beam_center.yaml", 'r') as yaml_file:
                beam_center_data = yaml.safe_load(yaml_file)
            self.beam_center = (beam_center_data["beamcenter"]['y'], beam_center_data["beamcenter"]['x'])
        
            sigma_x = beam_center_data["sigma"]['x']
            a_x = beam_center_data["amplitude"]['x']
            sigma_y = beam_center_data["sigma"]['y']
            a_y = beam_center_data["amplitude"]['y']
            print("Loaded beam center: ({}, {})".format(*self.beam_center))
            self.bc_sigma = sigma_y
            self.bc_amp = a_y
            self.angles, self.intensity_data = self.load()
        else:
            self.beam_center = self.find_beam_center(direct_beam_file_index)
            print("Direct beam is at ({}, {})".format(*self.beam_center))
            del(self.file_list[direct_beam_file_index])
            motor, self.intensity_data = self.load()

        if self.type == "om":
            self.angles = motor
        elif self.type == "z":
            self.z_motor = motor
        else:
            raise AttributeError("Scan type not established")

        self.process_data()
        self.fit()
    
    def determine_instrument(self):
        num_edf_files = len(list(self.data_directory.glob("*.edf")))
        num_tif_files = len(list(self.data_directory.glob("*.tif")))

        if num_edf_files > num_tif_files:
            self.instrument = "east"
            self.detector = detectors.Eiger2_1M()
            self.file_type = ".edf"
        else:
            self.instrument = "main"
            self.detector = Detector()
            self.file_type = ".tif"

    def find_direct_beam_file_index(self):
        possible_identifiers = ("db", "direct_beam")
        for ii, file in enumerate(self.file_list):
            if self.file_type == ".tif":
                if "direct_beam" in file.name:
                    return ii
            elif fabio.open(file).header["Comment"].replace("Base - ", "").lower() in possible_identifiers:
                return ii
        print("Did not find a direct beam file")
        return None
    
    def find_beam_center(self, direct_beam_index_file: int) -> tuple:
        intensity_db = self.load_image(self.file_list[direct_beam_index_file])
        rows = self.detector.shape[0]
        columns = self.detector.shape[1]
        
        px_x = np.arange(columns)
        db_x = np.sum(intensity_db, axis=0)
        px_y = np.arange(rows)
        db_y = np.sum(intensity_db, axis=1)
        def gaussian(x, x0, sigma, amplitude):
            arg = (x - x0) / sigma
            return amplitude * np.exp(-0.5 * arg * arg)
        (x0, sigma_x, a_x), _, = curve_fit(gaussian, px_x, db_x,
                                           p0=(0.5 * columns, 100, db_x.max()),
                                           nan_policy="omit")
        (y0, sigma_y, a_y), _ = curve_fit(gaussian, px_y, db_y,
                                           p0=(0.75 * rows, 100, db_y.max()),
                                           nan_policy="omit")
        beam_center = (y0, x0)
        self.bc_sigma = abs(sigma_y * self.detector.get_pixel1() * 1e3)
        print(self.bc_sigma)
        self.bc_amp = a_y
        to_write = {"beamcenter": {'x': float(x0), 'y': float(y0)},
                    "sigma": {'x': float(sigma_x), 'y': float(sigma_y)},
                    "amplitude": {'x': float(a_x), 'y': float(a_y)}}
        with open(self.data_directory.parent / "beam_center.yaml", 'w') as yaml_file:
            yaml.dump(to_write, yaml_file, default_flow_style=False)

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 7))
        fig.delaxes(ax4)
        fig.suptitle("Direct Beam", fontsize=12)

        ax3.set_title("Horizontal beam profile")
        ax3.scatter(px_x, db_x,
                    s=10, marker='o', edgecolors='k', lw=.75, facecolor='w')
        ax3.plot(px_x , gaussian(px_x, x0, sigma_x, a_x), "r")
        ax3.set_xlim(beam_center[1] - 30, beam_center[1] + 30)
        ax3.set_ylabel("Counts")
        
        ax2.set_title("Vertical beam profile")
        ax2.scatter(db_y, px_y,
                    s=10, marker='o', edgecolors='k', lw=.75, facecolor='w')
        ax2.plot(gaussian(px_y, y0, sigma_y, a_y), px_y, "r")
        ax2.set_ylim(beam_center[0] - 30, beam_center[0] + 30)
        ax2.set_xlabel("Counts")
        
        for ax in (ax2, ax3):
            ax.tick_params(axis='both', which='both', direction='in', right=True, top=True)
            # ax.set_xlabel("$\\omega$ motor position $(\\degree)$", fontsize=12)
            ax.grid(linestyle='dotted')
            ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

        ax1.set_facecolor('k')
        pos = ax1.imshow(intensity_db, norm=LogNorm(1, intensity_db.max()))
        ax1.set_xlim(beam_center[1] - 30, beam_center[1] + 30)
        ax1.set_ylim(beam_center[0] - 30, beam_center[0] + 30)
        for ii in range(1, 6):
            sig_x = ii * sigma_x
            sig_y = ii * sigma_y
            x = np.linspace(x0 - sig_x, x0 + sig_x, 1000)
            y = sig_y * np.sqrt(1 - (x - x0) ** 2 / (sig_x * sig_x))
            if ii < 3:
                color = "r"
            else:
                color = "w"
            ax1.plot(x, y + y0, color=color, linewidth=0.5)
            ax1.plot(x, -y + y0, color=color, linewidth=0.5)
        
        ax1.axvline(x0, color='r', linewidth=0.5)
        ax1.axhline(y0, color='r', linewidth=0.5)
        fig.colorbar(pos, ax=ax1)
        fig.tight_layout()

        print("")
        print("CLOSE THE WINDOW TO CONTINUE")
        plt.show()
        response = input("XRDpy>>> Was the beam center found? (Y/n): ")
        if not response:
            response = "y"
        self.save_fig(fig)
        return beam_center
    
    def load_image(self, filename: Path) -> np.ndarray:
        file_data = fabio.open(filename).data
        mask = self.base_mask.copy()
        if file_data.dtype == np.uint32:
            mask[np.where(file_data == self.MAX_INT)] = 0
        file_data *= mask
        return file_data
    
    def load(self):
        motor = np.empty(len(self.file_list))
        intensity_data = np.empty((len(self.file_list), *self.detector.shape))
        file_types = [""] * len(self.file_list)
        for ii, file in enumerate(self.file_list):
            if "om" in file.name:
                file_types[ii] = "om"
            elif "z" in file.name:
                file_types[ii] = "z"
            motor[ii] = self.get_angle(file)
            intensity_data[ii] = self.load_image(file)
        
        self.type = max(set(file_types), key=file_types.count)
        
        sorting_args = np.argsort(motor)
        motor = motor[sorting_args]
        intensity_data = intensity_data[sorting_args]

        return motor, intensity_data
    
    def get_angle(self, filename: Path) -> float:
        if self.file_type == ".edf":
            angle = fabio.open(filename).header["Comment"].replace("Base - ", "")
        elif self.file_type == ".tif":
            angle = self.angle_from_filename(filename.name)
        return float(angle)

    @staticmethod
    def angle_from_filename(filename: str):
        left_of_decimal = filename.split("_")[-3]
        angle = float(left_of_decimal)
        right_of_decimal = filename.split("_")[-2].replace(".tif", "")
        if left_of_decimal[0] == "-":
            angle -= float(right_of_decimal) / 10. ** len(right_of_decimal)
        else:
            angle += float(right_of_decimal) / 10. ** len(right_of_decimal)
        angle = round(angle, 3)
        return angle
       
    def process_data(self):
        if self.beam_center is not None:
            self.z -= (self.detector.shape[0] - self.beam_center[0]) * self.pixel_size
            bc_z = round(self.beam_center[0])
            z_lo = bc_z - self.pixels_above
            z_hi = bc_z + self.pixels_below

            self.z = self.z[z_lo:z_hi]

            if self.beam_width is None:
                data = self.intensity_data[:, z_lo:z_hi, :]
            else:
                half_width = 0.5 * self.beam_width / self.pixel_size
                x_lo = round(self.beam_center[1] - half_width)
                x_hi = round(self.beam_center[1] + half_width)
                x_lo += self.pixel_horiontal_offset
                x_hi += self.pixel_horiontal_offset
                data = self.intensity_data[:, z_lo:z_hi, x_lo:x_hi]
        else:
            data = self.intensity_data
        self.intensity_specular = np.sum(data, axis=2)
        # print(self.intensity_specular.max())

    def fit(self, z0=None, max_angle=1, pixel_cut=None, standard_deviations=None):
        if self.type == "om":
            self.fit_om(z0, max_angle, pixel_cut, standard_deviations)
        elif self.type == "z":
            self.fit_z(standard_deviations)
        else:
            raise AttributeError("Scan type was not established.")
        
    def check_goodness_of_om_fit(self):
        # func = self.specular_om_fit
        # x = self.where_max_angle
        # y = self.z_valid
        # p = (self.omega0, self.det_dist_fit)
        # err = self.perr
        fitted_curve = self.specular_om_fit(self.where_max_angle, self.omega0, self.det_dist_fit)
        # fitted_err = self.specular_om_error(self.where_max_angle, self.omega0, self.det_dist_fit, *self.perr)
        difference = np.abs(self.z_valid - fitted_curve)


    def fit_om(self, z0=None, max_angle=1, pixel_cut=None, standard_deviations=None):
        self.max_angle = max_angle
        if z0 is not None:
            self.z0 = z0
        if standard_deviations is not None:
            self.standard_deviations = standard_deviations

        self.counts_total = np.sum(self.intensity_specular.T, axis=0)
        self.counts_specular = np.sum(
            self.intensity_specular.T[np.where(self.z > self.standard_deviations * self.bc_sigma)],
            axis=0
        )
        self.counts_dbeam = np.sum(
            self.intensity_specular.T[np.where(self.z < self.standard_deviations * self.bc_sigma)],
            axis=0
        )
        
        print(f"z\u2080 = {self.z0}")
        max_ind = np.argmax(self.intensity_specular, axis=0)
        where_max_angle = self.angles[max_ind]

        # Remove data below the beam center
        valid = np.where(np.logical_and(
            self.z > self.z0,
            self.z < self.det_dist * np.radians(max_angle) + self.z0
            ))
        self.z_valid = self.z[valid]
        self.where_max_angle = where_max_angle[valid]
        if pixel_cut is not None:
            self.z_valid = self.z_valid[:-pixel_cut]
            self.where_max_angle = self.where_max_angle[:-pixel_cut]
        
        (self.omega0, self.det_dist_fit), pcov = curve_fit(self.specular_om_fit, self.where_max_angle, self.z_valid, p0=[0, self.det_dist])
        self.perr = np.sqrt(np.diag(pcov))
        print("Fit results:")
        print(f"    \u03C9\u2080 = ({self.omega0} \u00B1 {self.perr[0]})\u00B0")
        print(f"    d\u209B = ({self.det_dist_fit} \u00B1 {self.perr[1]}) mm")

    def fit_z(self, standard_deviations=None):
        if standard_deviations is not None:
            self.standard_deviations = standard_deviations
        self.counts_total = np.sum(self.intensity_specular.T, axis=0)
        self.counts_specular = np.sum(
            self.intensity_specular.T[np.where(self.z > self.standard_deviations * self.bc_sigma)],
            axis=0
        )
        self.counts_dbeam = np.sum(
            self.intensity_specular.T[np.where(self.z < self.standard_deviations * self.bc_sigma)],
            axis=0
        )

        print("Fit results:")
        fit_names = ["max", "z\u2080", "\u03C3\u2080"]
        unit_names = ["counts", "mm", "mm"]
        try:
            self.total_fit, pcov_total = curve_fit(self.occlusion_fit_single, self.z_motor, self.counts_total,
                                                   # p0=(1e6, 0.5, self.bc_sigma))
                                                   p0=(self.counts_total.max(), 0.5 * (self.z_motor.min() + self.z_motor.max()), self.bc_sigma))
            self.perr_total = np.sqrt(np.diag(pcov_total))
            print("  Total counts:")
            for name, fit_res, err, unit in zip(fit_names, self.total_fit, self.perr_total, unit_names):
                print(f"    {name} = ({fit_res:.5f} \u00B1 {err:.5f}) {unit}")
            self.success_total = True
        except RuntimeError:
            print("Failed to fit total counts")
            self.success_total = False

        try:
            self.dbeam_fit, pcov_dbeam = curve_fit(self.occlusion_fit_single, self.z_motor, self.counts_dbeam,
                                                   # p0=(1e6, 0.5, self.bc_sigma)
                                                   p0=(self.counts_dbeam.max(), 0.5 * (self.z_motor.min() + self.z_motor.max()), self.bc_sigma))
            self.perr_dbeam = np.sqrt(np.diag(pcov_dbeam))
            print("  Primary beam counts:")
            for name, fit_res, err, unit in zip(fit_names, self.dbeam_fit, self.perr_dbeam, unit_names):
                print(f"    {name} = ({fit_res:.5f} \u00B1 {err:.5f}) {unit}")
            self.success_dbeam = True
        except RuntimeError:
            print("Failed to fit direct beam counts")
            self.success_dbeam = False

        try:
            self.specz_fit, pcov_specz = curve_fit(self.specular_z_fit, self.z_motor, self.counts_specular,
                                                   # p0=(1e5, .4, .7, 0.5 * self.bc_sigma, 0.5 * self.bc_sigma))
                                                   p0=(self.counts_specular.max(), self.z_motor.min(), self.z_motor.max(), 0.5 * self.bc_sigma, 0.5 * self.bc_sigma))
            self.perr_specz = np.sqrt(np.diag(pcov_specz))
            fit_names = ["max", "z\u2081", "z\u2082", "\u03C3\u2081", "\u03C3\u2082"]
            unit_names = ["counts", "counts", "mm", "mm", "mm", "mm"]
            print("  Specular counts:")
            for name, fit_res, err, unit in zip(fit_names, self.specz_fit, self.perr_specz, unit_names):
                print(f"    {name} = ({fit_res:.5f} \u00B1 {err:.5f}) {unit}")
            self.success_specz = True
        except RuntimeError:
            print("Failed to fit specular counts")
            self.success_specz = False

        self.plot()

    def specular_om_fit(self, omega, omega0, det_dist):
        return det_dist * np.tan(2. * np.radians(omega - omega0)) + self.z0
    
    @staticmethod
    def specular_om_error(omega, omega0, det_dist, omega_err, dist_err):
        omega_center = omega - omega0
        dist_deriv = np.tan(2. * omega_center)
        sec_2omega = 1. / np.cos(2. * omega_center)
        omega_deriv = 2. * det_dist * sec_2omega * sec_2omega
        dist_term = dist_err * dist_deriv
        omega_term = omega_err * omega_deriv
        return np.sqrt(dist_term * dist_term + omega_term * omega_term)

    
    def specular_z_fit(self, z_motor, max_counts, z_lo, z_hi, sigma_lo, sigma_hi):
        return 0.5 * max_counts * (erf((z_motor - z_lo) / sigma_lo) - erf((z_motor - z_hi) / sigma_hi))
    
    def zero_angle(self, omega, omega0, det_dist):
        return det_dist * np.tan(np.radians(omega - omega0)) + self.z0
    
    def yoneda(self, omega, omega0, det_dist, critical_angle):
        return det_dist * np.tan(np.radians(omega - omega0 + critical_angle)) + self.z0
    
    def refraction_pos1(self, omega, omega0, det_dist, critical_angle):
        alpha = np.radians(omega - omega0)
        alpha_sq = alpha * alpha
        critical_angle = np.radians(critical_angle)
        crit_sq = critical_angle * critical_angle
        refraction_angle_sq = (alpha_sq - crit_sq) / (1 - 0.5 * crit_sq * alpha_sq)
        return self.z0 + det_dist * np.tan(alpha - np.sqrt(refraction_angle_sq))
    
    def refraction_pos2(self, omega, omega0, det_dist, critical_angle):
        alpha = np.radians(omega - omega0)
        alpha_sq = alpha * alpha
        critical_angle = np.radians(critical_angle)
        crit_sq = critical_angle * critical_angle
        refraction_angle_sq = (alpha_sq - crit_sq) / (1 - 0.5 * crit_sq * alpha_sq)
        return self.z0 + det_dist * np.tan(alpha + np.sqrt(refraction_angle_sq))
    
    def refraction_neg(self, omega, omega0, det_dist, critical_angle):
        alpha = np.radians(omega - omega0)
        critical_angle = np.radians(critical_angle)
        crit_sq = critical_angle * critical_angle
        refraction_angle = np.sqrt(alpha * alpha * (1 - 0.5 * crit_sq) + crit_sq)
        return self.z0 + det_dist * np.tan(alpha + refraction_angle)
    
    def plot(self, title="", critical_angle=None, horizon=False, det_dist=None, omega0=None):
        if self.type == "om":
            self.plot_om(title, critical_angle, horizon, det_dist, omega0)
        elif self.type == "z":
            self.plot_z()
        else:
            raise AttributeError("Type of specular not determined")

    def plot_om(self, title="", critical_angle=None, horizon=False, det_dist=None, omega0=None):
        if det_dist is None:
            det_dist = self.det_dist_fit
        if omega0 is None:
            omega0 = self.omega0

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 7))
        if self.plot_title:
            fig.suptitle(self.plot_title, fontsize=12)

        """FIT THROUGH MAX COUNT IN EACH ROW"""
        ax1.scatter(self.where_max_angle, self.z_valid, s=10, marker='o',
                    edgecolors='k', lw=.75, facecolor='w')
        omega = np.linspace(self.where_max_angle[-1] - 0.02, self.where_max_angle[0] + 0.02, 100)
        ax1.plot(omega, self.specular_om_fit(omega, self.omega0, self.det_dist_fit), "r")
        
        ax1.set_ylabel("$z$ (mm)", fontsize=12)
        ax1.set_title("where max pixel occurs", fontsize=12)
        # ax.legend(title="pixel")
        annotation_text = f"$\\omega_0 = {self.omega0:.4f} \\pm {self.perr[0]:.4f}^\\circ$\n$d_{{sd}} = {self.det_dist_fit:.2f} \\pm {self.perr[1]:.2f}$ mm"
        ax1.text(omega.min(), 0.95 * self.max_angle, annotation_text, transform=ax1.transAxes, fontsize=12,
                 verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
        
        """TRADITIONAL OM SCAN"""
        ax2.scatter(self.angles, self.counts_total,
                    s=10, marker='o', edgecolors='k', lw=.75, facecolor='w')
        ax2.set_title("Total Counts")

        ax3.scatter(self.angles, self.counts_specular,
                    s=10, marker='o', edgecolors='k', lw=.75, facecolor='w')
        ax3.set_title("Counts above beam")

        ax4.scatter(self.angles, self.counts_dbeam,
                    s=10, marker='o', edgecolors='k', lw=.75, facecolor='w')
        ax4.set_title("Counts in Direct Beam")
        
        for ax in (ax1, ax2, ax3, ax4):
            ax.tick_params(axis='both', which='both', direction='in', right=True, top=True)
            ax.set_xlabel("$\\omega$ motor position $(\\degree)$", fontsize=12)
            ax.grid(linestyle='dotted')
            ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        for ax in (ax2, ax3, ax4):
            ax.set_ylabel("Counts")
        fig.tight_layout()
        
        self.save_fig(fig)

        """SPECULAR"""

        fig, ax = plt.subplots(1, 1, figsize=(10, 7))
        ax.set_facecolor('k')
        if self.plot_title:
            fig.suptitle(self.plot_title, fontsize=12)

        # z = np.arange(intensity_x_total.shape[1])[::-1] * self.pixel_size
        ax.set_ylabel("$z$ (mm)", fontsize=12)

        color_map = ax.pcolormesh(self.angles, self.z, self.intensity_specular.T,
                                  norm=LogNorm(1, self.intensity_specular.max()), cmap="plasma")
        
        ax.axhline(self.z0, linestyle="--", color="#FF5349", linewidth=0.7)
        ax.axhline(self.standard_deviations * self.bc_sigma,
                   color="w", linewidth=0.5, linestyle="--")

        omega2 = np.linspace(omega0, omega0 + .75, 1000)
        ax.plot(omega2, self.specular_om_fit(omega2, omega0, det_dist), "white", linewidth=1, alpha=0.5)
        if horizon:
            ax.plot(omega2, self.zero_angle(omega2, omega0, det_dist), "white", linewidth=1, alpha=0.5)
        if critical_angle is not None:
            if isinstance(critical_angle, float):
                critical_angle = [critical_angle]
            last = None
            for crit in critical_angle:
                if crit == last:
                    omega1 = np.linspace(omega0 + crit, omega0 + crit + .05, 100)
                    ax.plot(omega1, self.refraction_pos2(omega1, omega0, det_dist, crit), "white", linewidth=1, alpha=0.5)
                else:
                    ax.plot(omega2, self.yoneda(omega2, omega0, det_dist, crit), "white", linewidth=1, alpha=0.5)
                    omega1 = np.linspace(omega0 + crit, self.angles[-1], 1000)
                    ax.plot(omega1, self.refraction_pos1(omega1, omega0, det_dist, crit), "white", linewidth=1, alpha=0.5)
                    omega1 = np.linspace(self.angles[0], omega0, 1000)
                    ax.plot(omega1, self.refraction_neg(omega1, omega0, det_dist, crit), "white", linewidth=1, alpha=0.5)
                    
                    # spec_at = self.specular_fit(crit + self.omega0, self.omega0, self.det_dist_fit)
                    # ax.plot([crit + self.omega0, crit + self.omega0], [spec_at + 0.6, spec_at + 1], "white", linewidth=.5, alpha=0.5)
                last = crit
        color_bar = fig.colorbar(color_map, ax=ax)

        ax.set_xlim(self.angles.min(), self.angles.max())

        ax.set_xlabel("$\\omega\\ (\\degree)$", fontsize=12)
        ax.set_ylabel("$z$ (mm)", fontsize=12)
        ax.set_title(title)
        fig.tight_layout()
        self.save_fig(fig)
        return fig, ax
    
    def plot_z(self, figsize=(10, 7)):
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        ax1.set_facecolor('k')
        if self.plot_title:
            fig.suptitle(self.plot_title, fontsize=12)

        z_m_fit = np.linspace(self.z_motor.min(), self.z_motor.max(), 500)

        annotation_texts = []
        locs = []
        fitted_axes = []

        # z = np.arange(intensity_x_total.shape[1])[::-1] * self.pixel_size
        ax1.set_ylabel("$z$ (mm)", fontsize=12)
        ax1.set_ylabel("$z$-motor (mm)", fontsize=12)
        #ax1.set_xticks(list(self.z_motor))
        color_map = ax1.pcolormesh(self.z_motor, self.z, self.intensity_specular.T,
                                   norm=LogNorm(1, self.intensity_specular.max()), cmap="plasma")
        color_bar = fig.colorbar(color_map, ax=ax1)
        ax1.axhline(self.standard_deviations * self.bc_sigma,
                    color="w", linewidth=0.5)

        # fig2, ax2 = plt.subplots(1, 1, figsize=(4, 3))
        ax2.scatter(self.z_motor, self.counts_total,
                   s=50,  # marker size
                   marker="o",  # marker shape
                   edgecolors="black",  # marker edge color
                   lw=2,  # marker edge width
                   alpha=1,  # transparency
                   facecolor='w')
        if self.success_total:
            ax2.plot(z_m_fit, self.occlusion_fit_single(z_m_fit, *self.total_fit), "r")
            annotation_texts.append(f"$z_0 = ({self.total_fit[1]:.3f} \\pm {self.perr_total[1]:.3f})$ mm")
            locs.append("lower left")
            fitted_axes.append(ax2)
        ax2.set_title("Total", fontsize=12)

        ax3.scatter(self.z_motor, self.counts_specular,
                   s=50,  # marker size
                   marker="o",  # marker shape
                   edgecolors="black",  # marker edge color
                   lw=2,  # marker edge width
                   alpha=1,  # transparency
                   facecolor='w')
        if self.success_specz:
            ax3.plot(z_m_fit, self.specular_z_fit(z_m_fit, *self.specz_fit), "r")
            annotation_texts.append(f"$z_2 = ({self.specz_fit[1]:.3f} \\pm {self.perr_specz[1]:.3f})$ mm\n$z_1 = ({self.specz_fit[2]:.3f} \\pm {self.perr_specz[2]:.3f})$ mm\n$z_{{ave}}={0.5*(self.specz_fit[1]+self.specz_fit[2]):.3f}$ mm")
            locs.append("lower center")
            fitted_axes.append(ax3)
        ax3.set_title("Above Beam", fontsize=12)
        
        ax4.scatter(self.z_motor, self.counts_dbeam,
                   s=50,  # marker size
                   marker="o",  # marker shape
                   edgecolors="black",  # marker edge color
                   lw=2,  # marker edge width
                   alpha=1,  # transparency
                   facecolor='w')
        if self.success_dbeam:
            ax4.plot(z_m_fit, self.occlusion_fit_single(z_m_fit, *self.dbeam_fit), "r")
            annotation_texts.append(f"$z_0 = ({self.dbeam_fit[1]:.3f} \\pm {self.perr_dbeam[1]:.3f})$ mm")
            locs.append("lower left")
            fitted_axes.append(ax4)
        ax4.set_title("Direct Beam", fontsize=12)
        
        for ax, ann, loc in zip(fitted_axes, annotation_texts, locs):
            ax.tick_params(axis='both', which='both', direction='in', right=True, top=True)
            ax.set_xlabel("$z$-motor (mm)", fontsize=12)
            ax.set_ylabel("Counts", fontsize=12)
            ax.grid(linestyle='dotted')
            ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            anchored_text = matplotlib.offsetbox.AnchoredText(ann, loc=loc, prop=dict(size=12), frameon=True)
            anchored_text.patch.set_boxstyle("round,pad=0.5")
            anchored_text.patch.set_facecolor('white')
            anchored_text.patch.set_alpha(0.8)
            ax.add_artist(anchored_text)
        fig.suptitle("Specular z-scan")
        fig.tight_layout()
        self.save_fig(fig)
    
    def animate_tiffs(self, fps: float = 6.0, scale: float = 1):
        scale *= 0.2
        bc_z = round(self.beam_center[0])
        z_lo = bc_z - self.pixels_above
        z_hi = bc_z + self.pixels_below
        self.z = self.z[z_lo:z_hi]
        if self.beam_width is None:
            data = self.intensity_data[:, z_lo:z_hi, :]
        else:
            half_width = 0.5 * self.beam_width / self.pixel_size
            x_lo = round(self.beam_center[1] - half_width)
            x_hi = round(self.beam_center[1] + half_width)
            x_lo += self.pixel_horiontal_offset
            x_hi += self.pixel_horiontal_offset
            data = self.intensity_data[:, z_lo:z_hi, x_lo:x_hi]

        fig, ax = plt.subplots(1, 1, figsize=(scale * (x_hi - x_lo), scale * 0.5 * (z_hi - z_lo)))
        fig.tight_layout()
        ax.set_facecolor("k")

        frame_delay = 1000. / fps

        im = ax.imshow(data[0], norm=LogNorm(1, data.max()))
        cbar = fig.colorbar(im, ax=ax)

        def _init():
            im.set_data(data[0])
            return im,

        def _update(ii):
            im.set_array(data[ii])
            return im,

        ani = FuncAnimation(fig, _update, frames=data.shape[0],
                            init_func=_init, interval=frame_delay, blit=True)

        return fig, ani

    def occlusion_fit_single(self, z_motor, max_counts, z0, sigma):
        return 0.5 * (max_counts - max_counts * erf((z_motor - z0) / sigma))

    def occlusion_fit_double(self, z_motor, max_counts, z_lo, z_hi, sigma_lo, sigma_hi):
        return 0.5 * max_counts * (erf((z_motor - z_hi) / sigma_hi) - erf((z_motor - z_lo) / sigma_lo)) + max_counts



    
    
    
