import matplotlib.offsetbox
import numpy as np
import fabio
import yaml
import pyFAI.detectors as detectors
from scipy.optimize import curve_fit, root_scalar
from scipy.special import erf
from pathlib import Path
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from matplotlib.animation import FuncAnimation
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

class SpecularOmega:

    MAX_INT = (1 << 32) - 1

    def __init__(self, data_directory, det_dist=150, anglular_range=1.5, beam_width=1):
        self.type = self.get_type()
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
        self.standard_deviations=4
        self.z_specular = None
        self.counts_total = None
        self.determine_instrument()

        self.file_list = list(data_directory.glob("*" + self.file_type))

        self.pixel_size = self.detector.get_pixel1() * 1e3  # mm/pixel
        self.pixels_above = int(det_dist * np.radians(anglular_range) / self.pixel_size)   # number of pixels above the beam to keep (for crop)
        self.pixels_below = 8   # number of pixels below the beam to keep (for crop)
        self.pixel_horiontal_offset = 0
        self.z = np.arange(self.detector.shape[0])[::-1] * self.pixel_size
        
        self.base_mask = np.logical_not(self.detector.calc_mask())

        self.angles = None
        self.z_motor = None
        self.intensity_data = None
        self.intensity_specular = None
        self.run()

    def get_type(self):
        return "om"

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

        
        self.z_motor = self.angles
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
            self.detector = detectors.Eiger1M()
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
        self.bc_sigma = sigma_y
        self.bc_amp = a_y
        to_write = {"beamcenter": {'x': float(x0), 'y': float(y0)},
                    "sigma": {'x': float(sigma_x), 'y': float(sigma_y)},
                    "amplitude": {'x': float(a_x), 'y': float(a_y)}}
        with open(self.data_directory.parent / "beam_center.yaml", 'w') as yaml_file:
            yaml.dump(to_write, yaml_file, default_flow_style=False)
        fig, ax = plt.subplots(1, 1)
        ax.set_facecolor('k')
        pos = ax.imshow(intensity_db, norm=LogNorm(1, intensity_db.max()))
        ax.set_xlim(beam_center[1] - 30, beam_center[1] + 30)
        ax.set_ylim(beam_center[0] - 30, beam_center[0] + 30)
        for ii in range(1, 6):
            sig_x = ii * sigma_x
            sig_y = ii * sigma_y
            x = np.linspace(x0 - sig_x, x0 + sig_x, 1000)
            y = sig_y * np.sqrt(1 - (x - x0) ** 2 / (sig_x * sig_x))
            if ii < 3:
                color = "r"
            else:
                color = "w"
            ax.plot(x, y + y0, color=color, linewidth=0.5)
            ax.plot(x, -y + y0, color=color, linewidth=0.5)
        
        ax.axvline(x0, color='r', linewidth=0.5)
        ax.axhline(y0, color='r', linewidth=0.5)
        fig.colorbar(pos, ax=ax)
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
        


        sorting_args = np.argsort(motor)
        motor = motor[sorting_args]
        intensity_data = intensity_data[sorting_args]
        self.angles = motor
        self.z_motor = motor

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

    def fit(self, z0=None, max_angle=1, pixel_cut=None):
        if z0 is not None:
            self.z0 = z0
        print(f"z\u2080 = {self.z0}")
        max_ind = np.argmax(self.intensity_specular, axis=0)
        where_max_angle = self.angles[max_ind]

        # Remove data below the beam center
        valid = np.where(np.logical_and(
            self.z > self.z0,
            self.z < self.det_dist * np.radians(max_angle) + self.z0
            ))
        z = self.z[valid]
        where_max_angle = where_max_angle[valid]
        if pixel_cut is not None:
            z = z[:-pixel_cut]
            where_max_angle = where_max_angle[:-pixel_cut]
        
        (self.omega0, self.det_dist_fit), pcov = curve_fit(self.specular_fit, where_max_angle, z, p0=[0, self.det_dist])
        perr = np.sqrt(np.diag(pcov))
        print("Fit results:")
        print(f"    \u03C9\u2080 = ({self.omega0} \u00B1 {perr[0]})\u00B0")
        print(f"    d\u209B = ({self.det_dist_fit} \u00B1 {perr[1]}) mm")

        fig, ax = plt.subplots(1, 1, figsize=(4, 3))

        ax.scatter(where_max_angle, z, s=10, marker='o',
                   edgecolors='k', lw=.75, facecolor='w')
        omega = np.linspace(where_max_angle[-1] - 0.02, where_max_angle[0] + 0.02, 100)
        ax.plot(omega, self.specular_fit(omega, self.omega0, self.det_dist_fit), "r")
        ax.tick_params(axis='both', which='both', direction='in', right=True, top=True)
        ax.set_ylabel("$z$ (mm)", fontsize=12)
        ax.set_xlabel("$\\omega$ motor position $(\\degree)$", fontsize=12)
        ax.grid(linestyle='dotted')
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ax.set_title("where max pixel occurs", fontsize=12)
        # ax.legend(title="pixel")
        annotation_text = f"$\\omega_0 = {self.omega0:.4f} \\pm {perr[0]:.4f}^\\circ$\n$d_{{sd}} = {self.det_dist_fit:.2f} \\pm {perr[1]:.2f}$ mm"
        ax.text(omega.min(), 0.95 * max_angle, annotation_text, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))

        fig.tight_layout()
        return fig, ax

    def specular_fit(self, omega, omega0, det_dist):
        return det_dist * np.tan(2. * np.radians(omega - omega0)) + self.z0
    
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
        if det_dist is None:
            det_dist = self.det_dist_fit
        if omega0 is None:
            omega0 = self.omega0
        fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
        ax.set_facecolor('k')

        # z = np.arange(intensity_x_total.shape[1])[::-1] * self.pixel_size
        ax.set_ylabel("$z$ (mm)", fontsize=12)

        color_map = ax.pcolormesh(self.angles, self.z, self.intensity_specular.T,
                                  norm=LogNorm(1, self.intensity_specular.max()), cmap="plasma")
        omega2 = np.linspace(omega0, omega0 + .75, 1000)
        ax.plot(omega2, self.specular_fit(omega2, omega0, det_dist), "white", linewidth=1, alpha=0.5)
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
        ax.set_ylabel("Counts in row", fontsize=12)
        ax.set_title(title)
        fig.tight_layout()
        return fig, ax
    
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


class SpecularZ(SpecularOmega):
    def __init__(self, data_directory, det_dist=150, standard_deviations=4, angular_range=1.5, beam_width=1.):
        self.fit_failed = False
        super().__init__(data_directory, det_dist, angular_range, beam_width)

    def get_type(self):
        return "z"

    def specular_fit(self, z_motor, max_counts, z_lo, z_hi, sigma_lo, sigma_hi):
        return 0.5 * max_counts * (erf((z_motor - z_lo) / sigma_lo) - erf((z_motor - z_hi) / sigma_hi))
    
    def occlusion_fit_single(self, z_motor, max_counts, z0, sigma):
        return 0.5 * (max_counts - max_counts * erf((z_motor - z0) / sigma))

    def occlusion_fit_double(self, z_motor, max_counts, z_lo, z_hi, sigma_lo, sigma_hi):
        return 0.5 * max_counts * (erf((z_motor - z_hi) / sigma_hi) - erf((z_motor - z_lo) / sigma_lo)) + max_counts

    def where_half(self, func, params):
        def where_is(z_0, params):
            return func(z_0, *params) - 0.5 * params[0]

        root = root_scalar(where_is, args=(params,), method="newton", x0=params[1])

        print(f"Found where counts are half-max at: {root.root}")

        return root.root


    def fit(self, standard_deviations=None):
        if standard_deviations is not None:
            self.standard_deviations = standard_deviations
        self.counts_total = np.sum(self.intensity_specular.T, axis=0)
        self.z_specular = np.sum(
            self.intensity_specular.T[np.where(self.z > abs(self.standard_deviations * self.bc_sigma * self.detector.get_pixel1() * 1e3))],
            axis=0
        )
        self.z_primary = np.sum(
            self.intensity_specular.T[np.where(self.z < abs(self.standard_deviations * self.bc_sigma * self.detector.get_pixel1() * 1e3))],
            axis=0
        )

        try:
            self.total_fit, pcov_total = curve_fit(self.occlusion_fit_single, self.z_motor, self.counts_total, p0=(1e6, .5, .01))
            self.dbeam_fit, pcov_dbeam = curve_fit(self.occlusion_fit_single, self.z_motor, self.z_primary, p0=(1e6, .5, .01))
            self.specz_fit, pcov_specz = curve_fit(self.specular_fit, self.z_motor, self.z_specular, p0=(1e6, .4, .7, .01, .01))

            self.perr_total = np.sqrt(np.diag(pcov_total))
            self.perr_dbeam = np.sqrt(np.diag(pcov_dbeam))
            self.perr_specz = np.sqrt(np.diag(pcov_specz))

            self.z_half_total = self.where_half(self.occlusion_fit_single, self.total_fit)
            self.z_half_dbeam = self.where_half(self.occlusion_fit_single, self.dbeam_fit)
    
            print("Fit results:")
            fit_names = ["max", "z\u2080", "\u03C3\u2080"]
            unit_names = ["counts", "mm", "mm"]
            print("  Total counts:")
            for name, fit_res, err, unit in zip(fit_names, self.total_fit, self.perr_total, unit_names):
                print(f"    {name} = ({fit_res:.5f} \u00B1 {err:.5f}) {unit}")
            print("  Primary beam counts:")
            for name, fit_res, err, unit in zip(fit_names, self.dbeam_fit, self.perr_dbeam, unit_names):
                print(f"    {name} = ({fit_res:.5f} \u00B1 {err:.5f}) {unit}")
            fit_names = ["max", "z\u2081", "z\u2082", "\u03C3\u2081", "\u03C3\u2082"]
            unit_names = ["counts", "counts", "mm", "mm", "mm", "mm"]
            print("  Specular counts:")
            for name, fit_res, err, unit in zip(fit_names, self.specz_fit, self.perr_specz, unit_names):
                print(f"    {name} = ({fit_res:.5f} \u00B1 {err:.5f}) {unit}")
        except RuntimeError:
            print("failed to fit")
            self.fit_failed = True

        self.plot()
    
    def plot(self, figsize=(10, 7)):
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        ax1.set_facecolor('k')

        z_m_fit = np.linspace(self.z_motor.min(), self.z_motor.max(), 500)

        # z = np.arange(intensity_x_total.shape[1])[::-1] * self.pixel_size
        ax1.set_ylabel("$z$ (mm)", fontsize=12)
        ax1.set_ylabel("$z$-motor (mm)", fontsize=12)
        #ax1.set_xticks(list(self.z_motor))
        color_map = ax1.pcolormesh(self.z_motor, self.z, self.intensity_specular.T,
                                   norm=LogNorm(1, self.intensity_specular.max()), cmap="plasma")
        color_bar = fig.colorbar(color_map, ax=ax1)
        ax1.axhline(abs(self.standard_deviations * self.bc_sigma * self.detector.get_pixel1() * 1e3),
                    color="w", linewidth=0.5)

        # fig2, ax2 = plt.subplots(1, 1, figsize=(4, 3))
        ax2.scatter(self.z_motor, self.counts_total,
                   s=50,  # marker size
                   marker="o",  # marker shape
                   edgecolors="black",  # marker edge color
                   lw=2,  # marker edge width
                   alpha=1,  # transparency
                   facecolor='w')
        if not self.fit_failed:
            ax2.plot(z_m_fit, self.occlusion_fit_single(z_m_fit, *self.total_fit), "r")
        ax2.set_title("Total", fontsize=12)

        ax3.scatter(self.z_motor, self.z_specular,
                   s=50,  # marker size
                   marker="o",  # marker shape
                   edgecolors="black",  # marker edge color
                   lw=2,  # marker edge width
                   alpha=1,  # transparency
                   facecolor='w')
        if not self.fit_failed:
            ax3.plot(z_m_fit, self.specular_fit(z_m_fit, *self.specz_fit), "r")
        ax3.set_title("Above Beam", fontsize=12)
        
        ax4.scatter(self.z_motor, self.z_primary,
                   s=50,  # marker size
                   marker="o",  # marker shape
                   edgecolors="black",  # marker edge color
                   lw=2,  # marker edge width
                   alpha=1,  # transparency
                   facecolor='w')
        if not self.fit_failed:
            ax4.plot(z_m_fit, self.occlusion_fit_single(z_m_fit, *self.dbeam_fit), "r")
        ax4.set_title("Direct Beam", fontsize=12)

        if not self.fit_failed:
            annotation_texts = (
                f"$z_0 = ({self.total_fit[1]:.3f} \\pm {self.perr_total[1]:.3f})$ mm\n$z_{{HM}} = ({self.z_half_total:.3f})$ mm",
                f"$z_2 = ({self.specz_fit[1]:.3f} \\pm {self.perr_specz[1]:.3f})$ mm\n$z_1 = ({self.specz_fit[2]:.3f} \\pm {self.perr_specz[2]:.3f})$ mm\n$z_{{ave}}={0.5*(self.specz_fit[1]+self.specz_fit[2]):.3f}$ mm",
                f"$z_0 = ({self.dbeam_fit[1]:.3f} \\pm {self.perr_dbeam[1]:.3f})$ mm\n$z_{{HM}} = ({self.z_half_dbeam:.3f})$ mm",
            )

            locs = ["lower left", "lower center", "lower left"]
            for ax, ann, loc in zip((ax2, ax3, ax4), annotation_texts, locs):
                ax.tick_params(axis='both', which='both', direction='in', right=True, top=True)
                ax.set_xlabel("$z$-motor (mm)", fontsize=12)
                ax.set_ylabel("Counts", fontsize=12)
                ax.grid(linestyle='dotted')
                ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                # ax.legend(title="pixel")

                anchored_text = matplotlib.offsetbox.AnchoredText(ann, loc=loc, prop=dict(size=12), frameon=True)
                anchored_text.patch.set_boxstyle("round,pad=0.5")
                anchored_text.patch.set_facecolor('white')
                anchored_text.patch.set_alpha(0.8)
                ax.add_artist(anchored_text)
                #ax.text(self.z_motor.min(), 0.5 * self.specz_fit[1], ann, transform=ax.transAxes, fontsize=12,
                #        verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
        fig.suptitle("Specular z-scan")
        fig.tight_layout()


if __name__ == "__main__":
    # om = Data(Path(r"C:\Users\Teddy\OneDrive - UCB-O365\Rogerslab3\Teddy\TPP Films\BTB-TPP\2024 Film Growth\Film 2\XRD\non-grazing on blank for comparison\tune\raw_tiffs"))
    # path = Path(r"C:\Users\Teddy\OneDrive - UCB-O365\Rogerslab3\Teddy\TPP Films\BTB-TPP\2024 Film Growth\Film 4 (east campus dsc)\GIWAXS\film1\alignment")
    # spec = SpecularOmega(path / "spec")
    # spec.plot(critical_angle=0.24)
    path = Path(r"C:\Users\Teddy\OneDrive - UCB-O365\Rogerslab3\Teddy\TPP Films\BAP-TPP\Film 20241014\GIWAXS\20241015\z")
    z = SpecularZ(path, 150.12, standard_deviations=5)
    fig, ani = z.animate_tiffs()
    plt.show()