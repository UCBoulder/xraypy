import numpy as np
from itertools import product
from pathlib import Path
from scipy.optimize import least_squares
import matplotlib.pylab as plt
import matplotlib
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams["mathtext.fontset"] = "cm"

print("To run fit on powder:")
print("hex = tpp_fitting.Hexagonal_WAXS(a, c, q_data, count_data, plot_title, weights)")
print("hex = tpp_fitting.Hexagonal_pXRD(a, c, q_data, count_data, plot_title, weights)")
print("")
print("to run for GI use HexagonalGI(a, c, q_data, count_data, plot_title, weights, det_dist_mm)")

def load_ras(filename: Path):
    try:
        return np.loadtxt(filename, comments="*")
    except UnicodeDecodeError:
        with open(filename, "rb")as file:
            data = []
            for line in file.readlines():
                if line.decode()[0] != "*":
                    data.append([float(x) for x in line.decode().strip('\r\n').split(" ")])
        return np.array(data)


LN2 = np.log(2.)
INV_ROOT_PI_LN2 = 1. / np.sqrt(LN2 * np.pi)
FWHM_SQ_TO_HALF_SIGMA_SQ = 1. / 8. * LN2
TWO_PI = 2. * np.pi


class Background:
    def __init__(self, centers, widths):
        self.centers = np.array(centers).reshape((len(centers), 1))
        self.widths = np.array(widths).reshape((len(widths), 1))


class NewBackground(Background):
    def __init__(self):
        # bkgd_peak_centers = [0.5747, 0.8418, 1.0727, 1.3745, 1.704, 2.0, 2.3]
        bkgd_peak_centers = [0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.3]
        bkgd_peak_widths = [.2, .2, .2, .2, .2, .2, .2]
        super().__init__(bkgd_peak_centers, bkgd_peak_widths)


class TPP:
    xlabel = r"$q\ (\mathregular{\AA}^{-1})$"
    ylabel = "Intensity"
    init_peak_height = None
    
    def __init__(self, a: float, c: float, q: np.ndarray, counts: np.ndarray, monoclinic=False,
                 weights:np.ndarray=None, det_dist:float=150.0, sample_size:float=0.,
                 wavelength:float=1.54185, name: str="", background:str="new"):
        self.giwaxs = False
        self.free_param_num = 0
        self.wavelength = wavelength
        self.fourpi_lambda = 4. * np.pi / wavelength
        self.q = q
        if isinstance(self.q, dict):
            self.theta = {}
            for sector_key in self.q.keys():
                self.theta[sector_key] = np.arcsin(self.q[sector_key] / self.fourpi_lambda)
        else:
            self.theta = np.arcsin(self.q / self.fourpi_lambda)  # half of scattering angle
        self.counts = counts
        self.name = name
        if weights is None:
            self.weights = np.sqrt(counts)
        else:
            self.weights = weights

        self.lorentz = None

        self.hex_strain_lock = False
        self.mono_strain_lock = False
        self.goni_lock = False

        """LATTICE PARAMETERS"""
        self.hex = np.array([a, c])
        if isinstance(monoclinic, np.ndarray):
            self.mono = monoclinic
        elif isinstance(monoclinic, list) or isinstance(monoclinic, tuple):
            self.mono = np.array(monoclinic)
        elif monoclinic:
            # a, b, c, beta
            self.mono = np.array([25.086, 25.913, 5.911, np.radians(95.97)])
        else:
            self.mono = None

        """PEAK HEIGHTS"""
        self.hex_peak_heights = None
        self.mono_peak_heights = None
        self.hex_multiplicity = None
        self.mono_multiplicity = None

        """MILLER"""
        self.hex_hkl = None
        self.mono_hkl = None
        self.hh_hk_kk_ll = None
        self.mono_hh = None
        self.mono_hk = None
        self.mono_kk = None
        self.mono_ll = None

        """OTHER FIT PARAMETERS"""
        self.hex_params = {
            "w0": 3e-3,         # inv-angstrom
            "grain": 1000,       # angstrom
            "voigt": 0.5,       # between 0 and 1
            "strain": 0.,       # unitless
            "goni": 0.,         # millimeter
            "size": sample_size,    # millimeter
            "dist": det_dist    # millimeter
        }
        if monoclinic:
            self.mono_params = {
                "w0": 3e-3,         # inv-angstrom
                "grain": 1000,       # angstrom
                "voigt": 0.5,       # between 0 and 1
                "strain": 0.,       # unitless
                # "goni": 0.,         # millimeter
                "size": sample_size,       # millimeter
                "dist": det_dist    # millimeter
            }
        else:
            self.mono_params = None

        """BACKGROUND(S)"""
        if background == "new":
            self.bkgd = NewBackground()
            if isinstance(q, dict):
                self.bkgd_heights = {}
                self.bkgd_const = {}
                self.bkgd_lin = 0
                for sector in self.q.keys():
                    self.bkgd_heights[sector] = np.zeros(self.bkgd.centers.shape)
                    self.bkgd_const[sector] = 1 
            else:
                self.bkgd_heights = np.zeros(self.bkgd.centers.shape)
                self.bkgd_const = 1
                self.bkgd_lin = 0
        
        self.show_data()
        print("In next cell:")
        if isinstance(q, dict):
            print(" - Set background constant with .set_background_constant(list_of_background_constants)")
            print(" - Set q_min and q_max to the desired range of q values for each sector to fit using `.set_range(q_min, q_max, sector_key)`.")
            print(" - INITIALIZE FIT with .initialize_fit()")
        else:
            print(" - Set background constant with .set_background_constant(background_constant)")
            print(" - Set q_min and q_max to the desired range of q values to fit using `.set_range(q_min, q_max)`.")

    def initialize_fit(self, hkl_override=None, q_buffer=1.0, azi_buffer=20.):
        self.lorentz = self.lorentz_polarization_factor()
        self.hex_multiplicity = []
        self.hex_hkl = []
        q_center = []
        if self.giwaxs:
            azi_list = []
        if hkl_override is None:
            for h, k, l in product(range(5), repeat=3):
                if h + k + l > 0 and h >= k:
                    q = 2 * np.pi * np.sqrt(4.0 / 3.0 * (h * h + h * k + k * k) / (self.hex[0] * self.hex[0]) + l * l / (self.hex[1] * self.hex[1]))
                    if self.q.min() <= q <= self.q.max() * q_buffer:
                        q_center.append(q)
                        self.hex_hkl.append([h, k, l])
                        self.hex_multiplicity.append(self.multiplicity_check(h, k, l))
        else:
            self.hex_hkl = hkl_override
            if self.giwaxs:
                for h, k, l in self.hex_hkl:
                    q_xy_sq = TWO_PI * TWO_PI * 4.0 / 3.0 * (h * h + h * k + k * k) / (self.hex[0] * self.hex[0])
                    q_z_sq = TWO_PI * TWO_PI * l * l / (self.hex[1] * self.hex[1])
                    q = np.sqrt(q_xy_sq + q_z_sq)
                    q_center.append(q)
                    # q_xy_list.append(q_xy)
                    # q_z_list.append(q_z)
                    self.hex_multiplicity.append(self.multiplicity_check(h, k, l))
                    azi = np.rad2deg(np.arctan2(np.sqrt(q_z_sq), np.sqrt(q_xy_sq)))
                    azi_list.append(azi)
                    # print(f"({h}{k}{l}) -- q = {q:.4f} -- psi = {azi:.2f}")
            else:
                for h, k, l in self.hex_hkl:
                    q = 2 * np.pi * np.sqrt(4.0 / 3.0 * (h * h + h * k + k * k) / (self.hex[0] * self.hex[0]) + l * l / (self.hex[1] * self.hex[1]))
                    q_center.append(q)
                    self.hex_multiplicity.append(self.multiplicity_check(h, k, l))

        self.hex_hkl = np.array(self.hex_hkl)
        q_center = np.array(q_center)
        self.hex_multiplicity = np.array(self.hex_multiplicity)
        self.hex_peak_heights = np.ones((len(self.hex_hkl), 1)) * self.init_peak_height
    
        sort_ind = np.argsort(q_center)
        q_center = q_center[sort_ind]
        self.hex_hkl = self.hex_hkl[sort_ind]
        self.hex_multiplicity = self.hex_multiplicity[sort_ind].reshape((len(self.hex_multiplicity), 1))

        print(self.hex_hkl.shape)
        hh_hk_kk = 4. / 3. * (self.hex_hkl[:, 0] * self.hex_hkl[:, 0] + self.hex_hkl[:, 0] * self.hex_hkl[:, 1] + self.hex_hkl[:, 1] * self.hex_hkl[:, 1])
        l_sq = (self.hex_hkl[:, 2] * self.hex_hkl[:, 2])
        self.hh_hk_kk_ll = np.column_stack((hh_hk_kk, l_sq))

        if self.giwaxs:
            azi_list = np.array(azi_list)[sort_ind]

            hkl_full = self.hex_hkl
            mult_full = self.hex_multiplicity
            hhhkkkll = self.hh_hk_kk_ll
            self.hex_hkl = {}
            self.hex_multiplicity = {}
            self.hh_hk_kk_ll = {}
            self.hex_peak_heights = {}
            for key in self.keys:
                self.hex_hkl[key] = []
                self.hex_multiplicity[key] = []      
                self.hh_hk_kk_ll[key] = []
                sector_start, sector_end = self.sectors[key]
                sector_start -= azi_buffer
                sector_end += azi_buffer
                print(f"Sector: '{key}' with range ({round(sector_start)}--{round(sector_end)})\u00B0:")
                for ii, (h, k, l) in enumerate(hkl_full):
                    if azi_list[ii] < sector_end and azi_list[ii] > sector_start:
                        print(f"({h}{k}{l}) -- q = {q_center[ii]:.4f} -- \u03C8 = {azi_list[ii]:2.2f}\u00B0")
                        self.hex_hkl[key].append(hkl_full[ii])
                        self.hex_multiplicity[key].append(mult_full[ii])
                        self.hh_hk_kk_ll[key].append(hhhkkkll[ii])
                self.hex_peak_heights[key] = np.ones_like(self.hex_multiplicity[key]) * self.init_peak_height
            return None

        print('Number of possible reflections within data: %d' % len(self.hex_hkl))
        for ii in range(len(self.hex_hkl)):
            print(f'hkl: ({self.hex_hkl[ii, 0]},{self.hex_hkl[ii, 1]},{self.hex_hkl[ii, 2]}), q: {q_center[ii]:.3f}, multiplicity: {self.hex_multiplicity[ii, 0]}')

        if self.mono is not None:
            self.mono_multiplicity = []
            self.mono_hkl = []
            q_center = []
            csc_beta = 1. / np.sin(self.mono[3])
            cos_beta = np.cos(self.mono[3])
            # # for h, k, l in ((2,0,0), (-1,2,0), (3,2,0), (4,0,0), (4,1,0), (0,1,1), (1,-1,1), (2,-4,0), (4,3,0), (1,-2,1), (2,-1,1), (4,1,1), (3,-4,1), (4,-3,1), (4,-1,1), (1,-4,2), (2,4,2)):
            # for h, k, l in ((1,-1,0), (2,0,0), (3,2,0), (4,0,0), (4,1,0), (0,1,1), (1,-1,1), (1,-2,1), (2,-1,1), (4,-3,1), (4,-1,1)):
            # # for h, k, l in product(range(-4, 5), repeat=3):
            # # for h, k, l in ((0,2,0), (2,0,0), (0,4,0), (4,0,0), (0,1,1), (2,-1,1), (2,-2,1)):
            # # for h, k, l in ((2,0,0), (1,2,0), (2,1,0), (4,0,0), (2,4,0), (0,4,1), (4,4,0), (3,4,1), (4,4,1)):
            #     q = 2 * np.pi * np.sqrt(
            #         (h * h / (self.mono[0] * self.mono[0])
            #          + 2 * h * k * cos_beta / (self.mono[0] * self.mono[1])
            #          + k * k / (self.mono[1] * self.mono[1])
            #          ) * csc_beta * csc_beta + l * l / (self.mono[2] * self.mono[2])
            #     )
            #     if self.q.min() <= q <= self.q.max():
            #         q_center.append(q)
            #         self.mono_hkl.append([h, k, l])
            #         if h + k == 0 or l == 0:
            #             self.mono_multiplicity.append([2])
            #         else:
            #             self.mono_multiplicity.append([4])
            for h, k, l in product(range(-4, 5)[::-1], repeat=3):
                q = 2 * np.pi * np.sqrt(
                    (h * h / (self.mono[0] * self.mono[0])
                     + 2 * h * k * cos_beta / (self.mono[0] * self.mono[1])
                     + k * k / (self.mono[1] * self.mono[1])
                     ) * csc_beta * csc_beta + l * l / (self.mono[2] * self.mono[2])
                )
                if (self.q.min() <= q <= self.q.max()) and q not in q_center:
                    q_center.append(q)
                    self.mono_hkl.append([h, k, l])
                    if h + k == 0 or l == 0:
                        self.mono_multiplicity.append([2])
                    else:
                        self.mono_multiplicity.append([4])
            self.mono_hkl = np.array(self.mono_hkl)
            q_center = np.array(q_center)
            self.mono_multiplicity = np.array(self.mono_multiplicity)
            self.mono_peak_heights = np.zeros((len(self.mono_hkl), 1)) * self.init_peak_height

            sort_ind = np.argsort(q_center)
            q_center = q_center[sort_ind]
            self.mono_hkl = self.mono_hkl[sort_ind]
            self.mono_multiplicity = self.mono_multiplicity[sort_ind].reshape((len(self.mono_multiplicity), 1))
            for ii, (h,k,l) in enumerate(self.mono_hkl):
                print(f"({h}{k}{l}): q = {q_center[ii]}")

            self.mono_hh = (self.mono_hkl[:, 0] * self.mono_hkl[:, 0]).reshape((len(self.mono_hkl), 1))
            self.mono_hk = (self.mono_hkl[:, 0] * self.mono_hkl[:, 1]).reshape((len(self.mono_hkl), 1))
            self.mono_kk = (self.mono_hkl[:, 1] * self.mono_hkl[:, 1]).reshape((len(self.mono_hkl), 1))
            self.mono_ll = (self.mono_hkl[:, 2] * self.mono_hkl[:, 2]).reshape((len(self.mono_hkl), 1))
            print('Number of possible monoclinic reflections within data: %d' % len(self.mono_hkl))
            for ii in range(len(self.mono_hkl)):
                print(f'hkl: ({self.mono_hkl[ii, 0]},{self.mono_hkl[ii, 1]},{self.mono_hkl[ii, 2]}), q: {q_center[ii]:.3f}, multiplicity: {self.mono_multiplicity[ii, 0]}')
        
    
    def calc_widths(self, params, theta1):
        """Return FWHM^2 wrt q-space"""
        # cos_theta = np.cos(theta1)
        # width_grain = 0.9 * self.wavelength / (params["grain"] * cos_theta)
        # width_strain = params["strain"] * np.tan(theta1)
        # width_size = params["size"] / params["dist"] * np.tan(2. * theta1)
        #fwhm = params["w0"] + width_grain + width_strain
        #fwhm_sq = fwhm * fwhm
        instrumental = params["w0"] # * np.cos(theta1)
        width_grain = TWO_PI * 0.9 / params["grain"]
        width_strain = TWO_PI * params["strain"] * np.sin(theta1) / self.wavelength
        fwhm_sq = instrumental * instrumental + width_grain * width_grain + width_strain * width_strain
        # convert_to_q = 0.5 * self.fourpi_lambda * cos_theta
        return fwhm_sq  # * convert_to_q * convert_to_q
    
    def current_fit(self, sector=None, mixture=False):
        if sector is None:
            if mixture:
                hexagonal = self.fitting_function(self.hex, self.hex_peak_heights, self.hex_params,
                                                  self.bkgd_heights, self.bkgd_const, self.bkgd_lin,
                                                  self.mono, np.zeros_like(self.mono_peak_heights), self.mono_params)
                monoclinic = self.fitting_function(self.hex, np.zeros_like(self.hex_peak_heights), self.hex_params,
                                                  self.bkgd_heights, self.bkgd_const, self.bkgd_lin,
                                                  self.mono, self.mono_peak_heights, self.mono_params)
                return hexagonal, monoclinic
            else:
                return self.fitting_function(self.hex, self.hex_peak_heights, self.hex_params,
                                             self.bkgd_heights, self.bkgd_const, self.bkgd_lin,
                                             self.mono, self.mono_peak_heights, self.mono_params)
        else:
            return self.fitting_function(self.hex, self.hex_peak_heights[sector], self.hex_params,
                                         self.bkgd_heights[sector], self.bkgd_const[sector], self.bkgd_lin, sector)

    def fitting_function(self, hex, hex_peak_heights, hex_params, bkgd_heights, bkgd_const, bkgd_lin, mono=None, mono_peak_heights=None, mono_params=None):
        inv_d = np.sqrt(np.sum(self.hh_hk_kk_ll / (hex * hex), axis=1, keepdims=True))       # inverse d from miller indices (column vector)
        theta1 = np.arcsin(0.5 * self.wavelength * inv_d)
        fwhm_sq = self.calc_widths(hex_params, theta1)
        sigma_sq = fwhm_sq * FWHM_SQ_TO_HALF_SIGMA_SQ
        width_gauss = 2 * sigma_sq
        width_lortz = 0.25 * fwhm_sq
        # hwhm = 0.5 * np.sqrt(fwhm_sq)
        
        if hex_params["goni"]:
            q_c = 4 * np.pi / self.wavelength * np.sin(theta1 + hex_params["goni"] * np.cos(theta1) / hex_params["dist"])
        else:
            q_c = 2 * np.pi * inv_d                                     # find q-centers from inv_d
        q_shift = self.q - q_c                                  # center q from the q-centers
        q_shift_sq = q_shift * q_shift

        hex_peaks = np.sum(
            # self.hex_multiplicity * hex_peak_heights * ((1 - hex_params["voigt"]) * np.exp(-0.5 * q_shift_sq / sigma_sq) / np.sqrt(2. * np.pi * sigma_sq) + hex_params["voigt"] * hwhm / (np.pi * (q_shift_sq / 0.25 * fwhm_sq + 1))),
            self.hex_multiplicity * hex_peak_heights * ((1 - hex_params["voigt"]) * np.exp(-q_shift_sq / width_gauss) + hex_params["voigt"] * INV_ROOT_PI_LN2 / (q_shift_sq / width_lortz + 1)),
            # self.hex_multiplicity * hex_peak_heights * ((1 - hex_params["voigt"]) * np.exp(-0.5 * q_shift_sq / fwhm_sq) + hex_params["voigt"] * fwhm_sq / (q_shift_sq + fwhm_sq)),
            # self.hex_multiplicity * hex_peak_heights * ((1 - hex_params["voigt"]) * np.exp(-0.5 * q_shift_sq / sigma_sq) + hex_params["voigt"] * fwhm_sq / (q_shift_sq + fwhm_sq)),
            axis=0
        )
        if mono is not None:
            csc_beta = 1. / np.sin(mono[3])
            inv_d = np.sqrt(
                (self.mono_hh / (mono[0] * mono[0]) 
                - 2. * self.mono_hk * np.cos(mono[3]) / (mono[0] * mono[1])
                + self.mono_kk / (mono[1] * mono[1])) * (csc_beta * csc_beta)
                + self.mono_ll / (mono[2] * mono[2])
            )
            theta1 = np.arcsin(0.5 * self.wavelength * inv_d)
            fwhm_sq = self.calc_widths(mono_params, theta1)
            # hwhm = 0.5 * np.sqrt(fwhm_sq)
            sigma_sq = fwhm_sq * FWHM_SQ_TO_HALF_SIGMA_SQ
            width_gauss = 2 * sigma_sq
            width_lortz = 0.25 * fwhm_sq

            if hex_params["goni"]:
                q_c = 4 * np.pi / self.wavelength * np.sin(theta1 + hex_params["goni"] * np.cos(theta1) / hex_params["dist"])
            else:
                q_c = 2 * np.pi * inv_d     # find q-centers from inv_d            
            q_shift = self.q - q_c          # center q from the q-centers
            q_shift_sq = q_shift * q_shift
            mono_peaks = np.sum(
                # self.mono_multiplicity * mono_peak_heights * ((1 - mono_params["voigt"]) * np.exp(-0.5 * q_shift_sq / sigma_sq) / np.sqrt(2. * np.pi * sigma_sq) + mono_params["voigt"] * hwhm / (np.pi * (q_shift_sq / 0.25 * fwhm_sq + 1))),
                self.mono_multiplicity * mono_peak_heights * ((1 - mono_params["voigt"]) * np.exp(-q_shift_sq / width_gauss) + mono_params["voigt"] * INV_ROOT_PI_LN2 / (q_shift_sq / width_lortz + 1)),
                # self.mono_multiplicity * mono_peak_heights * ((1 - mono_params["voigt"]) * np.exp(-0.5 * q_shift_sq / fwhm_sq) + mono_params["voigt"] * fwhm_sq / (q_shift_sq + fwhm_sq)),
                # self.mono_multiplicity * mono_peak_heights * ((1 - mono_params["voigt"]) * np.exp(-0.5 * q_shift_sq / sigma_sq) + mono_params["voigt"] * fwhm_sq / (q_shift_sq + fwhm_sq)),
                axis=0
            )
        else:
            mono_peaks = 0

        return self.lorentz * (hex_peaks + mono_peaks + self.calc_background(bkgd_heights)) + bkgd_const + bkgd_lin * self.q
    
    def lorentz_polarization_factor(self):
        print("ERROR: no LP factor")
        return None

    def calc_background(self, bkgd_heights, sector=None):
        """Calculate the background pre-lorentz factor"""
        if sector is None:
            q = self.q
        else:
            q = self.q[sector]
        arg = (q - self.bkgd.centers) / self.bkgd.widths
        return np.sum(bkgd_heights * np.exp(-0.5 * arg * arg), axis=0)

    def show_data(self, fig_size=None):
        if isinstance(self.q, dict):
            if len(self.q) > 1:
                num_sectors = len(self.q.keys())
                if fig_size is None:
                    fig_size = (10, 2 * num_sectors)
                fig, axes = plt.subplots(2, round(0.5 * num_sectors))
                for ii, sector in enumerate(self.q.keys()):
                    row = int(0.5 * ii)
                    col = ii % 2
                    axes[row, col].scatter(
                        self.q[sector], self.counts[sector],
                        s=5,  # marker size
                        marker="o",  # marker shape
                        edgecolors="black",  # marker edge color
                        lw=.75,  # marker edge width
                        alpha=1,  # transparency
                        facecolor='w'  # marker face color
                    )
                    axes[row, col].set_title(self.name + " " + sector, fontsize=12)
            else:
                if fig_size is None:
                    fig_size = (5, 4)
                fig, ax = plt.subplots(1, 1, figsize=fig_size)
                for ii, sector in enumerate(self.q.keys()):
                    ax.scatter(
                        self.q[sector], self.counts[sector],
                        s=5,  # marker size
                        marker="o",  # marker shape
                        edgecolors="black",  # marker edge color
                        lw=.75,  # marker edge width
                        alpha=1,  # transparency
                        facecolor='w'  # marker face color
                    )
                    ax.set_title(self.name + " " + sector, fontsize=12)
                axes = ((ax,),)
        else:
            if fig_size is None:
                fig_size = (5, 4)
            fig, ax = plt.subplots(1, 1, figsize=fig_size)
            ax.scatter(
                self.q, self.counts,
                s=5,  # marker size
                marker="o",  # marker shape
                edgecolors="black",  # marker edge color
                lw=.75,  # marker edge width
                alpha=1,  # transparency
                facecolor='w'  # marker face color
            )
            ax.set_title(self.name, fontsize=12)
            axes = ((ax,),)
        for ax_row in axes:
            for ax in ax_row:
                ax.grid(linestyle='dotted')
                ax.set_xlabel(self.xlabel, fontsize=12)
                ax.set_ylabel(self.ylabel, fontsize=12)
                ax.set_yscale("log")
                ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)
        fig.tight_layout()
        if len(axes) == 1:
            axes = ax
        return fig, axes
    
    def show_fit(self, fig_size=None, mixture=False):
        if self.giwaxs:
            fig, axes = self.show_data(fig_size)
            if len(self.q) > 1:
                for ii, sector in enumerate(self.keys):
                    fit = self.current_fit(sector)
                    row = int(0.5 * ii)
                    col = ii % 2
                    axes[row, col].plot(self.q[sector], fit, label="Fit", color="red", lw="0.5")
                return fig, axes
            else:
                for ii, sector in enumerate(self.keys):
                    fit = self.current_fit(sector)
                    axes.plot(self.q[sector], fit, label="Fit", color="red", lw="0.5")
                return fig, axes
        elif mixture:
            hexagonal, monoclinic = self.current_fit(mixture=mixture)
            fig, ax = self.show_data(fig_size=fig_size)
            ax.plot(self.q, hexagonal, label="Hexaganol", color="red", lw="0.5")
            ax.plot(self.q, monoclinic, label="Monoclinic", color="g", lw="0.5")
            return fig, ax
        else:
            fit = self.current_fit()
            fig, ax = self.show_data(fig_size=fig_size)
            ax.plot(self.q, fit, label="Fit", color="red", lw="0.5")
            return fig, ax
    
    def save(self, name, title=None, dpi=600, fig_size=None):
        fig, ax = self.show_fit(fig_size=fig_size)
        if title is not None:
            ax.set_title(title)
        fig.savefig(name, dpi=dpi, bbox_inches="tight")
    
    def set_range(self, q_min, q_max, sector=None, q_buffer=1.0):
        if sector == None:
            self.counts = self.counts[(self.q > q_min) & (self.q < q_max)]
            self.weights = self.weights[(self.q > q_min) & (self.q < q_max)]
            self.theta = self.theta[(self.q > q_min) & (self.q < q_max)]
            self.q = self.q[(self.q > q_min) & (self.q < q_max)]
            self.initialize_fit(hkl_override=None, q_buffer=q_buffer)
            return self.show_fit()
        else:
            self.counts[sector] = self.counts[sector][(self.q[sector] > q_min) & (self.q[sector] < q_max)]
            self.weights[sector] = self.weights[sector][(self.q[sector] > q_min) & (self.q[sector] < q_max)]
            self.theta[sector] = self.theta[sector][(self.q[sector] > q_min) & (self.q[sector] < q_max)]
            self.q[sector] = self.q[sector][(self.q[sector] > q_min) & (self.q[sector] < q_max)]
            # self.initialize_fit(q_buffer)
        
    
    @staticmethod
    def multiplicity_check(h, k, l):
        print("ERROR: Multiplicity function missing")
        return None
    
    def report_peaks(self):
        q = 2 * np.pi * np.sqrt(np.sum(self.hh_hk_kk_ll / (self.hex * self.hex), axis=1, keepdims=True))
        warnings = 0
        for ii, (h, k, l) in enumerate(self.hex_hkl):
            if self.hex_peak_heights[ii] < 1e-10:
                warnings += 1
                print(f"Warning: ({h},{k},{l}) at q = {q[ii, 0]:5e}: {self.hex_peak_heights[ii, 0]:.5e}")
        if warnings:
            print("Remove peaks considering forbidden conditions (h+k=3n & l is odd) near the same q of these.")
            print("Use .remove_peak(h, k, l) to remove a peak. Do this above .fit_peak_heights() and re-run Notebook.")
            print("")
        for ii, (h, k, l) in enumerate(self.hex_hkl):
            print(f'({h},{k},{l}) at q = {q[ii, 0]:5e} inv A: {self.hex_peak_heights[ii, 0]:.5e}')
        if self.mono_hkl is not None:
            print("\nMONOCLINIC\n")
            csc_beta = 1. / np.sin(self.mono[3])
            q = 2 * np.pi * np.sqrt(
                (self.mono_hh / (self.mono[0] * self.mono[0]) 
                - 2. * self.mono_hk * np.cos(self.mono[3]) / (self.mono[0] * self.mono[1])
                + self.mono_kk / (self.mono[1] * self.mono[1])) * (csc_beta * csc_beta)
                + self.mono_ll / (self.mono[2] * self.mono[2])
            )
            warnings = 0
            for ii, (h, k, l) in enumerate(self.mono_hkl):
                if self.mono_peak_heights[ii] < 1e-10:
                    warnings += 1
                    print(f"Warning: ({h},{k},{l}) at q = {q[ii, 0]:5e}: {self.mono_peak_heights[ii, 0]:.5e}")
            if warnings:
                print("Remove peaks considering forbidden conditions (h+k=3n & l is odd) near the same q of these.")
                print("Use .remove_peak(h, k, l) to remove a peak. Do this above .fit_peak_heights() and re-run Notebook.")
                print("")
            for ii, (h, k, l) in enumerate(self.mono_hkl):
                print(f'({h},{k},{l}) at q = {q[ii, 0]:5e} inv A: {self.mono_peak_heights[ii, 0]:.5e}')

    def report_fit(self, fit):
        dof = len(self.counts) - len(fit.x)
        reduced_chi_sq = 2 * fit.cost / dof
        if reduced_chi_sq > 1e5:
            print(f"Reduced chi-squared: {reduced_chi_sq:.3e}")
        else:
            print(f"Reduced chi-squared: {reduced_chi_sq:.3f}")
        print(f"message: {fit.message}")
        print(f"Number of function evaluations: {fit.nfev}")
        print(f"Number of Jacobian evaluations: {fit.njev}\n")
    
    def chi_sq(self):
        if self.giwaxs:
            chi_sq = 0
            for sector in self.q.keys():
                counts_hat = self.current_fit(sector)
                chi_sq += np.sum(((self.counts[sector] - counts_hat) / self.weights[sector]) ** 2)
        else:
            counts_hat = self.fitting_function(self.hex, self.hex_peak_heights, self.hex_params, self.bkgd_heights, self.bkgd_const, self.bkgd_lin,
                                               self.mono, self.mono_peak_heights, self.mono_params)
            chi_sq = np.sum(((self.counts - counts_hat) / self.weights) ** 2)
        return chi_sq

    def report_chi_sq(self):
        chi_sq = self.chi_sq()
        print(f"Chi-squared: {chi_sq:.3f}")
        dof = len(self.counts) - self.free_param_num
        print(f"Degrees of freedom: {dof}")
        print(f"Reduced chi-squared: {chi_sq / dof:.3f}")
        
    def remove_peak(self, h, k, l, mono=False, sector=None):
        if mono:
            ind = np.where((self.mono_hkl[:, 0] == h) & (self.mono_hkl[:, 1] == k) & (self.mono_hkl[:, 2] == l))[0]
            if len(ind) == 1:
                self.mono_hkl = np.delete(self.mono_hkl, ind, axis=0)
                self.mono_multiplicity = np.delete(self.mono_multiplicity, ind, axis=0)
                self.mono_peak_heights = np.delete(self.mono_peak_heights, ind, axis=0)
                self.mono_hh = np.delete(self.mono_hh, ind, axis=0)
                self.mono_hk = np.delete(self.mono_hk, ind, axis=0)
                self.mono_kk = np.delete(self.mono_kk, ind, axis=0)
                self.mono_ll = np.delete(self.mono_ll, ind, axis=0)
            else:
                print("({},{},{}) Hexagonal peak already is not present.".format(h, k, l))
        elif sector is None:
            ind = np.where((self.hex_hkl[:, 0] == h) & (self.hex_hkl[:, 1] == k) & (self.hex_hkl[:, 2] == l))[0]
            if len(ind) == 1:
                self.hex_hkl = np.delete(self.hex_hkl, ind, axis=0)
                self.hex_multiplicity = np.delete(self.hex_multiplicity, ind, axis=0)
                self.hex_peak_heights = np.delete(self.hex_peak_heights, ind, axis=0)
                self.hh_hk_kk_ll = np.delete(self.hh_hk_kk_ll, ind, axis=0)
            else:
                print("({},{},{}) Hexagonal peak already is not present.".format(h, k, l))
        else:
            ind = np.where((self.hex_hkl[sector][:, 0] == h) & (self.hex_hkl[sector][:, 1] == k) & (self.hex_hkl[sector][:, 2] == l))[0]
            if len(ind) == 1:
                self.hex_hkl[sector] = np.delete(self.hex_hkl[sector], ind, axis=0)
                self.hex_multiplicity[sector] = np.delete(self.hex_multiplicity[sector], ind, axis=0)
                self.hex_peak_heights[sector] = np.delete(self.hex_peak_heights[sector], ind, axis=0)
                self.hh_hk_kk_ll[sector] = np.delete(self.hh_hk_kk_ll[sector], ind, axis=0)
            else:
                print("({},{},{}) Hexagonal peak already is not present.".format(h, k, l))


    def add_peak(self, h, k, l, mono=False, sector=None):
        if mono:
            self.mono_hkl = np.vstack((self.mono_hkl, [h, k, l]))
            if h + k == 0 or l == 0:
                self.mono_multiplicity = np.vstack((self.mono_multiplicity, 2))
            else:
                self.mono_multiplicity = np.vstack((self.mono_multiplicity, 4))
            self.mono_hh = np.vstack((self.mono_hh, h * h))
            self.mono_hk = np.vstack((self.mono_hk, h * k))
            self.mono_kk = np.vstack((self.mono_kk, k * k))
            self.mono_ll = np.vstack((self.mono_ll, l * l))
            self.mono_peak_heights = np.vstack((self.mono_peak_heights, self.init_peak_height))
        else:
            if sector is None:
                self.hex_hkl = np.vstack((self.hex_hkl, [h, k, l]))
                self.hex_peak_heights = np.vstack((self.hex_peak_heights, self.init_peak_height))
                self.hex_multiplicity = np.vstack((self.hex_multiplicity, self.multiplicity_check(h, k, l)))
                hh_hk_kk = 4. / 3. * (h * h + h * k + k * k)
                l_sq =  l * l
                self.hh_hk_kk_ll = np.vstack((self.hh_hk_kk_ll, np.array([hh_hk_kk, l_sq])))
            else:
                self.hex_hkl[sector] = np.vstack((self.hex_hkl[sector], [h, k, l]))
                self.hex_peak_heights[sector] = np.vstack((self.hex_peak_heights[sector], self.init_peak_height))
                self.hex_multiplicity[sector] = np.vstack((self.hex_multiplicity[sector], self.multiplicity_check(h, k, l)))
                hh_hk_kk = 4. / 3. * (h * h + h * k + k * k)
                l_sq =  l * l
                self.hh_hk_kk_ll[sector] = np.vstack((self.hh_hk_kk_ll[sector], np.array([hh_hk_kk, l_sq])))
    
    def clear_forbidden_peaks(self):
        for h, k, l in self.hex_hkl:
            # if h-k = 3n and l is odd
            if not (h - k) % 3 and l & 1:
                if not (h == 1 and k == 1 and l == 1):
                    print("Forbidden: ({},{},{})".format(h, k, l))
                    self.remove_peak(h, k, l)
    
    # def clear_forbidden_mono_peaks(self):
    #     to_clear = (
    #         (1,0,0), (1,1,0),  # below 0.5
    #         (2,2,0), (0,3,0), (1,3,0), (3,0,0), (3,1,1), (2,3,0), (3,2,0), (0,4,0)  # below 1
    #     )
    #     for hkl in to_clear:
    #         self.remove_peak(*hkl, mono=True)

    def set_background_constant(self, background_constant, background_linear=0):
        if isinstance(background_constant, list) or isinstance(background_constant, tuple):
            for ii, sector_key in enumerate(self.keys):
                self.bkgd_const[sector_key] = background_constant[ii]
            self.bkgd_lin = background_linear
            print("Keep adjusting and re-running this cell until the fit looks good.\n")
            print("In next cell:")
            print(f" - Set background with `.set_background([len {len(self.bkgd_heights)} list of [len {len(self.bkgd_heights['full'])} lists]])`.\n")
        else:
            self.bkgd_const = background_constant
            self.bkgd_lin = background_linear
            print("Keep adjusting and re-running this cell until the fit looks good.\n")
            print("In next cell:")
            print(f" - Set background with `.set_background([len {len(self.bkgd_heights)} array])`.\n")
    
    def set_background(self, background_heights):
        if self.giwaxs:
            for ii, sector in enumerate(self.keys):
                self.bkgd_heights[sector] = np.array(background_heights[ii]).reshape((len(background_heights[ii]), 1))
        else:
            self.bkgd_heights = np.array(background_heights).reshape((len(background_heights), 1))
        print("Keep adjusting and re-running this cell until the fit looks good.")
        print(f"The peak centers are at: {self.bkgd.centers.flatten()}\n")
        print("In next cell:")
        if self.giwaxs:
            print(" - Fit peak heights with .fit_lattice_parameters().")
        else:
            print(" - Fit peak heights with .fit_peak_heights().")
        return self.show_fit()
    
    def fit_peak_heights(self):
        print("In next cell:")
        print(" - Fit lattice parameters and peak heights with `.fit_lattice_parameters()`.\n")
        free_parameters = list(self.hex_peak_heights.flatten())
        self.free_param_num = len(free_parameters)
        def residuals(params):
            peak_heights = np.array(params).reshape((len(params), 1))
            counts_hat = self.fitting_function(self.hex, peak_heights, self.hex_params,
                                               self.bkgd_heights, self.bkgd_const, self.bkgd_lin,
                                               mono=None, mono_peak_heights=None, mono_params=None)
            return (counts_hat - self.counts) / self.weights

        fit = least_squares(residuals, free_parameters, method="lm")
        
        self.report_fit(fit)
        std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
        for ii in range(len(self.hex_hkl)):
            print(f'({self.hex_hkl[ii, 0]},{self.hex_hkl[ii, 1]},{self.hex_hkl[ii, 2]}): {fit.x[ii]:.5e} \u00B1 {std[ii]:.5e}')
        self.hex_peak_heights = fit.x.reshape((len(fit.x), 1))
        return self.show_fit()
    
    def add_monoclinic(self):
        print("In next cell:")
        print(" - Fit lattice parameters and peak heights with `.fit_lattice_parameters_monoclinic()`.\n")
        self.mono_peak_heights = np.ones_like(self.mono_peak_heights) * self.init_peak_height
        free_parameters = list(self.mono_peak_heights.flatten())
        self.free_param_num = len(free_parameters)
        def residuals(params):
            peak_heights = np.array(params).reshape((len(params), 1))
            counts_hat = self.fitting_function(self.hex, self.hex_peak_heights, self.hex_params,
                                               self.bkgd_heights, self.bkgd_const, self.bkgd_lin,
                                               self.mono, peak_heights, self.mono_params)
            return (counts_hat - self.counts) / self.weights
        
        fit = least_squares(residuals, free_parameters, method="lm")

        self.report_fit(fit)
        std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
        csc_beta = 1. / np.sin(self.mono[3])
        q = 2 * np.pi * np.sqrt(
                (self.mono_hh / (self.mono[0] * self.mono[0]) 
                - 2. * self.mono_hk * np.cos(self.mono[3]) / (self.mono[0] * self.mono[1])
                + self.mono_kk / (self.mono[1] * self.mono[1])) * (csc_beta * csc_beta)
                + self.mono_ll / (self.mono[2] * self.mono[2])
            )
        for ii in range(len(self.mono_hkl)):
            print(f'({self.mono_hkl[ii, 0]},{self.mono_hkl[ii, 1]},{self.mono_hkl[ii, 2]}): q={q[ii]} invA  -- --  {fit.x[ii]:.5e} \u00B1 {std[ii]:.5e}')
        self.mono_peak_heights = fit.x.reshape((len(fit.x), 1))
        return self.show_fit()
    
    def fit_free_params(self, keys=None):
        free_parameters = list(self.hex)
        if keys is not None:
            for key in keys:
                free_parameters.append(self.hex_params[key])
        N = len(free_parameters)
        free_parameters.extend(list(self.hex_peak_heights.flatten()))
        self.free_param_num = len(free_parameters)
        def residuals(params):
            hex_lattice = np.array(params[:2])
            hex_params = self.hex_params.copy()
            for ii in range(N-2):
                hex_params[keys[ii]] = params[2 + ii]
            peak_heights = np.array(params[N:]).reshape((len(self.hex_peak_heights), 1))
            counts_hat = self.fitting_function(hex_lattice, peak_heights, hex_params,
                                               self.bkgd_heights, self.bkgd_const, self.bkgd_lin,
                                               mono=None, mono_peak_heights=None, mono_params=None)
            return (counts_hat - self.counts) / self.weights
        
        fit = least_squares(residuals, free_parameters, method="lm")

        try:
            std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
            print(f"a = ({fit.x[0]:.4f} \u00B1 {std[0]:.4f}) Angstroms")
            print(f"c = ({fit.x[1]:.4f} \u00B1 {std[1]:.4f}) Angstroms")
            if keys is not None:
                for ii, key in enumerate(keys):
                    print(f"{key} = ({fit.x[2 + ii]:.4f} \u00B1 {std[2 + ii]:.4f})")
                    if key == "strain":
                        if fit.x[2 + ii] < 0.5 * std[2 + ii]:
                            print("Rejecting strain")
                            self.hex_strain_lock = True
                            return None
        except:
            print("failed to calculate error")
            print(f"a = {fit.x[0]:.4f} Angstroms")
            print(f"c = {fit.x[1]:.4f} Angstroms")
            if keys is not None:
                for ii, key in enumerate(keys):
                    print(f"{key} = {fit.x[2 + ii]:.4f}")
                    if key == "strain":
                        if fit.x[2 + ii] < 1e-10:
                            print("Rejecting strain")
                            self.hex_strain_lock = True
                            return None
        self.report_fit(fit)
        self.hex = np.array(fit.x[:2])
        if keys is not None:
            for ii, key in enumerate(keys):
                self.hex_params[key] = fit.x[2 + ii]
        self.hex_peak_heights = fit.x[N:].reshape((len(self.hex_peak_heights), 1))
        self.report_peaks()
        return self.show_fit()
    
    def fit_free_params_monoclinic(self, keys=None):
        free_parameters = list(self.hex)            # 2
        free_parameters.extend(list(self.mono))     # 4
        hex_keys = ["grain", "voigt"]
        if not self.hex_strain_lock:
            hex_keys.append("strain")
        if not self.goni_lock:
            hex_keys.append("goni")
        for k in hex_keys:
            free_parameters.append(self.hex_params[k])
        N = len(free_parameters)
        if keys is not None:
            for key in keys:
                free_parameters.append(self.hex_params[key])
        M = len(free_parameters)
        free_parameters.extend(list(self.hex_peak_heights.flatten()))
        P = len(free_parameters)
        free_parameters.extend(list(self.mono_peak_heights.flatten()))
        self.free_param_num = len(free_parameters)
        def residuals(params):
            hex_lattice = np.array(params[:2])
            mono_lattice = np.array(params[2:6])
            hex_params = self.hex_params.copy()
            mono_params = self.mono_params.copy()
            for ii in range(N - 6):
                hex_params[hex_keys[ii]] = params[6 + ii]
            for ii in range(M - N):
                mono_params[keys[ii]] = params[N + ii]
            hex_peak_heights = np.array(params[M:P]).reshape((len(self.hex_peak_heights), 1))
            mono_peak_heights = np.array(params[P:]).reshape((len(self.mono_peak_heights), 1))
            counts_hat = self.fitting_function(hex_lattice, hex_peak_heights, hex_params,
                                               self.bkgd_heights, self.bkgd_const, self.bkgd_lin,
                                               mono_lattice, mono_peak_heights, mono_params)
            return (counts_hat - self.counts) / self.weights
        
        fit = least_squares(residuals, free_parameters, method="lm")

        try:
            std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
            print("Hexagonal lattice parameters:")
            print(f"a = ({fit.x[0]:.4f} \u00B1 {std[0]:.4f}) Angstroms")
            print(f"c = ({fit.x[1]:.4f} \u00B1 {std[1]:.4f}) Angstroms")
            print("Monoclinic lattice parameters")
            print(f"a = ({fit.x[2]:.4f} \u00B1 {std[2]:.4f}) Angstroms")
            print(f"b = ({fit.x[3]:.4f} \u00B1 {std[3]:.4f}) Angstroms")
            print(f"c = ({fit.x[4]:.4f} \u00B1 {std[4]:.4f}) Angstroms")
            print(f"\u03B2 = ({np.rad2deg(fit.x[5]):.4f} \u00B1 {np.rad2deg(std[5]):.4f})\u00B0")
            for ii, key in enumerate(hex_keys):
                print(f"Hexagonal {key} = ({fit.x[6 + ii]:.4f} \u00B1 {std[6 + ii]:.4f})")
            if keys is not None:
                for ii, key in enumerate(keys):
                    print(f"Monoclinic {key} = ({fit.x[N + ii]:.4f} \u00B1 {std[N + ii]:.4f})")
                    if key == "strain":
                        if fit.x[N + ii] < 0.5 * std[N + ii]:
                            print("Rejecting strain")
                            self.mono_strain_lock = True
                            return None
        except:
            print("failed to calculate error")
            print("Hexagonal lattice parameters:")
            print(f"a = {fit.x[0]:.4f} Angstroms")
            print(f"c = {fit.x[1]:.4f} Angstroms")
            print("Monoclinic lattice parameters")
            print(f"a = {fit.x[2]:.4f} Angstroms")
            print(f"b = {fit.x[3]:.4f} Angstroms")
            print(f"c = {fit.x[4]:.4f} Angstroms")
            print(f"\u03B2 = {np.rad2deg(fit.x[5]):.4f}\u00B0")
            
            for ii, key in enumerate(hex_keys):
                print(f"Hexagonal {key} = {fit.x[6 + ii]:.4f}")
            if keys is not None:
                for ii, key in enumerate(keys):
                    print(f"Monoclinic {key} = {fit.x[N + ii]:.4f}")
                    if key == "strain":
                        if fit.x[N + ii] < 1e-10:
                            print("Rejecting strain")
                            self.mono_strain_lock = True
                            return None
        self.report_fit(fit)
        self.hex = np.array(fit.x[:2])
        self.mono = np.array(fit.x[2:6])
        for ii, key in enumerate(hex_keys):
            self.hex_params[key] = fit.x[6 + ii]
        if keys is not None:
            for ii, key in enumerate(keys):
                self.mono_params[key] = fit.x[N + ii]
        self.hex_peak_heights = fit.x[M:P].reshape((len(self.hex_peak_heights), 1))
        self.mono_peak_heights = fit.x[P:].reshape((len(self.mono_peak_heights), 1))
        self.report_peaks()
        return self.show_fit()
    
    def fit_background(self, keys=None):
        free_parameters = list(self.hex)
        if keys is not None:
            for key in keys:
                free_parameters.append(self.hex_params[key])
        N = len(free_parameters)
        free_parameters.extend(list(self.hex_peak_heights.flatten()))
        PH = len(free_parameters)
        free_parameters.extend(list(self.bkgd_heights.flatten()))
        free_parameters.append(self.bkgd_const)
        self.free_param_num = len(free_parameters)
        def residuals(params):
            hex_lattice = np.array(params[:2])
            hex_params = self.hex_params.copy()
            keys = params[2:N]
            for ii, key in enumerate(keys):
                hex_params[key] = params[2 + ii]
            peak_heights = np.array(params[N:PH]).reshape((len(self.hex_peak_heights), 1))
            bkgd_heights = np.array(params[PH:-1]).reshape((len(self.bkgd_heights), 1))
            bkgd_const = params[-1]
            counts_hat = self.fitting_function(hex_lattice, peak_heights, hex_params,
                                               bkgd_heights, bkgd_const, self.bkgd_lin,
                                               mono=None, mono_peak_heights=None, mono_params=None)
            return (counts_hat - self.counts) / self.weights
        
        fit = least_squares(residuals, free_parameters, method="lm")

        try:
            std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
            print(f"a = ({fit.x[0]:.4f} \u00B1 {std[0]:.4f}) Angstroms")
            print(f"c = ({fit.x[1]:.4f} \u00B1 {std[1]:.4f}) Angstroms")
            if keys is not None:
                for ii, key in enumerate(keys):
                    print(f"{key} = ({fit.x[2 + ii]:.4f} \u00B1 {std[2 + ii]:.4f})")
                    if key == "strain":
                        if fit.x[2 + ii] < 0.5 * std[2 + ii]:
                            print("Rejecting strain")
                            self.hex_strain_lock = True
                            return None
        except:
            print("failed to calculate error")
            print(f"a = {fit.x[0]:.4f} Angstroms")
            print(f"c = {fit.x[1]:.4f} Angstroms")
            if keys is not None:
                for ii, key in enumerate(keys):
                    print(f"{key} = ({fit.x[2 + ii]:.4f}")
                    if key == "strain":
                        if fit.x[2 + ii] < 1e-10:
                            print("Rejecting strain")
                            self.hex_strain_lock = True
                            return None
        self.report_fit(fit)
        self.hex = np.array(fit.x[:2])
        if keys is not None:
            for ii, key in enumerate(keys):
                self.hex_params[key] = fit.x[2 + ii]
        self.hex_peak_heights = fit.x[N:PH].reshape((len(self.hex_peak_heights), 1))
        self.bkgd_heights = fit.x[PH:-1].reshape((len(self.bkgd_heights), 1))
        self.bkgd_const = fit.x[-1]
        self.report_peaks()
        return self.show_fit()

    def fit_background_monoclinic(self, keys=None):
        free_parameters = list(self.hex)            # 2
        free_parameters.extend(list(self.mono))     # 4
        hex_keys = [self.hex_params["grain"], self.hex_params["voigt"]]
        if not self.hex_strain_lock:
            hex_keys.append(self.hex_params["strain"])
        if not self.goni_lock:
            hex_keys.append(self.hex_params["goni"])
        free_parameters.extend(hex_keys)
        N = len(free_parameters)
        if keys is not None:
            for key in keys:
                free_parameters.append(self.hex_params[key])
        M = len(free_parameters)
        free_parameters.extend(list(self.hex_peak_heights.flatten()))
        P = len(free_parameters)
        free_parameters.extend(list(self.mono_peak_heights.flatten()))
        P2 = len(free_parameters)
        free_parameters.extend(list(self.bkgd_heights.flatten()))
        free_parameters.append(self.bkgd_const)
        self.free_param_num = len(free_parameters)
        def residuals(params):
            hex_lattice = np.array(params[:2])
            mono_lattice = np.array(params[2:6])
            hex_params = self.hex_params.copy()
            mono_params = self.mono_params.copy()
            for ii in range(N - 6):
                hex_params[hex_keys[ii]] = params[6 + ii]
            for ii in range(M - N):
                mono_params[keys[ii]] = params[N + ii]
            hex_peak_heights = np.array(params[M:P]).reshape((len(self.hex_peak_heights), 1))
            mono_peak_heights = np.array(params[P:P2]).reshape((len(self.mono_peak_heights), 1))
            bkgd_heights = np.array(params[P2:-1]).reshape((len(self.bkgd_heights), 1))
            bkgd_const = params[-1]
            counts_hat = self.fitting_function(hex_lattice, hex_peak_heights, hex_params,
                                               bkgd_heights, bkgd_const, self.bkgd_lin,
                                               mono_lattice, mono_peak_heights, mono_params)
            return (counts_hat - self.counts) / self.weights
        
        fit = least_squares(residuals, free_parameters, method="lm")

        try:
            std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
            print("Hexagonal lattice parameters:")
            print(f"a = ({fit.x[0]:.4f} \u00B1 {std[0]:.4f}) Angstroms")
            print(f"c = ({fit.x[1]:.4f} \u00B1 {std[1]:.4f}) Angstroms")
            print("Monoclinic lattice parameters")
            print(f"a = ({fit.x[2]:.4f} \u00B1 {std[2]:.4f}) Angstroms")
            print(f"b = ({fit.x[3]:.4f} \u00B1 {std[3]:.4f}) Angstroms")
            print(f"c = ({fit.x[4]:.4f} \u00B1 {std[4]:.4f}) Angstroms")
            print(f"\u03B2 = ({np.rad2deg(fit.x[5]):.4f} \u00B1 {np.rad2deg(std[5]):.4f})\u00B0")
            for ii, key in enumerate(hex_keys):
                print(f"Hexagonal {key} = ({fit.x[6 + ii]:.4f} \u00B1 {std[6 + ii]:.4f})")
            if keys is not None:
                for ii, key in enumerate(keys):
                    print(f"Monoclinic {key} = ({fit.x[N + ii]:.4f} \u00B1 {std[N + ii]:.4f})")
                    if key == "strain":
                        if fit.x[2 + ii] < 0.5 * std[2 + ii]:
                            print("Rejecting strain")
                            self.mono_strain_lock = True
                            return None
        except:
            print("failed to calculate error")
            print(f"a = {fit.x[0]:.4f} Angstroms")
            print(f"c = {fit.x[1]:.4f} Angstroms")
            if keys is not None:
                for ii, key in enumerate(keys):
                    print(f"{key} = ({fit.x[2 + ii]:.4f}")
                    if key == "strain":
                        if fit.x[2 + ii] < 1e-10:
                            print("Rejecting strain")
                            self.mono_strain_lock = True
                            return None
        self.report_fit(fit)
        self.hex = np.array(fit.x[:2])
        self.mono = np.array(fit.x[2:6])
        for ii, key in enumerate(hex_keys):
            self.hex_params[key] = fit.x[6 + ii]
        if keys is not None:
            for ii, key in enumerate(keys):
                self.mono_params[key] = fit.x[N + ii]
        self.hex_peak_heights = fit.x[M:P].reshape((len(self.hex_peak_heights), 1))
        self.mono_peak_heights = fit.x[P:P2].reshape((len(self.mono_peak_heights), 1))
        self.bkgd_heights = fit.x[P2:-1].reshape((len(self.bkgd_heights), 1))
        self.bkgd_const = fit.x[-1]
        self.report_peaks()
        return self.show_fit()
    
    def fit_lattice_parameters(self):
        print("In next cell:")
        print(" - Fit peak widths with .fit_peak_widths().\n")
        self.fit_free_params()
    
    def fit_peak_widths(self):
        print("In next cell:")
        print(" - Fit peak shapes with `.fit_strain()`.\n")
        return self.fit_free_params(["grain"])

    def fit_strain(self):
        self.hex_strain_lock = False
        print("In next cell:")
        if self.goni_lock:
            print(" - Fit voigt parameter with `.fit_voigt()`.\n")
        else:
            print(" - Fit goniometer offset with `.fit_goni()`\n")
        return self.fit_free_params(["grain", "strain"])
    
    def fit_voigt(self):
        print("In next cell:")
        if self.mono is None:
            print(" - finish fit with `.fit_hexagonal()`.\n")
        else:
            print(" - add monoclinic peaks with `.add_monoclinic()`.\n")
        keys = ["grain", "voigt"]
        if not self.hex_strain_lock:
            keys.append("strain")
        if not self.goni_lock:
            keys.append("goni")
        return self.fit_free_params(keys)
    
    def fit_hexagonal(self):
        print("In next cell:")
        print(" - Save with  `.save(file_name, title=None, dpi=600)` or `fig, ax = self.show_fit(fig_size)`.")
        keys = ["grain", "voigt"]
        if not self.hex_strain_lock:
            keys.append("strain")
        if not self.goni_lock:
            keys.append("goni")
        return self.fit_background(keys)


class Powder(TPP):

    xlabel = r"$q\ (\mathregular{\AA}^{-1})$"
    ylabel = "Intensity (counts / second)"
    init_peak_height = 1

    @staticmethod
    def multiplicity_check(h, k, l) -> int:
        if h + k == 0:
            multiplicity = 2
        elif l == 0 and (k == 0 or h == k):
            multiplicity = 6
        elif h == k or l == 0 or k == 0:
            multiplicity = 12
        else:
            multiplicity = 24
        return multiplicity

    def lorentz_polarization_factor(self):
        """calculate the lorentz polarization factor"""
        sintheta = np.sin(self.theta)
        cos2theta = np.cos(2 * self.theta)
        return (1 + cos2theta * cos2theta) / (sintheta * sintheta * np.cos(self.theta))
    
    def fit_goni(self):
        self.goni_lock = False
        print("In next cell:")
        print(" - Fit voigt parameter with `.fit_voigt()`.\n")
        keys = ["grain", "goni"]
        if not self.hex_strain_lock:
            keys.append("strain")
        return self.fit_free_params(keys)
    
    def fit_lattice_parameters_monoclinic(self):
        print("In next cell:")
        print(" - Fit monoclinic widths with `.fit_widths_monoclinic().")
        return self.fit_free_params_monoclinic()
    
    def fit_widths_monoclinic(self):
        print("In next cell:")
        print(" - Fit monoclinic strain with `.fit_strain_monoclinic().")
        return self.fit_free_params_monoclinic(["grain"])
    
    def fit_strain_monoclinic(self):
        print("In next cell:")
        print(" - Fit monoclinic peak shapes with `.fit_voigt_monoclinic().")
        return self.fit_free_params_monoclinic(["grain", "strain"])
    
    def fit_voigt_monoclinic(self):
        print("In next cell:")
        print(" - Fit background with `.fit_monoclinic().")
        keys = ["grain", "voigt"]
        if not self.mono_strain_lock:
            keys.append("strain")
        result =  self.fit_free_params_monoclinic(keys)
        if result is None:
            print("Redoing fit with strain locked")
            result =  self.fit_free_params_monoclinic(keys)
        return result
    
    def fit_monoclinic(self):
        print("In next cell:")
        print(" - Save with  `.save(file_name, title=None, dpi=600)` or `fig, ax = self.show_fit(fig_size)`.")
        keys = ["grain", "voigt"]
        if not self.mono_strain_lock:
            keys.append("strain")
        return self.fit_background_monoclinic(keys)
    

class Monoclinic(Powder):

    init_peak_height = 5e-3

    def __init__(self, q: np.ndarray, counts: np.ndarray, lattice_parameters=None,
                 weights:np.ndarray=None, det_dist:float=150.0, sample_size:float=0.,
                 wavelength:float=1.54185, name: str="", background:str="new"):
        if lattice_parameters is None:
            lattice_parameters = True
        super().__init__(12, 10, q, counts, monoclinic=lattice_parameters, weights=weights, det_dist=det_dist,
                         sample_size=sample_size, wavelength=wavelength, name=name, background=background)
        self.mono_params["goni"] = 0.

    def initialize_fit(self, hkl_override=None, q_buffer=1.0, azi_buffer=20.):
        super().initialize_fit(hkl_override=None, q_buffer=1.0, azi_buffer=20.)
        self.mono_peak_heights = np.ones_like(self.mono_peak_heights) * self.init_peak_height
        self.hex_peak_heights = np.zeros_like(self.hex_peak_heights)

    def fit_free_params(self, keys=None):
        free_parameters = list(self.mono)
        if keys is not None:
            for key in keys:
                free_parameters.append(self.mono_params[key])
        N = len(free_parameters)
        free_parameters.extend(list(self.mono_peak_heights.flatten()))
        self.free_param_num = len(free_parameters)
        def residuals(params):
            mono_lattice = np.array(params[:4])
            mono_params = self.mono_params.copy()
            for ii in range(N-4):
                mono_params[keys[ii]] = params[4 + ii]
            peak_heights = np.array(params[N:]).reshape((len(self.mono_peak_heights), 1))
            counts_hat = self.fitting_function(self.hex, self.hex_peak_heights, self.hex_params,
                                               self.bkgd_heights, self.bkgd_const, self.bkgd_lin,
                                               mono=mono_lattice, mono_peak_heights=peak_heights, mono_params=mono_params)
            return (counts_hat - self.counts) / self.weights
        
        fit = least_squares(residuals, free_parameters, method="lm")

        try:
            std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
            for ii, lp in enumerate(("a", "b", "c")):
                print(f"{lp} = ({fit.x[ii]:.4f} \u00B1 {std[ii]:.4f}) Angstroms")
            print(f"\u03B2 = ({np.rad2deg(fit.x[3]):.4f} \u00B1 {np.rad2deg(std[3]):.4f})\u00B0")
            if keys is not None:
                for ii, key in enumerate(keys):
                    print(f"{key} = ({fit.x[4 + ii]:.4f} \u00B1 {std[4 + ii]:.4f})")
                    if key == "strain":
                        if fit.x[4 + ii] < 0.5 * std[4 + ii]:
                            print("Rejecting strain")
                            self.hex_strain_lock = True
                            return None
        except:
            print("failed to calculate error")
            for ii, lp in ("a", "b", "c"):
                print(f"{lp} = {fit.x[ii]:.4f} Angstroms")
            print(f"\u03B2 = {np.rad2deg(fit.x[3]):.4f}\u00B0")
            if keys is not None:
                for ii, key in enumerate(keys):
                    print(f"{key} = {fit.x[4 + ii]:.4f}")
                    if key == "strain":
                        if fit.x[4 + ii] < 1e-10:
                            print("Rejecting strain")
                            self.hex_strain_lock = True
                            return None
        self.report_fit(fit)
        self.mono = np.array(fit.x[:4])
        if keys is not None:
            for ii, key in enumerate(keys):
                self.mono_params[key] = fit.x[4 + ii]
        self.mono_peak_heights = fit.x[N:].reshape((len(self.mono_peak_heights), 1))
        self.report_peaks()
        return self.show_fit()
    

class pXRD(Powder):
    init_peak_height = 1

    def __init__(self, a: float, c: float, q: np.ndarray, counts: np.ndarray, monoclinic=False,
                 weights:np.ndarray=None, det_dist:float=150.0, sample_size:float=0.,
                 wavelength:float=1.54185, name: str="", background:str="new"):
        super().__init__(a, c, q, counts, monoclinic, weights, det_dist,
                         sample_size, wavelength, name, background)
        self.hex_params["w0"] = 0.003
        self.mono_params["w0"] = 0.003


class WAXS(Powder):
    init_peak_height = 5e-3

    def __init__(self, a: float, c: float, q: np.ndarray, counts: np.ndarray, monoclinic=False,
                 weights:np.ndarray=None, det_dist:float=150.0, sample_size:float=0.,
                 wavelength:float=1.54185, name: str="", background:str="new"):
        super().__init__(a, c, q, counts, monoclinic, weights, det_dist,
                         sample_size, wavelength, name, background)
        self.hex_params["w0"] = 0.003
        self.mono_params["w0"] = 0.003
        self.goni_lock = True


class Film(TPP):

    xlabel = r"$q_z\ (\mathregular{\AA}^{-1})$"
    ylabel = "Intensity (counts / second)"
    init_peak_height = 1

    @staticmethod
    def multiplicity_check(h, k, l) -> int:
        if h + k == 0:
            multiplicity = 2
        elif h == k or h == 0:
            multiplicity = 6
        else:
            multiplicity = 12
        return multiplicity

    def lorentz_polarization_factor(self, sector=None):
        """calculate the lorentz polarization factor"""
        if sector is None:
            twotheta = 2. * self.theta
        else:
            twotheta = 2. * self.theta[sector]
        cos2theta = np.cos(twotheta)
        sin2theta = np.sin(twotheta)
        return (1 + cos2theta * cos2theta) / sin2theta

class GIWAXS(Film):
    ylabel = "Intensity"
    xlabel = r"$q\ (\mathregular{\AA}^{-1})$"
    init_peak_height = 5e-4

    def __init__(self, a: float, c: float, q: dict, counts: dict,
                 weights: dict, sectors: dict, det_dist:float=150.0, sample_size:float=5.,
                 wavelength:float=1.54185, name: str="", background:str="new"):
        super().__init__(a, c, q, counts, False, weights, det_dist, sample_size, wavelength, name, background)
        self.giwaxs = True
        self.hex_params["w0"] = 0.
        self.hex_params["voigt"] = 0
        self.goni_lock = True
        self.keys = q.keys()
        self.sectors = dict(zip(self.keys, sectors))
        self.two_pi_over_wavelength = TWO_PI / self.wavelength
    
    def lorentz_polarization_factor(self):
        """calculate the lorentz polarization factor"""
        lorentz = {}
        for key in self.keys:
            lorentz[key] = super().lorentz_polarization_factor(key)
        return lorentz

    def initialize_fit(self, azi_buffer=20, hkl_override=None, q_buffer=1.):
        if hkl_override is None:
            hkl_override = [
                [1,0,0],
                [1,0,1],
                [2,0,0],
                [1,1,1],
                [0,0,2],
                [2,0,1],
                [1,0,2],
                [2,1,0],
                [1,1,2],
                [2,1,1],
                [2,0,2],
                [3,0,1],
                [2,1,2],
                [2,2,0],
                [3,1,0],
            ]
        super().initialize_fit(hkl_override, q_buffer, azi_buffer)
        self.show_fit()

    
    def calc_widths(self, params, theta1):
        """Return FWHM^2 wrt q-space"""
        instrumental = params["w0"]

        width_grain = TWO_PI * 0.9 / params["grain"]
        width_strain = self.two_pi_over_wavelength * params["strain"] * np.sin(theta1)
        width_size = self.two_pi_over_wavelength * params["size"] / params["dist"] * np.tan(2. * theta1) * np.cos(theta1)
        fwhm_sq = instrumental * instrumental + width_grain * width_grain + width_strain * width_strain + width_size * width_size
        return fwhm_sq
    
    def fitting_function(self, hex, hex_peak_heights, hex_params, bkgd_heights, bkgd_const, bkgd_lin, sector):
        inv_d = np.sqrt(np.sum(self.hh_hk_kk_ll[sector] / (hex * hex), axis=1, keepdims=True))       # inverse d from miller indices (column vector)
        theta1 = np.arcsin(0.5 * self.wavelength * inv_d)
        fwhm_sq = self.calc_widths(hex_params, theta1)
        sigma_sq = fwhm_sq * FWHM_SQ_TO_HALF_SIGMA_SQ
        width_gauss = 2 * sigma_sq
        width_lortz = 0.25 * fwhm_sq
        # hwhm = 0.5 * np.sqrt(fwhm_sq)
        
        q_c = TWO_PI * inv_d                                     # find q-centers from inv_d
        q_shift = self.q[sector] - q_c                                  # center q from the q-centers
        q_shift_sq = q_shift * q_shift

        hex_peaks = np.sum(
            # self.hex_multiplicity * hex_peak_heights * ((1 - hex_params["voigt"]) * np.exp(-0.5 * q_shift_sq / sigma_sq) / np.sqrt(2. * np.pi * sigma_sq) + hex_params["voigt"] * hwhm / (np.pi * (q_shift_sq / 0.25 * fwhm_sq + 1))),
            self.hex_multiplicity[sector] * hex_peak_heights * ((1 - hex_params["voigt"]) * np.exp(-q_shift_sq / width_gauss) + hex_params["voigt"] * INV_ROOT_PI_LN2 / (q_shift_sq / width_lortz + 1)),
            # self.hex_multiplicity * hex_peak_heights * ((1 - hex_params["voigt"]) * np.exp(-0.5 * q_shift_sq / fwhm_sq) + hex_params["voigt"] * fwhm_sq / (q_shift_sq + fwhm_sq)),
            # self.hex_multiplicity * hex_peak_heights * ((1 - hex_params["voigt"]) * np.exp(-0.5 * q_shift_sq / sigma_sq) + hex_params["voigt"] * fwhm_sq / (q_shift_sq + fwhm_sq)),
            axis=0
        )

        return self.lorentz[sector] * (hex_peaks + self.calc_background(bkgd_heights, sector)) + bkgd_const + bkgd_lin * self.q[sector]

    def report_peaks(self):
        for key in self.keys:
            print(f"Sector: {key}")
            q_xy_q_z = self.hh_hk_kk_ll[key] / (self.hex * self.hex)
            q = 2 * np.pi * np.sqrt(np.sum(q_xy_q_z, axis=1))
            azi = np.rad2deg(np.arctan2(*q_xy_q_z.T))
            warnings = 0
            for ii, (h, k, l) in enumerate(self.hex_hkl[key]):
                if self.hex_peak_heights[key][ii] < 1e-10:
                    warnings += 1
                    print(f"Warning: ({h}{k}{l}) at q = {q[ii]:4f}: {self.hex_peak_heights[key][ii, 0]:.5e}")
            if warnings:
                print("Remove peaks considering forbidden conditions (h+k=3n & l is odd) near the same q of these.")
                print("Use .remove_peak(h, k, l) to remove a peak. Do this above .fit_peak_heights() and re-run Notebook.")
                print("")
            for ii, (h, k, l) in enumerate(self.hex_hkl[key]):
                print(f'({h}{k}{l}) at q = {q[ii]:.4f} inv A -- {azi[ii]:.1f}\u00B0 -- {self.hex_peak_heights[key][ii, 0]:.5e}')
            print("")

    def fit_free_params(self, param_keys=[]):
        free_parameters = list(self.hex)
        for key in param_keys:
            free_parameters.append(self.hex_params[key])
        N = {}
        M = {}
        total_data = 0
        sector_start = {}
        sector_end = {}
        for sector in self.keys:
            N[sector] = len(free_parameters)
            free_parameters.extend(list(self.hex_peak_heights[sector].flatten()))
            M[sector] = len(free_parameters)
            sector_start[sector] = total_data
            total_data += len(self.counts[sector])
            sector_end[sector] = total_data
        
        def residuals(params):
            hex_lattice = np.array(params[:2])
            hex_params = self.hex_params.copy()
            for ii in range(len(param_keys)):
                hex_params[param_keys[ii]] = params[2 + ii]
            res = np.empty(total_data)
            for sector in self.keys:
                peak_heights = np.array(params[N[sector]:M[sector]]).reshape((len(self.hex_peak_heights[sector]), 1))
                counts_hat = self.fitting_function(hex_lattice, peak_heights, hex_params,
                                                   self.bkgd_heights[sector], self.bkgd_const[sector], self.bkgd_lin, sector)
                res[sector_start[sector]:sector_end[sector]] = (counts_hat - self.counts[sector]) / self.weights[sector]
            return res
        
        fit = least_squares(residuals, free_parameters, method="lm")

        try:
            std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
            print(f"a = ({fit.x[0]:.4f} \u00B1 {std[0]:.4f}) Angstroms")
            print(f"c = ({fit.x[1]:.4f} \u00B1 {std[1]:.4f}) Angstroms")
            if param_keys is not None:
                for ii, key in enumerate(param_keys):
                    print(f"{key} = ({fit.x[2 + ii]:.4f} \u00B1 {std[2 + ii]:.4f})")
                    if key == "strain":
                        if fit.x[2 + ii] < 0.5 * std[2 + ii]:
                            print("Rejecting strain")
                            self.hex_strain_lock = True
                            return None
        except:
            print("failed to calculate error")
            print(f"a = {fit.x[0]:.4f} Angstroms")
            print(f"c = {fit.x[1]:.4f} Angstroms")
            if param_keys is not None:
                for ii, key in enumerate(param_keys):
                    print(f"{key} = {fit.x[2 + ii]:.4f}")
                    if key == "strain":
                        if fit.x[2 + ii] < 1e-10:
                            print("Rejecting strain")
                            self.hex_strain_lock = True
                            return None
        self.report_fit(fit)
        self.hex = np.array(fit.x[:2])
        if param_keys is not None:
            for ii, key in enumerate(param_keys):
                self.hex_params[key] = fit.x[2 + ii]
        for sector in self.keys:
            self.hex_peak_heights[sector] = fit.x[N[sector]:M[sector]].reshape((len(self.hex_peak_heights[sector]), 1))
        self.report_peaks()
        return self.show_fit()
    
    def fit_background(self, param_keys=None):
        free_parameters = list(self.hex)
        if param_keys is not None:
            for key in param_keys:
                free_parameters.append(self.hex_params[key])
        N = {}
        M = {}
        P1 = {}
        P2 = {}
        total_data = 0
        sector_start = {}
        sector_end = {}
        for sector in self.keys:
            N[sector] = len(free_parameters)
            free_parameters.extend(list(self.hex_peak_heights[sector].flatten()))
            M[sector] = len(free_parameters)
            sector_start[sector] = total_data
            total_data += len(self.counts[sector])
            sector_end[sector] = total_data
        for sector in self.keys:
            P1[sector] = len(free_parameters)
            free_parameters.extend(list(self.bkgd_heights[sector].flatten()))
            P2[sector] = len(free_parameters)
        self.free_param_num = len(free_parameters)

        def residuals(params):
            hex_lattice = np.array(params[:2])
            hex_params = self.hex_params.copy()
            for ii in range(len(param_keys)):
                hex_params[param_keys[ii]] = params[2 + ii]
            res = np.empty(total_data)
            for sector in self.keys:
                peak_heights = np.array(params[N[sector]:M[sector]]).reshape((len(self.hex_peak_heights[sector]), 1))
                bkgd_heights = np.array(params[P1[sector]:P2[sector]]).reshape((len(self.bkgd_heights[sector]), 1))
                counts_hat = self.fitting_function(hex_lattice, peak_heights, hex_params,
                                                   bkgd_heights, self.bkgd_const[sector], self.bkgd_lin, sector)
                res[sector_start[sector]:sector_end[sector]] = (counts_hat - self.counts[sector]) / self.weights[sector]
            return res

        fit = least_squares(residuals, free_parameters, method="lm")

        try:
            std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
            print(f"a = ({fit.x[0]:.4f} \u00B1 {std[0]:.4f}) Angstroms")
            print(f"c = ({fit.x[1]:.4f} \u00B1 {std[1]:.4f}) Angstroms")
            if param_keys is not None:
                for ii, key in enumerate(param_keys):
                    print(f"{key} = ({fit.x[2 + ii]:.4f} \u00B1 {std[2 + ii]:.4f})")
                    if key == "strain":
                        if fit.x[2 + ii] < 0.5 * std[2 + ii]:
                            print("Rejecting strain")
                            self.hex_strain_lock = True
                            return None
        except:
            print("failed to calculate error")
            print(f"a = {fit.x[0]:.4f} Angstroms")
            print(f"c = {fit.x[1]:.4f} Angstroms")
            if param_keys is not None:
                for ii, key in enumerate(param_keys):
                    print(f"{key} = {fit.x[2 + ii]:.4f}")
                    if key == "strain":
                        if fit.x[2 + ii] < 1e-10:
                            print("Rejecting strain")
                            self.hex_strain_lock = True
                            return None
        self.report_fit(fit)
        self.hex = np.array(fit.x[:2])
        if param_keys is not None:
            for ii, key in enumerate(param_keys):
                self.hex_params[key] = fit.x[2 + ii]
        for sector in self.keys:
            self.hex_peak_heights[sector] = fit.x[N[sector]:M[sector]].reshape((len(self.hex_peak_heights[sector]), 1))
        self.report_peaks()
        return self.show_fit()



"""
class Hexagonal_GIWAXS_Sectors(Hexagonal_GIWAXS):
    def __init__(self, a: float, c: float,
                 q_full: np.ndarray, counts_full: np.ndarray, weights_full:np.ndarray,
                 q_oop: np.ndarray, counts_oop: np.ndarray, weights_oop:np.ndarray,
                 q_ip: np.ndarray, counts_ip: np.ndarray, weights_ip:np.ndarray,
                 q_dia: np.ndarray, counts_dia: np.ndarray, weights_dia:np.ndarray,
                 name: str, det_dist:float=150.0, wavelength:float=1.54185, background:str="new"):
        super().__init__(a=a, c=c, q=q_full, counts=counts_full, name=name, weights=weights_full, det_dist=det_dist,
                         wavelength=wavelength, background=background)
        self.q = {
            "full": q_full,
            "oop": q_oop,
            "ip": q_ip,
            "dia": q_dia
        }
        self.counts = {
            "full": counts_full,
            "oop": counts_oop,
            "ip": counts_ip,
            "dia": counts_dia
        }
        self.weights = {
            "full": weights_full,
            "oop": weights_oop,
            "ip": weights_ip,
            "dia": weights_dia
        }

        if background == "new":
            self.bkgd = {
                "full": NewBackground(),
                "oop": NewBackground(),
                "ip": NewBackground(),
                "dia": NewBackground()
            }
            self.bkgd_heights = {
                "full": np.zeros(self.bkgd.centers.shape),
                "oop": np.zeros(self.bkgd.centers.shape),
                "ip": np.zeros(self.bkgd.centers.shape),
                "dia": np.zeros(self.bkgd.centers.shape)
            }
            self.bkgd_const = {
                "full": 1,
                "oop": 1,
                "ip": 1,
                "dia": 1
            }
            self.bkgd_lin = {
                "full": 0,
                "oop": 0,
                "ip": 0,
                "dia": 0
            }

        labels = ["full", "oop", "ip", "dia"]
        for label in labels:
            print(" - `.set_background_constant(background_constant, '{}')`".format(label))
            print(" - `.set_range(q_min, q_max, '{}')`.".format(label))


    def initialize_fit(self, sample_length_in_millimeter):
        # self.sample_length_sq = sample_length_in_millimeter * 
        self.lorentz_polarization_factor()

        self.grain_size = 660
        self.voigt = 0.5
        self.strain = 0
        self.goniometer_offset = 0

        self.hkl = {
            "full": np.array([
                [1,0,0],
                [1,0,1],
                [2,0,0],
                [1,1,1],
                [0,0,2],
                [2,0,1],
                [1,0,2],
                [2,1,0],
                [1,1,2],
                [2,1,1],
                [2,0,2],
                [3,0,1],
                [2,1,2],
                [2,2,0],
                [3,1,0],
            ]),
            "oop": np.array([
                [1,0,1],
                [0,0,2],
                [1,0,2],
                [1,1,2],
                [2,0,2]
            ]),
            "ip": np.array([
                [1,0,0],
                [2,0,0],
                [2,0,1]
                [2,1,0],
                [2,1,1],
                [3,0,1],
                [2,2,0],
                [3,1,0],
            ]),
            "dia": np.array([
                [1,0,1],
                [1,1,1],
                [2,0,1],
                [2,1,0],
                [1,1,2],
                [2,1,1],
                [2,0,2],
                [3,0,1],
                [2,1,2]
            ]),
        }
        self.multiplicity = {}
        self.peak_heights = {}
        self.hh_hk_kk = {}
        self.l_sq = {}
        for sector in self.hkl.keys():
            self.multiplicity[sector] = np.empty(len(self.hkl[sector]))
            for ii, (h, k, l) in enumerate(self.hkl[sector]):
                if h + k == 0:
                    self.multiplicity[sector][ii] = 1
                elif h == k or h == 0:
                    self.multiplicity[sector][ii] = 6
                else:
                    self.multiplicity[sector][ii] = 12
            self.peak_heights[sector] = np.ones((len(self.hkl[sector]), 1)) * self.init_peak_height

            self.hh_hk_kk[sector] = 4. / 3. * (self.hkl[sector][:, 0] * self.hkl[:, 0] + self.hkl[sector][:, 0] * self.hkl[sector][:, 1] + self.hkl[sector][:, 1] * self.hkl[sector][:, 1]).reshape((self.hkl.shape[0], 1))
            self.l_sq[sector] = (self.hkl[sector][:, 2] * self.hkl[sector][:, 2]).reshape((self.hkl[sector].shape[0], 1))
        print("Initialized GIWAXS with 4 sectors")

    def fit_peak_heights(self):
        print("In next cell:")
        print(" - Fit lattice parameters and peak heights with `.fit_lattice_parameters()`.\n")
        free_parameters = []
        for sector in self.peak_heights.keys():
            free_parameters.extend(list(self.peak_heights[sector].flatten()))
        self.free_param_num = len(free_parameters)
        def residuals(params):
            peak_heights = np.array(params).reshape((len(params), 1))
            counts_hat = self.fitting_function(self.a, self.c, self.w0, peak_heights,
                                               self.grain_size, self.voigt, self.strain,
                                               self.goniometer_offset, self.det_dist,
                                               self.bkgd_heights, self.bkgd_const, self.bkgd_lin)
            return (counts_hat - self.counts) / self.weights

        fit = least_squares(residuals, free_parameters, method="lm")
        
        self.report_fit(fit)
        std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
        for ii in range(len(self.hkl)):
            print(f'({self.hkl[ii, 0]},{self.hkl[ii, 1]},{self.hkl[ii, 2]}): {fit.x[ii]:.5e} \u00B1 {std[ii]:.5e}')
        self.peak_heights = fit.x.reshape((len(fit.x), 1))
        return self.show_fit()

    def fit_lattice_parameters(self):
        print("In next cell:")
        print(" - Fit peak widths with .fit_peak_widths().\n")
        N = 2
        free_parameters = [self.a, self.c] + list(self.peak_heights.flatten())
        self.free_param_num = len(free_parameters)
        def residuals(params):
            a, c = params[:N]
            peak_heights = np.array(params[N:]).reshape((len(self.peak_heights), 1))
            counts_hat = self.fitting_function(a, c, self.w0, peak_heights,
                                               self.grain_size, self.voigt, self.strain,
                                               self.goniometer_offset, self.det_dist,
                                               self.bkgd_heights, self.bkgd_const, self.bkgd_lin)
            return (counts_hat - self.counts) / self.weights
        
        fit = least_squares(residuals, free_parameters, method="lm")

        try:
            std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
            print(f"a = ({fit.x[0]:.4f} \u00B1 {std[0]:.4f}) Angstroms")
            print(f"c = ({fit.x[1]:.4f} \u00B1 {std[1]:.4f}) Angstroms")
        except:
            print("failed to calculate error")
            print(f"a = {fit.x[0]:.4f} Angstroms")
            print(f"c = {fit.x[1]:.4f} Angstroms")
        
        self.report_fit(fit)
        self.a, self.c = fit.x[:N]
        self.peak_heights = fit.x[N:].reshape((len(self.peak_heights), 1))
        self.report_peaks()
        return self.show_fit()

    def fit_peak_widths(self):
        print("In next cell:")
        print(" - Fit peak shapes with `.fit_strain()`.\n")
        N = 3
        free_parameters = [self.a, self.c, self.grain_size] + list(self.peak_heights.flatten())
        self.free_param_num = len(free_parameters)
        def residuals(params):
            a, c, grain_size = params[:N]
            peak_heights = np.array(params[N:]).reshape((len(self.peak_heights), 1))
            counts_hat = self.fitting_function(a, c, self.w0, peak_heights,
                                               grain_size, self.voigt, self.strain,
                                               self.goniometer_offset, self.det_dist,
                                               self.bkgd_heights, self.bkgd_const, self.bkgd_lin)
            return (counts_hat - self.counts) / self.weights
        
        fit = least_squares(residuals, free_parameters, method="lm")

        try:
            std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
            print(f"a = ({fit.x[0]:.4f} \u00B1 {std[0]:.4f}) Angstroms")
            print(f"c = ({fit.x[1]:.4f} \u00B1 {std[1]:.4f}) Angstroms")
            print(f"grain size = ({fit.x[2]:.4f} \u00B1 {std[2]:.4f})")
        except:
            print("failed to calculate error")
            print(f"a = {fit.x[0]:.4f} Angstroms")
            print(f"c = {fit.x[1]:.4f} Angstroms")
            print(f"grain size = {fit.x[2]:.4f}")

        self.report_fit(fit)
        self.a, self.c, self.grain_size = fit.x[:N]
        self.peak_heights = fit.x[N:].reshape((len(self.peak_heights), 1))
        self.report_peaks()
        return self.show_fit()
    
    def fit_strain(self):
        self.strain_lock = False
        print("In next cell:")
        print(" - Fit voigt parameter with `.fit_voigt()`.\n")
        # free_parameters = [self.a, self.c, self.grain_size, self.strain] + list(self.peak_heights.flatten())
        N = 4
        free_parameters = [self.a, self.c, self.grain_size, self.strain] + list(self.peak_heights.flatten())
        self.free_param_num = len(free_parameters)
        def residuals(params):
            # a, c, grain_size, strain = params[:4]
            a, c, grain_size, strain = params[:N]
            peak_heights = np.array(params[N:]).reshape((len(self.peak_heights), 1))
            counts_hat = self.fitting_function(a, c, self.w0, peak_heights,
                                               grain_size, self.voigt, strain,
                                               # grain_size, self.voigt, strain,
                                               self.goniometer_offset, self.det_dist,
                                               self.bkgd_heights, self.bkgd_const, self.bkgd_lin)
            return (counts_hat - self.counts) / self.weights
        
        fit = least_squares(residuals, free_parameters, method="lm")

        try:
            std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
            print(f"a = ({fit.x[0]:.4f} \u00B1 {std[0]:.4f}) Angstroms")
            print(f"c = ({fit.x[1]:.4f} \u00B1 {std[1]:.4f}) Angstroms")
            print(f"grain size = ({fit.x[2]:.4f} \u00B1 {std[2]:.4f})")
            print(f"strain = ({fit.x[3]:.6f} \u00B1 {std[3]:.6f})")
        except:
            print("errors failed to be calculated")
            print(f"a = {fit.x[0]:.4f} Angstroms")
            print(f"c = {fit.x[1]:.4f} Angstroms")
            print(f"grain size = {fit.x[2]:.4f} Angstroms")
            print(f"strain = {fit.x[3]:.6f}")

        if fit.x[3] < std[3]:
            print("Strain is negative or less than its uncertainty. Rejecting this fit and locking strain to 0.")
            self.strain_lock = True
        else:
            self.report_fit(fit)
            self.a, self.c, self.grain_size, self.strain = fit.x[:N]
            self.peak_heights = fit.x[N:].reshape((len(self.peak_heights), 1))
            self.report_peaks()
            return self.show_fit()
    
    def fit_voigt(self):
        print("In next cell:")
        print(" - Fit everything with the background free `.fit_full()`.\n")
        free_parameters = [self.a, self.c, self.grain_size, self.voigt]

        if not self.strain_lock:
            free_parameters.append(self.strain)
        N = len(free_parameters)

        free_parameters.extend(list(self.peak_heights.flatten()))
        self.free_param_num = len(free_parameters)
        if self.strain_lock:
            def residuals(params):
                a, c, grain_size, voigt = params[:N]
                peak_heights = np.array(params[N:]).reshape((len(self.peak_heights), 1))
                counts_hat = self.fitting_function(a, c, self.w0, peak_heights,
                                                   grain_size, voigt, 0,
                                                   self.goniometer_offset, self.det_dist,
                                                   self.bkgd_heights, self.bkgd_const, self.bkgd_lin)
                return (counts_hat - self.counts) / self.weights
        else:
            def residuals(params):
                a, c, grain_size, voigt, strain = params[:N]
                peak_heights = np.array(params[N:]).reshape((len(self.peak_heights), 1))
                counts_hat = self.fitting_function(a, c, self.w0, peak_heights,
                                                   grain_size, voigt, strain,
                                                   self.goniometer_offset, self.det_dist,
                                                   self.bkgd_heights, self.bkgd_const, self.bkgd_lin)
                return (counts_hat - self.counts) / self.weights
        
        fit = least_squares(residuals, free_parameters, method="lm", ftol=1e-9, xtol=1e-9, gtol=1e-9)
        
        try:
            std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
            print(f"a = ({fit.x[0]:.4f} \u00B1 {std[0]:.4f}) Angstroms")
            print(f"c = ({fit.x[1]:.4f} \u00B1 {std[1]:.4f}) Angstroms")
            print(f"grain size = ({fit.x[2]:.4f} \u00B1 {std[2]:.4f})")
            print(f"voigt = ({fit.x[3]:.4f} \u00B1 {std[3]:.4f})")
            if self.strain_lock:
                print(f"strain locked at 0")
            else:
                print(f"strain = ({fit.x[4]:.4f} \u00B1 {std[4]:.4f})")
        except:
            print(f"a = {fit.x[0]:.4f} Angstroms")
            print(f"c = {fit.x[1]:.4f} Angstroms")
            print(f"grain size = {fit.x[2]:.4f}")
            print(f"voigt = {fit.x[3]:.4f}")
            if self.strain_lock:
                print(f"strain locked at 0")
            else:
                print(f"strain = {fit.x[4]:.4f}")


        if not (fit.x[3] >= 0 and fit.x[3] <= 1):
            print("voigt parameter outside of bounds. Rejecting fit and setting voigt to 0.5.")
            self.voigt = 0.5
        else:
            self.report_fit(fit)
            self.a, self.c, self.grain_size, self.voigt = fit.x[:N]
            if not self.strain_lock:
                self.strain = fit.x[4]
            self.peak_heights = fit.x[N:].reshape((len(self.peak_heights), 1))
            self.report_peaks()
            return self.show_fit()
    
    def fit_full(self):
        free_parameters = [self.a, self.c, self.grain_size, self.voigt]
        if not self.strain_lock:
            free_parameters.append(self.strain)
        N = len(free_parameters)
        free_parameters.extend(list(self.peak_heights.flatten()))
        free_parameters.extend(list(self.bkgd_heights.flatten()))
        free_parameters.append(self.bkgd_const)
        self.free_param_num = len(free_parameters)

        if self.strain_lock:
            def residuals(params):
                a, c, grain_size, voigt = params[:N]
                peak_heights = np.array(params[N:N+len(self.peak_heights)]).reshape((len(self.peak_heights), 1))
                bkgd_heights = np.array(params[N+len(self.peak_heights):N+len(self.peak_heights)+len(self.bkgd_heights)]).reshape((len(self.bkgd_heights), 1))
                bkgd_const = params[N+len(self.peak_heights)+len(self.bkgd_heights)]
                counts_hat = self.fitting_function(a, c, self.w0, peak_heights,
                                                   grain_size, voigt, 0,
                                                   self.goniometer_offset, self.det_dist,
                                                   bkgd_heights, bkgd_const, self.bkgd_lin)
                return (counts_hat - self.counts) / self.weights
        else:
            def residuals(params):
                a, c, grain_size, voigt, strain = params[:N]
                peak_heights = np.array(params[N:N+len(self.peak_heights)]).reshape((len(self.peak_heights), 1))
                bkgd_heights = np.array(params[N+len(self.peak_heights):N+len(self.peak_heights)+len(self.bkgd_heights)]).reshape((len(self.bkgd_heights), 1))
                bkgd_const = params[N+len(self.peak_heights)+len(self.bkgd_heights)]
                counts_hat = self.fitting_function(a, c, self.w0, peak_heights,
                                                   grain_size, voigt, strain,
                                                   self.goniometer_offset, self.det_dist,
                                                   bkgd_heights, bkgd_const, self.bkgd_lin)
                return (counts_hat - self.counts) / self.weights

        fit = least_squares(residuals, free_parameters, method="lm", ftol=1e-9, xtol=1e-9, gtol=1e-9)

        std = np.sqrt(np.diagonal(np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun.T @ fit.fun / (fit.fun.size - fit.x.size))))
        print(f"a = ({fit.x[0]:.4f} \u00B1 {std[0]:.4f}) Angstroms")
        print(f"c = ({fit.x[1]:.4f} \u00B1 {std[1]:.4f}) Angstroms")
        print(f"grain size = ({fit.x[2]:.4f} \u00B1 {std[2]:.4f})")
        print(f"voigt = ({fit.x[3]:.4f} \u00B1 {std[3]:.4f})")
        if self.strain_lock:
            print("strain locked at 0")
        else:
            print(f"strain = ({fit.x[4]:.4f} \u00B1 {std[4]:.4f})")
        for ii in range(len(self.hkl)):
            print(f'({self.hkl[ii, 0]},{self.hkl[ii, 1]},{self.hkl[ii, 2]}): {fit.x[ii+N]:.5e} \u00B1 {std[ii+4]:5e}')
        for ii in range(len(self.bkgd_heights)):
            print(f"bkgd peak {ii}: {fit.x[ii+N+len(self.hkl)]:.5e} \u00B1 {std[ii+N+len(self.hkl)]:.5e}")
        print(f"background constant = {fit.x[N+len(self.hkl)+len(self.bkgd_heights)]:.5e} \u00B1 {std[N+len(self.hkl)+len(self.bkgd_heights)]:5e}")

        self.report_fit(fit)
        self.a, self.c, self.grain_size, self.voigt = fit.x[:4]
        if not self.strain_lock:
            self.strain = fit.x[4]

        self.peak_heights = fit.x[N:N+len(self.peak_heights)].reshape((len(self.peak_heights), 1))
        self.bkgd_heights = fit.x[N+len(self.peak_heights):N+len(self.peak_heights)+len(self.bkgd_heights)].reshape((len(self.bkgd_heights), 1))
        self.bkgd_const = fit.x[N+len(self.peak_heights)+len(self.bkgd_heights)]
        
        self.report_peaks()

        print(".save(file_name, title=None, dpi=900)")
        return self.show_fit()
    
    def fitting_function(self, a, c, w0, peak_heights, grain_size, voigt, strain, goniometer_offset, det_dist, bkgd_heights, bkgd_const, bkgd_lin):
        result = np.array([])
        for sector in self.q.keys():
            inv_d = np.sqrt(self.hh_hk_kk[sector] / (a * a) + self.l_sq[sector] / (c * c))              # inverse d from miller indices
            theta1 = np.arcsin(0.5 * self.wavelength * inv_d)                           # 2theta corresponding to inv_d
            # width_grain = 0.9 * self.wavelength / (grain_size * np.cos(theta1))         # width due to grain size
            # width_strain = strain * np.tan(theta1)                                      # width due to strain
            geometric = 0.5 * self.fourpioverlambda * np.cos(theta1)
            tan2theta = np.tan(2. * self.theta)
            wq_sq = geometric * geometric * (w0 * w0 + self.sample * tan2theta * tan2theta )   # combine widths in quadrature
            q_c = 2 * np.pi * inv_d                                     # find q-centers from inv_d
            q_shift = self.q[sector] - q_c                                      # center q from the q-centers
            arg_sq = q_shift * q_shift / wq_sq

            peaks = np.sum(
                self.multiplicity[sector] * peak_heights[sector] * ((1 - voigt) * np.exp(-0.5 * arg_sq) + (voigt / (arg_sq + 1))),
                axis=0
            )
            result = np.concatenate((result, self.lorentz[sector] * (peaks + self.calc_background(bkgd_heights[sector])) + bkgd_const[sector] + bkgd_lin[sector] * self.q[sector]))
        return result

def q_from_2th(two_theta, wavelength=1.54185):
    return 4 * np.pi * np.sin(0.5 * np.radians(two_theta)) / wavelength



# def fitting_function(self, a, c, w0, peak_heights, grain_size, voigt, strain, goniometer_offset, det_dist, bkgd_heights, bkgd_const, bkgd_lin):
#         inv_d = np.sqrt(self.hh_hk_kk / (a * a) + self.l_sq / (c * c))              # inverse d from miller indices
#         theta1 = np.arcsin(0.5 * self.wavelength * inv_d)                           # 2theta corresponding to inv_d
#         width_grain = 0.9 * self.wavelength / (grain_size * np.cos(theta1))         # width due to grain size
#         width_strain = strain * np.tan(theta1)                                      # width due to strain
#         wq_sq = w0 * w0 + width_grain * width_grain + width_strain * width_strain   # combine widths in quadrature
#         q_c = 4 * np.pi / self.wavelength * np.sin(theta1)          # find q-centers from 2theta (this could be calculated more simply from inv_d)
#         q_shift = self.q - q_c                                      # center q from the q-centers
#         arg_sq = 0.5 * q_shift * q_shift / wq_sq                    
# 
#         peaks = np.sum(
#             self.multiplicity * peak_heights * ((1 - voigt) * np.exp(-arg_sq) + (voigt * ROOT_LN2__PI / (arg_sq + LN2))),
#             axis=0
#         )
#         return self.lorentz * (peaks + self.calc_background(bkgd_heights)) + bkgd_const + bkgd_lin * self.q"""