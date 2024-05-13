import numpy as np
import argparse
import yaml
import fabio
import pyFAI
from pathlib import Path
import XRDpy.package_params as package
import shutil
from XRDpy.tiff_loader import load_from
from XRDpy.transform import TransformGIX
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

class ProcessParse:

    FIG_SIZE = (7.5, 3.75)

    PARAMS_FILE = "params.yaml"
    
    FLAGS_VAL = [
        ('-E', '--exposure'),
        ('-I', '--incident'),
        ('-T', '--tilt'),
        ('-D', '--dir'),
    ]
    FLAGS_BOOl = [
        ('-P', '--plot'),
        ('-O', '--override'),
        ('-C', '--cake'),
        ('-R', '--reduce'),
    ]

    def __init__(self, program_name, description):
        self.parser = argparse.ArgumentParser(prog=program_name, description=description, epilog="Author: Teddy Tortorici edward.tortorici@colorado.edu")
        for flag in self.FLAGS_VAL:
            self.parser.add_argument(*flag)
        for flag in self.FLAGS_BOOl:
            self.parser.add_argument(*flag, action='store_true')
        self.dir = None
        self.exposure = None
        self.incident = None
        self.tilt = None

    def stitch(self):
        self.args = self.parser.parse_args()

        if self.args.dir is None:
            self.dir = Path.cwd()
        else:
            self.dir = Path(self.args.dir)
        if self.args.exposure is not None:
            self.exposure = int(self.args.exposure)
        else:
            self.exposure = None
        if self.args.incident is not None:
            self.incident = float(self.args.incident)
        else:
            self.incident = None
        if self.args.tilt is not None:
            self.tilt = float(self.args.tilt)
        else:
            self.tilt = None

        if self.args.override:
            stitched_data_file = self.dir / "raw-stitched-data.tif"
            flat_field_file = self.dir / "flat-field.tif"
            if stitched_data_file.is_file():
                stitched_data_file.unlink()
            if flat_field_file.is_file():
                flat_field_file.unlink()

        self.load_params()
        data, flat_field = load_from(self.dir, self.exposure)

        if self.args.plot:
            adjuster = self.exposure / flat_field
            adjuster[np.where(adjuster == np.infty)] = 0.0
            data_adj = data * adjuster
            fig = plt.figure(figsize=self.FIG_SIZE, facecolor="w")
            ax1 = plt.subplot()
            pos = ax1.imshow(data_adj+1, norm=LogNorm(1, np.max(data_adj)))
            ax1.set_title("Stitched image")
            ax1.set_xlabel("column (pixels)")
            ax1.set_ylabel("row (pixels)")
            fig.colorbar(pos, ax=ax1, shrink=0.7)

        return data, flat_field
        
    def make_params(self):
        params = {"exposure": self.exposure,
                  "incident": self.incident,
                  "tilt": self.tilt}
        return params

    def save_tiff(self, image, filename):
        fabio.tifimage.tifimage(image).write(self.dir / filename)

    def save_yaml(self, data, filename):
        with open(self.dir / filename, "w") as f:
            yaml.dump(data, f)

    def load_params(self):
        exposure_load, incident_load, tilt_load = False, False, False
        if self.exposure is None:
            exposure_load = True
        if self.args.incident is None:
            incident_load = True
        if self.args.tilt is None:
            tilt_load = True

        params_file = self.dir / self.PARAMS_FILE

        if params_file.is_file():
            with open(params_file, "r") as f:
                params = yaml.safe_load(f)
            if not isinstance(params, dict):
                params = self.make_params()
            if "exposure" not in params:
                params["exposure"] = self.exposure
            if "incident" not in params:
                params["incident"] = self.incident
            if "tilt" not in params:
                params["tilt"] = self.tilt
        else:
            params = self.make_params()

        if exposure_load:
            self.exposure = params["exposure"]
        if incident_load:
            self.incident = params["incident"]
        if tilt_load:
            self.tilt = params["tilt"]

        self.save_yaml(self.make_params(), params_file)
    
    def copy_notebook(self, notebook):
        if (self.dir / notebook).is_file():
            print(f'File {notebook} already exists in "{self.dir}"')
        elif not (package.directory / "Jupyter-Notebooks" / notebook).is_file():
            print(f'File {notebook} does not exist in "{package.directory / 'Jupyter-Notebooks'}"')
        else:
            shutil.copyfile(package.directory / "Jupyter-Notebooks" / notebook, self.dir / notebook)
            print(f'Copied {notebook} to "{self.dir}"')
        

class StitchParse(ProcessParse):
    FLAGS_VAL = ProcessParse.FLAGS_VAL[1:]
    FLAGS_BOOl = ProcessParse.FLAGS_BOOl[:2]

    def __init__(self):
        super().__init__(
            "XRD-stitch",
            "Stitch together raw images from Eiger R 1M that eiger-stitch produced. These should all be together in the same directory with names X_X_X",
        )
        self.parser.add_argument("exposure")
        self.stitch()
        self.copy_notebook("reduce-WAXS.ipynb")
        


class FilmParse(ProcessParse):
    def __init__(self):
        super().__init__(
            "XRD-film",
            "Transform stitched images into reciprocal space. Assumes the sample is an aligned thin film. Requires an incident angle, but this can be given through 'params.yaml'"
        )
        data, flat_field = self.stitch()

        data_t_file = self.dir / "image-transformed.tif"
        flat_t_file = self.dir / "flat-field-transformed.tif"
        cal_t_file = self.dir / "cal-transformed.poni"

        if self.args.override or not (data_t_file.is_file() and
                                      flat_t_file.is_file() and
                                      cal_t_file.is_file()):
            transformer = TransformGIX(self.incident, self.tilt)
            transformer.load(self.dir / "cal.poni")
            # print("Start image transform")
            (data_t, flat_t), beam_center_t = transformer.transform_image(data, flat_field)

            ai = pyFAI.load(self.dir / "cal.poni")
            pixel1 = ai.get_pixel1()
            pixel2 = ai.get_pixel2()
            ai.poni1 = (data_t.shape[0] - beam_center_t[0]) * pixel1
            ai.poni2 = beam_center_t[1] * pixel2
            ai.detector = pyFAI.detectors.Detector(pixel1=pixel1, pixel2=pixel2, max_shape=data_t.shape, orientation=2)
            ai.save(cal_t_file)
            self.save_tiff(data_t, "image-transformed.tif")
            self.save_tiff(flat_t, "flat-field-transformed.tif")
        else:
            print("Found preexisting transform to load. To override, run again with --override or -O")
            data_t = fabio.open(self.dir / "image-transformed.tif")
            flat_t = fabio.open(self.dir / "flat-field-transformed.tif")
            ai = pyFAI.load(cal_t_file)

        if self.args.cake or self.args.reduce or self.args.plot:
            adjuster = self.exposure / flat_t
            adjuster[np.where(adjuster == np.infty)] = 0.0
            data_adj = data_t * adjuster

        if self.args.plot:
            fig = plt.figure(figsize=self.FIG_SIZE, facecolor="w")
            ax1 = plt.subplot()
            pos = ax1.imshow(data_adj+1, norm=LogNorm(1, np.max(data_adj)))
            ax1.set_title("Transformed image")
            ax1.set_xlabel("column (pixels)")
            ax1.set_ylabel("row (pixels)")
            fig.colorbar(pos, ax=ax1, shrink=0.7)
        
        if self.args.cake or self.args.reduce:
            mask = np.logical_not(flat_t)
        
        if self.args.cake:
            cake_file = self.dir / "cake.edf"
            if cake_file.is_file():
                cake_file.unlink()
            res2 = ai.integrate2d_ng(data_adj, 1000, 1800, mask=mask, unit="q_A^-1",
                                     filename=str(cake_file))
            if self.args.plot:
                fig = plt.figure(figsize=self.FIG_SIZE, facecolor="w")
                ax1 = plt.subplot()
                pos = ax1.imshow(res2[0]+1, norm=LogNorm(1, np.max(res2[0])),
                                 extent=(np.min(res2[1]), np.max(res2[1]), np.min(res2[2]), np.max(res2[2])),
                                 aspect='auto')
                ax1.set_title("Caked image")
                ax1.set_xlabel(r"$q\ (\mathregular{\AA}^{-1})$")
                ax1.set_ylabel(r"$\Omega$ (degrees)")
                fig.colorbar(pos, ax=ax1, shrink=0.7)
        if self.args.reduce:
            reduced_file = self.dir / "1D-reduction.dat"
            if reduced_file.is_file():
                reduced_file.unlink()
            res = ai.integrate1d_ng(data_adj, 1000, mask=mask, unit="q_A^-1",
                                    filename=str(reduced_file))
            if self.args.plot:
                plt.figure()
                plt.plot(res[0], res[1])
                plt.title("1D reduction")
                plt.xlabel(r"$q\ (\mathregular{\AA}^{-1})$")
                plt.ylabel("Intensity")
        self.copy_notebook("reduce-GIWAXS.ipynb")
        
        
    
# class DefaultParse(ProcessParse):
#     # FLAGS_BOOl = ProcessParse.FLAGS_BOOl[-1]
#     user_config_dir = Path.home() / "Documents" / package.name
#     
#     def __init__(self):
#         super().__init__("XRD-default", "Change user defaults")
#         self.args = self.parser.parse_args()
# 
# 
#     def save_user_config(self):
#         print("Saving new config file: {}".format(self.user_config_dir / package.name))
#         if not self.user_config_dir.exists():
#             self.user_config_dir.mkdir()
#         with open(self.user_config_dir / package.config_name, 'w') as f:
#             yaml.dump(self.config, f)
# 
#     @classmethod
#     def load(cls):
#         user_config = cls.user_config_dir / package.config_name
#         if user_config.is_file():
#             with open(user_config, 'r') as f:
#                 config = yaml.safe_load(f)
#         else:
#             config = DefaultParse.load_default()
#         return config
#     
#     @staticmethod
#     def load_default():
#         config = package.default_config
#         return config
#     
#     def start(self):
#         while True:
#             try:
#                 self.next()
#             except StopIteration:
#                 break
#         self.save_user_config()
        

def stitch():
    parser = StitchParse()
    if parser.args.plot:
        plt.show()

def film():
    parser = FilmParse()
    if parser.args.plot:
        plt.show()

# def default():
#     parser = DefaultParse()