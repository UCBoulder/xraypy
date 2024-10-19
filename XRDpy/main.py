import numpy as np
from datetime import datetime
import argparse
import yaml
import fabio
import pyFAI
from pathlib import Path
import XRDpy.package_params as package
import shutil
from XRDpy.tiff_loader import load_from, Stitcher
from XRDpy.transform import TransformGIX
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

DATA_PATH = Path("DATA")

class ProcessParse:

    FIG_SIZE = (7.5, 3.75)

    PARAMS_FILE = "params.yaml"
    
    FLAGS_VAL = [
        # ('-E', '--exposure'),
        # ('-I', '--incident'),
        ('-T', '--tilt'),
        ('-D', '--dir'),
        ('-J', '--jupyter-notebook'),
        ('-F', '--file')
    ]
    FLAGS_BOOl = [
        # ('-P', '--plot'),
        # ('-O', '--override'),
        # ('-C', '--cake'),
        # ('-R', '--reduce'),
    ]

    def __init__(self, program_name: str, description: str):
        self.parser = argparse.ArgumentParser(prog=program_name,
                                              description=description,
                                              epilog="Author: Teddy Tortorici edward.tortorici@colorado.edu")
        for flag in self.FLAGS_VAL:
            self.parser.add_argument(*flag)
        for flag in self.FLAGS_BOOl:
            self.parser.add_argument(*flag, action='store_true')
        self.dir = None
        self.exposure = None
        self.incident = None
        self.tilt = None
        self.waveguiding = None
        self.jupyter_notebook = None
        self.file = None

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
        try:
            if self.args.incident is not None:
                self.incident = float(self.args.incident)
            else:
                self.incident = None
            if self.args.tilt is not None:
                self.tilt = float(self.args.tilt)
            else:
                self.tilt = None
            if self.args.waveguide is not None:
                self.waveguiding = self.args.waveguide
            else:
                self.waveguiding = None
        except AttributeError:
            self.incident = None
            self.tilt = None
            self.waveguiding = None
        if self.args.jupyter_notebook is not None:
            self.jupyter_notebook = self.args.jupyter_notebook.replace("'", "").replace('"', '').split(",")
            for ii, notebook in enumerate(self.jupyter_notebook):
                if not notebook.endswith(".ipynb"):
                    self.jupyter_notebook[ii] = notebook.strip(" ") + ".ipynb"
        if self.args.file is not None:
            self.file = self.args.file.replace("'", "").replace('"', '').split(",")
            for ii, file in enumerate(self.file):
                self.file[ii] = file.strip(" ")

        stitched_data_file = self.dir / "raw-stitched-data.tif"
        flat_field_file = self.dir / "flat-field.tif"
        if stitched_data_file.is_file():
            stitched_data_file.unlink()
        if flat_field_file.is_file():
            flat_field_file.unlink()

        self.save_params()
        data, flat_field = load_from(self.dir, self.exposure)

        return data, flat_field
        
    def make_params(self):
        params = {"exposure": self.exposure,
                  "incident": self.incident,
                  "tilt": self.tilt}
        return params

    def save_tiff(self, image: np.ndarray, filename: str):
        fabio.tifimage.tifimage(image).write(self.dir / filename)

    def save_yaml(self, data, filename: str):
        with open(self.dir / filename, "w") as f:
            yaml.dump(data, f)

    def save_params(self):
        params_file = self.dir / self.PARAMS_FILE
        self.save_yaml(self.make_params(), params_file)
    
    def copy_notebook(self, notebook: str):
        if (self.dir / notebook).is_file():
            print(f'File {notebook} already exists in "{self.dir}"')
        elif not (package.directory / "Jupyter-Notebooks" / notebook).is_file():
            directory = package.directory / "Jupyter-Notebooks"
            print(f'File {notebook} does not exist in "{directory}"')
        else:
            shutil.copyfile(package.directory / "Jupyter-Notebooks" / notebook, self.dir / notebook)
            print(f'Copied {notebook} to "{self.dir}"')
    
    def copy_file(self, filename: str):
        if (self.dir / filename).is_file():
            print(f'File {filename} already exists in "{self.dir}"')
        elif not (package.directory / "User-files" / filename).is_file():
            directory = package.directory / 'User-files'
            print(f'File {filename} does not exist in "{directory}"')
        else:
            shutil.copyfile(package.directory / "User-files" / filename, self.dir / filename)
            print(f'Copied {filename} to "{self.dir}"')

    def copy_notebooks(self):
        if self.args.jupyter_notebook is not None:
            for notebook in self.args.jupyter_notebook:
                self.copy_notebook(notebook)
    
    def copy_files(self):
        if self.args.file is not None:
            for file in self.args.file:
                self.copy_file(file)
        

class StitchParse(ProcessParse):
    FLAGS_VAL = ProcessParse.FLAGS_VAL[1:]
    # FLAGS_BOOl = ProcessParse.FLAGS_BOOl[:2]

    def __init__(self):
        super().__init__(
            "WAXS/SAXS",
            "Stitch together raw images from Eiger R 1M that eiger-stitch produced. These should all be together in the same directory with names X_X_X",
        )
        self.parser.add_argument("exposure")
        self.stitch()
        self.copy_notebook("reduce-WAXS.ipynb")
        self.copy_notebooks()
        self.copy_files()


class StitchParse2(ProcessParse):
    FLAGS_VAL = ProcessParse.FLAGS_VAL[1:]
    DATA = ""

    def __init__(self):
        super().__init__(
            "stitch",
            "Stitch together raw images from Eiger R 1M that eiger-stitch produced.",
        )
        self.parser.add_argument("experiment", help="experiment type. Example: WAXS, SAXS, GIWAXS, GISAXS")
        self.parser.add_argument("rows")
        self.parser.add_argument("columns")
        self.parser.add_argument("exposure")
        self.parser.add_argument("-I", "--incident", help="incident angle in degrees (for GIXS)")
        self.parser.add_argument("-O", "--om", help="omega motor position in degrees (for GIXS)")
        self.parser.add_argument("-D", "--dist", help="sample detector distance")
        self.parser.add_argument("-T", "--tif", action="store_true")
        self.stitch()
    
    def stitch(self, name_override=None):
        self.args = self.parser.parse_args()

        if self.args.dir is None:
            self.dir = Path.cwd()
        else:
            self.dir = Path(self.DATA)
        
        if name_override is None:
            newest_file = self.get_newest_file(self.dir)
            filename_base = newest_file.name.strip(".tif")
        else:
            filename_base = name_override.strip(".tif").strip(".edf")

        stitcher = Stitcher(int(self.args.rows), int(self.args.columns))
        data, flat_field = stitcher.load_data(DATA_PATH)

        date = datetime.now()
        header = {
            "StitchingRows": self.args.rows,
            "StitchingColumns": self.args.cols,
            "ExposureTime(s)": self.args.exposure,
            "ExperimentType": self.args.experiment,
            "IncidentAngle": self.args.incident,
            "OmegaMotorAngle": self.args.om,
            "DetectorDistance": self.args.dist,
            "StitchDate": f"{date.year}-{date.month:02d}-{date.day:02d}"
        }

        edf_data_obj = fabio.edfimage.EdfImage(data=data, header=header)
        edf_flat_obj = fabio.edfimage.EdfImage(data=flat_field, header=header)
        edf_data_obj.write(self.dir / (filename_base + "_data.edf"))
        edf_flat_obj.write(self.dir / (filename_base + "_flat-field.edf"))

        if self.args.tif:
            data_im = fabio.tifimage.tifimage(data)
            flat_field = fabio.tifimage.tifimage(flat_field)
            data_im.write(self.dir / (filename_base + "_data.tif"))
            flat_field.write(self.dir / (filename_base + "_flat-field.tif"))

    @staticmethod
    def get_newest_file(directory: Path) -> Path:
        files = list(directory.glob("*.tif"))
        return max(files, key=lambda f: f.stat().st_mtime)


class FilmParse(ProcessParse):
    def __init__(self):
        super().__init__(
            "GIWAXS/GISAXS",
            "Transform stitched images into reciprocal space. Assumes the sample is an aligned thin film. Requires an incident angle."
        )
        self.parser.add_argument("exposure")
        self.parser.add_argument("incident")
        self.parser.add_argument('-W', '--waveguide', action='store_true')
        data, flat_field = self.stitch()

        data_t_filename = "image-transformed.tif"
        flat_t_filename = "flat-field-transformed.tif"
        cal_t_filename = "cal-transformed.poni"

        transformer = TransformGIX(self.incident, self.tilt, self.waveguiding)
        transformer.load(self.dir / "cal.poni")
        # print("Start image transform")
        (data_t, flat_t), beam_center_t = transformer.transform_image(data, flat_field)
        
        ai = pyFAI.load(self.dir / "cal.poni")
        pixel1 = ai.get_pixel1()
        pixel2 = ai.get_pixel2()
        ai.poni1 = (data_t.shape[0] - beam_center_t[0]) * pixel1
        ai.poni2 = beam_center_t[1] * pixel2
        ai.detector = pyFAI.detectors.Detector(pixel1=pixel1, pixel2=pixel2, max_shape=data_t.shape, orientation=2)
        ai.save(self.dir / cal_t_filename)
        self.save_tiff(data_t, data_t_filename)
        self.save_tiff(flat_t, flat_t_filename)

        self.copy_notebook("reduce-GIWAXS.ipynb")
        self.copy_notebooks()
        self.copy_files()
        

def stitch():
    parser = StitchParse()


def stitch2():
    parser = StitchParse2()
    

def film():
    parser = FilmParse()
