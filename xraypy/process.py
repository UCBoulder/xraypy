import sys
import yaml
import fabio
import matplotlib.pylab as plt
import matplotlib
from pathlib import Path
from XRDpy.data_loading import load_from
from XRDpy.process_default import ConfigHandler
from XRDpy.transform import TransformGIX
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

class Handler:
    def __init__(self, to_do: str, args: list) -> None:
        self.to_do = to_do
        self.config = ConfigHandler.load()
        self.iterator = iter(args)
        self.flag = None
        self.item = None
        self.directory = Path.cwd()
        self.process = self.config["process"].lower()
        self.incident_angle = self.config["incident"]
        self.incident_load = True
        self.tilt_angle = self.config["tilt"]
        self.tilt_load = True
        self.exposure_time = self.config["exposure"]
        self.exposure_load = True
        self.plot = self.config["plot"]
        self.parse_flags()
        self.start()
    
    def make_params(self):
        params = {"exposure": self.exposure_time,
                  "incident": self.incident_angle,
                  "tilt": self.tilt_angle}
        return params

    def start(self):
        params_file = self.directory / "params.yaml"
        if params_file.is_file():
            with open(params_file, "r") as f:
                params = yaml.safe_load(f)
            if not isinstance(params, dict):
                params = self.make_params()
            if "exposure" not in params:
                print(params)
                params["exposure"] = self.exposure_time
                print(params)
            if "incident" not in params:
                params["incident"] = self.incident_angle
            if "tilt" not in params:
                params["tilt"] = self.tilt_angle
        else:
            params = self.make_params()
        if self.exposure_load:
            self.exposure_time = params["exposure"]
        if self.exposure_time is None:
            raise ValueError("You must give an expsoure time (this is the same number you give eiger_stitch or eiger_run). If a 'params.yaml' does not exist in this directory, you can set the exposure time (s) with `--exposure=x`.")
        data, weight = load_from(self.directory, self.exposure_time)
        
        if self.to_do == "stitch":
            with open(params_file, "w") as f:
                yaml.dump(self.make_params(), f)
            return data, weight
        else:
            if not (self.directory / "cal.poni").is_file():
                raise FileExistsError("You must have a PONI file titled 'cal.poni' in your working directory.\nYou can generate this with `pyFAI-calib2`.")
        
        if self.to_do == "gix":
            if self.incident_load:
                self.incident_angle = params["incident"]
            if self.tilt_load:
                self.tilt_angle = params["tilt"]
            if self.incident_angle is None:
                raise ValueError("You must give an incident angle. If a 'params.yaml' does not exist in this directory, you can set incident angle and tilt angle with `--incident=x` and `--tilt=y`")
            with open(params_file, "w") as f:
                yaml.dump(self.make_params(), f)

            transformer = TransformGIX(self.incident_angle, self.tilt_angle)
            transformer.load(self.directory / "cal.poni")
            print("Start image transform")

            data_t, weight_t, extent_t = transformer.transform_image(data, self.exposure_time, weight,
                                                                 scale=0.9, plot=self.plot)
            self.save_tiff(data_t, "transformed_data.tif")
            self.save_tiff(weight_t, "transformed_image_weights.tif")
            self.save_yaml(data={"qxy": (extent_t[0], extent_t[1], data_t.shape[1]),
                                 "qy": (extent_t[2], extent_t[3], data_t.shape[0])},
                           filename="transformed_extent.yaml")
            
            data_c, weight_c, extent_c = transformer.transform_cake(data, self.exposure_time, weight)
            self.save_tiff(data_c, "caked_data.tif")
            self.save_tiff(weight_c, "transformed_image_weights.tif")
            self.save_yaml(data={"azimuthal": (extent_c[0], extent_c[1], data_c.shape[1]),
                                 "q_magnitude": (extent_c[2], extent_c[3], data_c.shape[0])},
                           filename="caked_extent.yaml")

        if self.plot:
            plt.show()

    def save_tiff(self, image, filename):
        fabio.tifimage.tifimage(image).write(self.directory / filename)

    def save_yaml(self, data, filename):
        with open(self.directory / filename, "w") as f:
            yaml.dump(data, f)
    
    def next(self):
        self.flag = next(self.iterator).lower()
        self.handle_flag()

    def handle_flag(self):
        if "=" in self.flag:
            flag, param = self.flag.split("=")
            if flag in ["--alpha", "--incident", "--in", "-i"]:
                self.incident_load = False
                self.incident_angle = float(param)
            elif flag in ["--tilt", "-t"]:
                self.tilt_load = False
                self.tilt_angle = float(param)
            elif flag in ["--expose", "--exposure", "--exposuretime", "--time", '-e']:
                self.exposure_load = False
                self.exposure_time = int(param)
            elif flag in ["--directory", "--dir", "-d"]:
                self.directory = Path(param)
        elif self.flag in ["--help", "-h"]:
            with open("help.txt", "r") as f:
                print(f.read())
    
    def parse_flags(self):
        while True:
            try:
                self.next()
            except StopIteration:
                break
    

def main(args:list=None):
    if args is None:
        args = sys.argv[1:]
    handler = Handler("gix", args)
    
    
def stitch(args:list=None):
    if args is None:
        args = sys.argv[1:]
    handler = Handler("stitch", args)