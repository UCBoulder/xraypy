from datetime import datetime
import argparse
import yaml
import fabio
from pathlib import Path
import XRDpy.package_params as package
import shutil
from XRDpy.tiff_loader import Stitcher
import XRDpy.macro as macro
import XRDpy.specular as specular
import matplotlib.pylab as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


class StitchParse:

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            prog="stitch",
            description="Stitch together raw images from Eiger R 1M that eiger-stitch produced.",
            epilog="Author: Teddy Tortorici edward.tortorici@colorado.edu"
        )
        self.parser.add_argument("experiment", help="experiment type. Example: WAXS, SAXS, GIWAXS, GISAXS")
        self.parser.add_argument("rows")
        self.parser.add_argument("columns")
        self.parser.add_argument("exposure")
        self.parser.add_argument("-P", "--dir", help="where to find data")
        self.parser.add_argument("-I", "--incident", help="incident angle in degrees (for GIXS)")
        self.parser.add_argument("-O", "--om", help="omega motor position in degrees (for GIXS)")
        self.parser.add_argument("-Z", "--z", help="z motor position in mm (for GIXS)")
        self.parser.add_argument("-D", "--dist", help="sample detector distance in meters")
        self.parser.add_argument("-T", "--tif", action="store_true")
        self.parser.add_argument("-J", "--jupyter", action="store_true")

        self.args = self.parser.parse_args()
        
        if self.args.dir is None:
            self.dir = Path.cwd()
        elif self.args.dir == "DATA":
            with open(package.directory / "config.yaml") as f:
                self.dir = Path(yaml.safe_load(f)["data_path"])
        else:
            self.dir = Path(self.args.dir)

        newest_file = self.get_newest_file(self.dir)
        filename_base = newest_file.name.strip(".tif")

        stitcher = Stitcher(int(self.args.rows), int(self.args.columns))
        data, flat_field = stitcher.load_data(self.dir)

        date = datetime.now()
        header = {
            "StitchingRows": int(self.args.rows),
            "StitchingColumns": int(self.args.columns),
            "ExposureTime(s)": float(self.args.exposure),
            "ExperimentType": self.args.experiment,
            "IncidentAngle(deg)": self.args.incident,
            "OmegaMotorAngle(deg)": self.args.om,
            "ZMotorPosition(mm)": self.args.z,
            "DetectorDistance(m)": self.args.dist,
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

        if self.args.jupyter:
            if "GI" in self.args.experiment.upper() or self.args.incident is not None or self.args.om is not None:
                notebooks_to_copy = ("transform-GIWAXS.ipynb")
            else:
                notebooks_to_copy = ("reduce_WAXS.ipynb",)
            for notebook in notebooks_to_copy:
                self.copy_notebook(notebook)
        

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
    
    def copy_files(self):
        if self.args.file is not None:
            for file in self.args.file:
                self.copy_file(file)

    @staticmethod
    def get_newest_file(directory: Path) -> Path:
        files = list(directory.glob("*.tif"))
        return max(files, key=lambda f: f.stat().st_mtime)
        

def stitch():
    parser = StitchParse()


def plot():
    parser = argparse.ArgumentParser(
        prog="plot",
        description="Plot om scan as an image with horizontal integration over pixels, rows: y-pixels, columns: om, color: intensity\nThis will look for data in the latest directory of the user",
        epilog="author: Teddy Tortorici <edward.tortorici@colorado.edu>"
    )
    parser.add_argument("dir", help="Specify a specific directory, or CWD for current working directory")
    parser.add_argument("-A", "--animate", help="animate: specify frame-rate in FPS")
    parser.add_argument("-T", "--title")
    parser.add_argument("-P", "--zposition", help="move the z-position of the sample relative to the beam center")
    parser.add_argument("-C", "--crit", help="add critical angles (comma separated)")
    parser.add_argument("-B", "--beamheight", help="change where the beam cutoff is, in standard deviations.")
    parser.add_argument("-R", "--range", help="set angular range in degrees")
    parser.add_argument("-W", "--beamwidth", help="set beam width in mm")
    parser.add_argument("-S", "--save", help="save the plot at a certain DPI")
    parser.add_argument("-N", "--name", help="Title plots")
    args = parser.parse_args()

    if args.dir.lower() == "cwd" or args.dir.lower() == "here":
        directory = Path.cwd()
    else:
        directory = Path(args.dir)
    
    if args.range is None:
        angular_range = 1.5
    else:
        angular_range = args.range
    if args.beamwidth is None:
        beamwidth = 1.
    else:
        beamwidth = args.beamwidth

    if args.beamheight:
        std = args.beamheight
    else:
        std = 3

    directory_name_split = directory.name.split("_")
    scan_type = directory_name_split[1]
    other_position = float(directory_name_split[2].split("-")[1])
    if scan_type == "om-scan":
        other_unit = "um"
    elif scan_type == "z-scan":
        other_unit = "millidegree"
    else:
        other_unit = ""
    
    plot_name = f"{scan_type}_at-{other_position}-{other_unit}"

    spec = specular.SpecularScan(
        directory,
        anglular_range=angular_range,
        beam_width=beamwidth,
        standard_deviations=std,
        plot_name=plot_name,
        plot_dpi=args.save,
        plot_title=args.name,
    )

    if spec.type == "om":
        if args.crit is None:
            crit = None
        else:
            crit = [float(c) for c in args.crit.split(",")]
        if args.zposition:
            spec.fit(z0=float(args.mod))
        spec.plot(critical_angle=crit)

    if args.animate is not None:
        fps = args.animate
        fig, ani = spec.animate_tiffs(fps)
    
    spec.save(directory)
    print("Saved data in: " + directory.as_posix())
    
    plt.show()
    

def move():
    date = datetime.now()
    parser = argparse.ArgumentParser(
        prog="move",
        description="Move files from Data to a user folder.",
        epilog="author: Teddy Tortorici <edward.tortorici@colorado.edu>"
    )
    parser.add_argument("type", help="Either om or z")
    parser.add_argument("user", help="This will be the name of the directory the files will be moved to")
    parser.add_argument("other_position", help="Specify the position of the other motor (the one note being scanned).")
    parser.add_argument("-A", "--append", action="store_true", help="Add this tag to append the latest set of data")
    args = parser.parse_args()

    if args.type.lower() not in ["om", "z", "omr", "zr"]:
        raise ValueError("Type must be 'om', 'z', 'omr', or 'zr'.")
    if "om" in args.type.lower():
        other_unit = "um"
    else:
        other_unit = "mdeg"
    
    with open(package.directory / "config.yaml") as f:
        DATA = Path(yaml.safe_load(f)["data_path"])
    
    usr_dir = package.directory / args.user / "scans"   # / "{}_scans".format(args.type.lower())
    if not usr_dir.exists():
        usr_dir.mkdir(parents=True, exist_ok=True)
    
    directory_name = f"{date.year}-{date.month:02d}-{date.day:02d}_{args.type.lower()}-scan_at-{args.other_position * 1000}-{other_unit}"
    
    data_path = usr_dir / directory_name
    
    if data_path.exists():
        print("Directory has already been created today")
        current_number = 1
        for sister_dir in usr_dir.glob(f"{directory_name}_*"):
            try:
                sister_num = int(sister_dir.name.split("_")[-1])
                if sister_num > current_number:
                    current_number = sister_num
            except ValueError:
                pass
        current_dir = usr_dir / f"{directory_name}_{current_number}".rstrip("_1")
        print(f"Latest directory found: {current_dir.as_posix()}")
        if not args.append:
            data_path = usr_dir / f"{directory_name}-{current_number + 1}"
            data_path.mkdir(parents=True, exist_ok=True)
            print(f"Making new directory: {data_path.as_posix()}")
        else:
            data_path = current_dir
            print("Appending to this directory")
    else:
        data_path.mkdir(parents=True, exist_ok=True)
        print(f"Making new directory: {data_path.as_posix()}")

    scan_glob = DATA.glob(f"{args.type.lower()}_scan*.tif")
    print(f"Found {len(list(scan_glob))} files to move")

    for tif in DATA.glob(f"{args.type.lower()}_scan*.tif"):
        new_tif = data_path / tif.name
        print(new_tif)
        if new_tif.is_file():
            new_tif.unlink()
        shutil.move(tif, new_tif)


def make_scan():
    parser = argparse.ArgumentParser(
        prog="macro_om",
        description="Create a macro for SAXS to scan through motor angles to calibrate for GIXS.",
        epilog="author: Teddy Tortorici <edward.tortorici@colorado.edu>"
    )
    parser.add_argument("type")
    parser.add_argument("start")
    parser.add_argument("end")
    parser.add_argument("step")
    parser.add_argument("-R", "--relative", action='store_true', help="relative scan from current position")
    # parser.add_argument("-T", "-thickness", help="sample thickness")
    parser.add_argument("-N", "--name", help="add a short tag to the TIFF files generated")
    parser.add_argument("-C", "--clear", action='store_true', help="remove all macro files saved")
    dir = Path.home() / "XRDpy" / "Macros"
    args = parser.parse_args()
    if args.type.lower() not in ["om", "z"]:
        raise ValueError("Invalid type. Must be 'om' or 'z'.")
    if args.clear:
        for macro_file in dir.glob("*.txt"):
            macro_file.unlink()
    if args.name is None:
        args.name = ""
    if args.relative:
        macro.create_rel_file(dir, args.type.lower(), args.start, args.end, args.step, args.name)
    else:
        list_to_scan = macro.arange_list(args.start, args.end, args.step)
        if args.type.lower() == "om":
            macro.create_om_file(dir, list_to_scan, args.name)
        elif args.type.lower() == "z":
            macro.create_z_file(dir, list_to_scan, args.name)

