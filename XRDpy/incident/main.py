import argparse
import shutil
import yaml
import pickle
import matplotlib.pylab as plt
from pathlib import Path
from datetime import datetime
import XRDpy.package_params as package
import XRDpy.incident.macro as imacro
import XRDpy.incident.plot as iplot


def plot():
    parser = argparse.ArgumentParser(
        prog="XRD-iplot",
        description="Plot om scan as an image with horizontal integration over pixels, rows: y-pixels, columns: om, color: intensity\nThis will look for data in the latest directory of the user",
        epilog="author: Teddy Tortorici <edward.tortorici@colorado.edu>"
    )
    parser.add_argument("user", help="This will be the name of the directory the files will be looked for")
    parser.add_argument("-A", "--animate", action="store_true")
    parser.add_argument("-T", "--title")
    parser.add_argument("-D", "--date", help="Override current date with format YYYYMMDD")
    parser.add_argument("-N", "--number", help="Override latest directory of the given (or current) date")
    parser.add_argument("-CW", "--crop_width", help="Override number of pixels in horizontal crop (default=13)")
    parser.add_argument("-CA", "--crop_above", help="Override number of pixels above the direct beam to keep in crop (default=100).")
    parser.add_argument("-CB", "--crop_below", help="Override number of pixels below the direct beam to keep in crop (default=20).")
    parser.add_argument("-CO", "--crop_offset", help="Move the center of the crop to the right (in pixels)")
    parser.add_argument("-C", "--color", help="Integer value corresponding to different color schemes")
    args = parser.parse_args()
    
    usr_dir = package.directory / args.user

    if args.date is None:
        date = datetime.now()
        dir_date = f"{date.year}{date.month}{date.day}"
    else:
        dir_date = args.date
        if len(dir_date) != 8:
            raise ValueError(f"The given date must be in the form YYYYMMDD.")
        try:
            int(dir_date)
        except ValueError:
            raise ValueError(f"The given date must be in the form YYYYMMDD with all integer values.")
    
    if args.number is None:
        sister_dirs = usr_dir.glob(f"{dir_date}")
        trailing_numbers = [int(d.name.split("-")[-1]) for d in sister_dirs]
        if trailing_numbers:
            dir_name = f"{dir_date}-{max(trailing_numbers)}"
        else:
            dir_name = dir_date
    else:
        dir_name = f"{dir_date}-{int(args.number)}"
        dir_name.rstrip("-1")
    
    directory = usr_dir / dir_name

    angles, intensity_data, direct_beam = iplot.load_tiff_data(directory)

    print(f"Loaded {len(angles) + 1} files from {directory.as_posix()}")
    print(f"Found {len(angles)} angles in the scan")
    print(f"Data shape is: {intensity_data.shape}")
    print(f"Direct beam data shape is: {direct_beam.shape}")

    x_pixel_width = 13
    pixels_above = 100
    pixels_below = 20
    x_offset = 0
    if args.crop_width is not None:
        x_pixel_width = int(args.crop_width)
    if args.crop_above is not None:
        pixels_above = int(args.crop_above)
    if args.crop_below is not None:
        pixels_below = int(args.crop_below)
    if args.crop_offset is not None:
        x_offset = int(args.crop_offset)
    
    id_c, db_c = iplot.crop_data(intensity_data, direct_beam, x_pixel_width, pixels_above, pixels_below, x_offset)

    color_code = 0
    if args.color is not None:
        color_code = int(args.color) % len(iplot.color_scheme_choices)
    
    if args.animate:
        fps = 48
        fig, ani = iplot.animate_tiffs(id_c, fps, log=True, color_scheme=color_code)
    
    iplot.plot_tuning(angles, id_c, pixel_size=0.75, log=True, color_scheme=color_code)
    if args.title is not None:
        plt.title(args.title)

    data_to_save = {"angles": angles,
                    "data": id_c}
    with open(directory / "data.pkl", "wb") as f:
        pickle.dump(data_to_save, f)
    print("Saved data in: " + (directory / "data.pkl").as_posix())
    
    plt.show()
    

def move():
    date = datetime.now()
    parser = argparse.ArgumentParser(
        prog="XRD-imove",
        description="Move files from Data to a user folder.",
        epilog="author: Teddy Tortorici <edward.tortorici@colorado.edu>"
    )
    parser.add_argument("user", help="This will be the name of the directory the files will be moved to")
    parser.add_argument("-A", "--append", help="Add this tag to append the latest set of data")
    args = parser.parse_args()
    with open(package.directory / "config.yaml") as f:
        data_path = Path(yaml.safe_load(f)["data_path"])

    usr_dir = package.directory / args.user
    if not usr_dir.exists():
        usr_dir.mkdir(parents=True, exist_ok=True)
    
    date_name = f"{date.year}{date.month}{date.day}"
    
    directory = usr_dir / date_name

    if directory.exists():
        if not args.append:
            sister_dirs = usr_dir.glob(f"{date_name}*")
            trailing_numbers = [int(d.name.split("-")[-1]) for d in sister_dirs]
            if trailing_numbers:
                next_number = max(trailing_numbers) + 1
            else:
                next_number = 2
            directory = usr_dir / f"{date_name}-{next_number}"
            print(f"Making new directory: {directory.as_posix()}")
    else:
        directory.mkdir(parents=True, exist_ok=True)
        print(f"Making new directory: {directory.as_posix()}")

    scan_glob = data_path.glob("om_scan*.tif")
    print(f"Found {len(list(scan_glob))} files to move")

    for tif in scan_glob:
        new_tif = directory / tif.name
        if not new_tif.is_file():
            new_tif.unlink()
        shutil.move(tif, new_tif)
    scan_glob = directory.glob("om_scan*.tif")
    print(f"Moved {len(list(scan_glob))} files")


def make():
    parser = argparse.ArgumentParser(
        prog="XRD-imacro",
        description="Create a macro for SAXS to scan through om motor angles to calibrate incident angle to a flat sample.",
        epilog="author: Teddy Tortorici <edward.tortorici@colorado.edu>"
    )
    parser.add_argument("start_angle")
    parser.add_argument("end_angle")
    parser.add_argument("angle_step")
    # parser.add_argument("-T", "-thickness", help="sample thickness")
    parser.add_argument("-N", "--name", help="add a short tag to the TIFF files generated")
    parser.add_argument("-C", "--clear", action='store_true', help="remove all macro files saved")
    dir = Path.home() / "XRDpy" / "Macros"
    args = parser.parse_args()
    if args.clear:
        for macro in dir.glob("*.txt"):
            macro.unlink()
    angles = imacro.arange_angles(args.start_angle, args.end_angle, args.angle_step)
    if args.name is None:
        args.name = ""
    imacro.create_file(dir, angles, args.name)
