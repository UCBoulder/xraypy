import argparse
import shutil
import yaml
import matplotlib.pylab as plt
from pathlib import Path
from datetime import datetime
import XRDpy.package_params as package
import XRDpy.incident.macro as macro
import XRDpy.incident.specular as specular


def plot():
    parser = argparse.ArgumentParser(
        prog="plot",
        description="Plot om scan as an image with horizontal integration over pixels, rows: y-pixels, columns: om, color: intensity\nThis will look for data in the latest directory of the user",
        epilog="author: Teddy Tortorici <edward.tortorici@colorado.edu>"
    )
    parser.add_argument("user", help="This will be the name of the directory the files will be looked for")
    parser.add_argument("-A", "--animate", help="animate: specify frame-rate in FPS")
    parser.add_argument("-T", "--title")
    parser.add_argument("-M", "--mod", help="modify fit or plot (z0 adjust for om or change std for z)")
    parser.add_argument("-C", "--crit", help="add critical angles (comma separated)")
    parser.add_argument("-Z", "--z", action="store_true", help="change to a z-specular scan")
    parser.add_argument("-D", "--dir", help="Specify a specific directory, or CWD for current working directory")
    parser.add_argument("-R", "--range", help="set angular range in degrees")
    parser.add_argument("-B", "--beamwidth", help="set beam width in mm")
    parser.add_argument("-S", "--save", help="save the plot at a certain DPI")
    args = parser.parse_args()

    if args.z:
        spec_type = "z"
    else:
        spec_type = "om"
    usr_dir = package.directory / args.user / "{}_scans".format(spec_type)

    if args.date is None:
        date = datetime.now()
        date_name = f"{date.year}{date.month:02d}{date.day:02d}"
    else:
        date_name = args.date
        if len(date_name) != 8:
            raise ValueError(f"The given date must be in the form YYYYMMDD.")
        try:
            int(date_name)
        except ValueError:
            raise ValueError(f"The given date must be in the form YYYYMMDD with all integer values.")
    
    if args.number is None:
        directory = Path(max(usr_dir.glob(f"{date_name}-*"), key=lambda d: int(d.name.split("-")[-1]), default=1).rstrip("-1"))
        date_name = f"{date_name}-{dir_num}"
        directory = usr_dir / date_name.rstrip("-1")
    else:
        dir_num = int(args.number)
        date_name = f"{date_name}-{dir_num}"
        directory = usr_dir / date_name.rstrip("-1")

    if args.dir is not None:
        if args.dir.upper() == "CWD":
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

    if args.z:
        if args.mod is None:
            std = float(args.mod)
        else:
            std = 4.
        spec = specular.SpecularZ(directory, angular_range=angular_range, beam_width=beamwidth, standard_deviations=std)
    else:
        spec = specular.SpecularOmega(directory, anglular_range=angular_range, beam_width=beamwidth)
        if args.crit is None:
            crit = None
        else:
            crit = [float(c) for c in args.mod.split(",")]
        if args.mod is not None:
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
    parser.add_argument("-A", "--append", action="store_true", help="Add this tag to append the latest set of data")
    args = parser.parse_args()

    if args.type.lower() not in ["om", "z", "omr", "zr"]:
        raise ValueError("Type must be 'om', 'z', 'omr', or 'zr'.")
    
    with open(package.directory / "config.yaml") as f:
        data_path = Path(yaml.safe_load(f)["data_path"])
    
    usr_dir = package.directory / args.user / "{}_scans".format(args.type.lower())
    if not usr_dir.exists():
        usr_dir.mkdir(parents=True, exist_ok=True)
    
    date_name = f"{date.year}{date.month:02d}{date.day:02d}"
    
    directory = usr_dir / date_name
    
    if directory.exists():
        print("Directory has already been created today")
        current_number = 1
        for sister_dir in usr_dir.glob(f"{date_name}-*"):
            sister_num = int(sister_dir.name.split("-")[-1])
            if sister_num > current_number:
                current_number = sister_num
        current_dir = usr_dir / f"{date_name}-{current_number}".rstrip("-1")
        print(f"Latest directory found: {current_dir.as_posix()}")
        if not args.append:
            directory = usr_dir / f"{date_name}-{current_number + 1}"
            directory.mkdir(parents=True, exist_ok=True)
            print(f"Making new directory: {directory.as_posix()}")
        else:
            directory = current_dir
            print("Appending to this directory")
    else:
        directory.mkdir(parents=True, exist_ok=True)
        print(f"Making new directory: {directory.as_posix()}")

    scan_glob = data_path.glob(f"{args.type.lower()}_scan*.tif")
    print(f"Found {len(list(scan_glob))} files to move")

    for tif in data_path.glob(f"{args.type.lower()}_scan*.tif"):
        new_tif = directory / tif.name
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
    if args.type.lower not in ["om", "z"]:
        raise ValueError("Invalid type. Must be 'om' or 'z'.")
    if args.clear:
        for macro in dir.glob("*.txt"):
            macro.unlink()
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




