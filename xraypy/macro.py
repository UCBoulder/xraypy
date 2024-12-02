from datetime import datetime
from pathlib import Path
import tomllib
import xraypy.package_params as package

with open(package.config, "rb") as file:
    config = tomllib.load(file)["macro"]

set_om = config["set_omega"]
set_vert = config["set_vertical"]
set_hori = config["set_horizontal"]
move_bs = config["move_beamstop"]
move_om = config["move_omega"]
move_vert = config["move_vertical"]
move_hori = config["move_horizontal"]
expose = config["expose"]


def create_om_file(directory: Path, angles: list, tag: str = "", final: int = -4) -> None:
    """
    Create a macro file for scanning angle for GIWAXS
    :param angles: list of angles to scan through
    :param tag: optional identifier to put in filename
    :return: None
    """
    if not directory.exists():
        directory.mkdir(parents=True, exist_ok=True)
    print(angles)
    if tag:
        tag += "_"
    date = datetime.now()
    macroname = f'Incident_angle_tuning_macro-{date.year:02d}{date.month:02d}{date.day:02d}-{date.hour:02d}.txt'
    print("Writing Macro...")
    with open(directory / macroname, 'w') as f:
        f.write(move_bs.format(5))  # move beam stop out of the way

        for om in angles:
            f.write(set_om.format(om))
            formatted_angle = "{}_{}".format(*str(om).split("."))
            img_name = f"om_scan_{tag}{formatted_angle}_degrees"
            f.write(expose.format(img_name))

        f.write(move_vert.format(-10))  # move sample out of the way
        f.write(expose.format("om_scan_direct_beam"))  # take direct beam exposure
        f.write(move_vert.format(10))  # move sample back into beam
        f.write(move_bs.format(-5))
        f.write(set_om.format(final))
    num = len(angles) + 1
    time_min = float(num) * 0.1
    minutes = int(time_min)
    seconds = round((time_min - minutes) * 60)
    print(f"Macro written with {num} images. Estimated time (min:sec): {minutes}:{seconds:02d}")
    print("Copy and paste the following into SAXS to run the macro:")
    print("do " + (directory / macroname).as_posix())
    print(f"WARNING: will leave om at {final} degrees")
    return None


def create_z_file(directory: Path, zs: list, tag: str = "", final: int = -5) -> None:
    """
    Create a macro file for scanning angle for GIWAXS
    :param zs: list of z-positions to scan through
    :param tag: optional identifier to put in filename
    :return: None
    """
    if not directory.exists():
        directory.mkdir(parents=True, exist_ok=True)
    print(zs)
    if tag:
        tag += "_"
    date = datetime.now()
    macroname = f'Specular_z_macro-{date.year:02d}{date.month:02d}{date.day:02d}-{date.hour:02d}.txt'
    print("Writing Macro...")
    with open(directory / macroname, 'w') as f:
        f.write(move_bs.format(5))  # move beam stop out of the way

        for z in zs:
            f.write(set_vert.format(z))
            formatted_angle = "{}_{}".format(*str(z).split("."))
            img_name = f"z_scan_{tag}{formatted_angle}_mm"
            f.write(expose.format(img_name))
        f.write(set_vert.format(final))
        f.write(expose.format("z_scan_direct_beam"))  # take direct beam exposure
        f.write(move_bs.format(-5))
    num = len(zs) + 1
    time_min = float(num) * 0.1
    minutes = int(time_min)
    seconds = round((time_min - minutes) * 60)
    print(f"Macro written with {num} images. Estimated time (min:sec): {minutes}:{seconds:02d}")
    print("Copy and paste the following into SAXS to run the macro:")
    print("do " + (directory / macroname).as_posix())
    print(f"WARNING: will leave z at {final} mm")
    return None  


def arange_list(start, finish, step):
    """
    Make a list of values similar to np.arange, but with values rounded to avoid floating point precision issues
    :param start: first element of the list
    :param finish: last element of the list
    :param step: difference between sequential elements
    :return: list of values
    """
    start = float(start)
    finish = float(finish)
    step = float(step)
    # Try to determine what digit to round to
    step_decimal = str(step).split(".")  # list[str]: [left of decimal, right of decimal]
    if step_decimal[1] == 0:        # then step is an integer
        rounding = 0
        # find lowest order non-zero digit
        while True:
            if step_decimal[0][::-1].find('0') == 0:
                step_decimal[0] = step_decimal[0][:-1]
                rounding -= 1
            else:
                break
    else:                           # then step is not an integer
        rounding = len(step_decimal[1])     # number of digits right of the decimal
    return [round(x * step + start, rounding) for x in list(range(int((finish + step - start) / step)))]