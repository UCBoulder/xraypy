from datetime import datetime
from pathlib import Path


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
        f.write("umvr wbs 5\n")  # move beam stop out of the way

        for om in angles:
            f.write(f"umv om {om}\n")
            formatted_angle = "{}_{}".format(*str(om).split("."))
            f.write(f"eiger_run 0.1 om_scan_{tag}{formatted_angle}_degrees.tif\n")

        f.write("umvr z -10\n")  # move sample out of the way
        f.write("eiger_run 0.1 om_scan_direct_beam.tif\n")  # take direct beam exposure
        f.write("umvr z 10\n")  # move sample back into beam
        f.write("umvr wbs -5\n")
        f.write(f"umv om {final}\n")
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
        f.write("umvr wbs 5\n")  # move beam stop out of the way

        for z in zs:
            f.write(f"umv z {z}\n")
            formatted_angle = "{}_{}".format(*str(z).split("."))
            f.write(f"eiger_run 0.1 z_scan_{tag}{formatted_angle}_mm.tif\n")
        f.write(f"umv z {final}\n")
        f.write("eiger_run 0.1 z_scan_direct_beam.tif\n")  # take direct beam exposure
        f.write("umvr wbs -5\n")
    num = len(zs) + 1
    time_min = float(num) * 0.1
    minutes = int(time_min)
    seconds = round((time_min - minutes) * 60)
    print(f"Macro written with {num} images. Estimated time (min:sec): {minutes}:{seconds:02d}")
    print("Copy and paste the following into SAXS to run the macro:")
    print("do " + (directory / macroname).as_posix())
    print(f"WARNING: will leave z at {final} degrees")
    return None

def create_rel_file(directory: Path, motor: str, below: float, above: float, step: float, tag: str = "") -> None:
    if motor.lower() not in ["om", "z"]:
        raise ValueError("Motor type must be 'om' or 'z'.")
    below = -abs(below)
    above = abs(above)
    if not directory.exists():
        directory.mkdir(parents=True, exist_ok=True)
    date = datetime.now()
    macroname = f'Specular_{motor}-relative_macro-{date.year:02d}{date.month:02d}{date.day:02d}-{date.hour:02d}.txt'
    print("Writing Macro...")
    n = 0
    with open(directory / macroname, 'w') as f:
        f.write("umvr wbs 5\n")  # move beam stop out of the way
        f.write(f"umvr {motor} {below:.5f}")
        where = below
        while where < above:
            formatted_position = "{}_{}".format(*str(where).split("."))
            f.write(f"eiger_run 0.1 {motor}r_scan_{tag}{formatted_position}_mm.tif\n")
            where += step
            f.write(f"umvr {motor} {step:.5f}")
            n += 1
        to_move_back = below - n * step 
        f.write(f"umvr {motor} {to_move_back - 5}")
        # f.write("eiger_run 0.1 om_scan_direct_beam.tif\n")  # take direct beam exposure
        f.write(f"umvr {motor} 5")
        f.write("umvr wbs -5\n")
    num = below + above
    time_min = float(num) * 0.1
    minutes = int(time_min)
    seconds = round((time_min - minutes) * 60)
    print(f"Macro written with {num} images. Estimated time (min:sec): {minutes}:{seconds:02d}")
    print("Copy and paste the following into SAXS to run the macro:")
    print("do " + (directory / macroname).as_posix())
    


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