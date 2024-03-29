from datetime import datetime
from pathlib import Path


def create_file(directory: Path, angles: list[float], tag: str = "") -> None:
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
        f.write("umvr wbs -5\n")  # move beam stop out of the way
        f.write("umvr z -10\n")  # move sample out of the way
        f.write("eiger_run 0.1 om_scan_direct_beam.tif\n")  # take direct beam exposure
        f.write("umvr z 10\n")  # move sample back into beam

        for om in angles:
            f.write(f"umv om {om}\n")
            formatted_angle = "{}_{}".format(*str(om).split("."))
            f.write(f"eiger_run 0.1 om_scan_{tag}{formatted_angle}_degrees.tif\n")

        f.write("umvr wbs 5\n")
        f.write("umv om 0\n")
    print("Macro written")
    print("Copy and paste the following into SAXS to run the macro:")
    print("do " + str(directory / macroname))
    return None


def arange_angles(start, finish, step):
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