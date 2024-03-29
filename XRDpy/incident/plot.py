import os
import tkinter.filedialog
import numpy as np
import fabio
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib
from XRDpy.tiff_loader import Detector
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
from matplotlib.animation import FuncAnimation


__color_scheme_choices = ["plasma", "hot", "PuBu_r", "afmhot", "gist_heat", "copper", "viridis", "inferno", "magma", "cividis",
                          "gray", "pink", "winter"]


def load_tiff_data(directory: Path) -> tuple:
    """
    Load all relevant data from a directory
    :param files:
    :return: angles, intensity_data, intensity_db
    """
    base_mask = np.ones((Detector.ROWS, Detector.COLS), dtype=np.uint32)
    base_mask[Detector.DEAD_PIXELS] = 0
    
    for row_start_quad in Detector.ISSUE_QUADS[0]:
        for col_start_quad in Detector.ISSUE_QUADS[1]:
            base_mask[row_start_quad:(row_start_quad + 4), col_start_quad:(col_start_quad + 4)] = 0

    direct_beam_file = directory / "om_scan_direct_beam.tif"
    intensity_db = __load_image(direct_beam_file, base_mask)
    
    tiff_glob = directory.glob("om_scan*degrees.tif")

    num_angles = len(list(tiff_glob))

    angles = np.empty(num_angles)
    intensity_data = np.empty((num_angles, *intensity_db.shape))

    for ii, tiff in enumerate(tiff_glob):
        angles[ii] = float(__extract_angle_from(tiff))
        intensity_data[ii] = __load_image(tiff, base_mask)
    
    sorting_args = np.argsort(angles)
    angles = angles[sorting_args]
    intensity_data = intensity_data[sorting_args]

    return angles, intensity_data, intensity_db


def crop_data(exposure_datas: np.ndarray, direct_beam_data: np.ndarray, width: int, above: int = 50, below: int = 20,
              x_offset: int = 0) -> tuple[np.ndarray]:
    """
    Crop TIFF data around where the direct beam appears.
    :param exposure_datas: 3D numpy array of counts of all the exposures (exposure, z, x).
    :param direct_beam_data: 2D numpy array of counts of the direct beam exposure (z, x).
    :param width: Number of horizontal pixels to keep centered around the direct beam center. Can estimate with
                  (beam horizontal shutter size [s2hg]) / (pixel width): e.g. (0.8 mm) / (0.075 mm/px) ~ 11 px.
    :param above: The number of pixels above the direct beam to keep (~50 recommended).
    :param below: The number of pixels below the direct beam to keep.
    :param x_offset: optional offset for the width center if you find that the beam center isn't being properly found.
    :return: cropped exposure data (3D array), cropped direct beam data (2D array).
    """

    db_loc_x, db_loc_z = __find_direct_beam_pixels(direct_beam_data)
    z_lo = db_loc_z[0] - above
    z_hi = db_loc_z[1] + below
    x_ave = (db_loc_x[0] + db_loc_x[1]) / 2

    x_lo = int(x_ave - width * 0.5)
    x_hi = round(x_ave + width * 0.5)
    x_lo += x_offset
    x_hi += x_offset
    exposure_datas = exposure_datas[:, z_lo:z_hi, x_lo:x_hi]
    direct_beam_data = direct_beam_data[z_lo:z_hi, x_lo:x_hi]
    # print(exposure_datas.shape)
    return exposure_datas, direct_beam_data


"""Plotting Functions"""


def show() -> None:
    """
    Show all plots.
    :return: None
    """
    plt.show()


def display_tiff(tiff_data: np.ndarray, clip: int | float | None = None, log: bool = False) -> None:
    """
    Display a TIFF image.
    :param tiff_data: 2D numpy array containing pixel counts data (z, x)
    :param clip: an optional counts scale data clipping.
    :param log: optionally set log scale to counts.
    :return: None
    """
    if clip is not None:
        tiff_data[tiff_data > clip] = clip
    if log:
        tiff_data[tiff_data == 0] = 1
        tiff_data = np.log10(tiff_data)
    plt.imshow(tiff_data)
    return None


def animate_tiffs(intensity: np.ndarray, fps: float = 24.0, clip: int | float | None = None, log: bool = False,
                  color_scheme: int = 0):
    """
    Animate TIFFs where each frame of the animation is a different exposure.
    :param intensity: 3D numpy array (exposure, z, x).
    :param fps: frames per second.
    :param clip: an optional counts scale data clipping.
    :param log: True will result in a log scale display of intensity.
    :param color_scheme: choice of color scheme (given with an int value).
    :return: figure, animation (need to return these to successfully display)
    """
    fig = plt.figure(figsize=(1.2 + intensity.shape[2] * 0.2, 6))
    fig.tight_layout()

    frame_delay = 1000. / fps
    
    color_scheme = color_scheme % len(__color_scheme_choices)

    extension = None
    if clip is not None:
        intensity[intensity > clip] = clip
        extension = "max"
    if log:
        intensity[intensity < 1] = 1
        color_norm = colors.LogNorm(vmin=intensity.min(), vmax=intensity.max())
    else:
        color_norm = colors.Normalize(vmin=intensity.min(), vmax=intensity.max())

    ax = plt.axes(xlim=(0, intensity[0].shape[1]),
                  ylim=(0, intensity[0].shape[0]))
    im = plt.imshow(intensity[0][::-1], interpolation='none', norm=color_norm, cmap=__color_scheme_choices[color_scheme])
    cbar = fig.colorbar(im, ax=ax, extend=extension)

    def _init():
        im.set_data(intensity[0][::-1])
        return im,

    def _update(ii):
        im.set_array(intensity[ii][::-1])
        return im,

    ani = FuncAnimation(fig, _update, frames=intensity.shape[0],
                        init_func=_init, interval=frame_delay, blit=True)

    return fig, ani


def plot_tuning(angle: np.ndarray, intensity: np.ndarray, pixel_size: float | None = None,
                clip: float | int | None = None, log: bool = False, color_scheme: int = 0) -> tuple:
    """
    Plot angle tuning data.
    :param angle: motor angle at each exposure.
    :param intensity: intensity data (exposure, z, x).
    :param pixel_size: optionally scale data to mm with pixel size in mm.
    :param clip: optionally clip intensity data to a given value.
    :param log: True will result in a log scale display of intensity.
    :param color_scheme: choice of color scheme (given with an int value).
    :return: fig, ax, color_map (pcolor), color_bar
    """
    color_scheme = color_scheme % len(__color_scheme_choices)

    intensity_x_ave = np.sum(intensity, axis=2).astype(float)

    extension = None

    if clip is not None:
        intensity_x_ave[intensity_x_ave > clip] = clip
        extension = "max"
    if log:
        intensity_x_ave[intensity_x_ave < 1] = 1
        color_norm = colors.LogNorm(vmin=intensity_x_ave.min(), vmax=intensity_x_ave.max())
    else:
        color_norm = colors.Normalize(vmin=intensity_x_ave.min(), vmax=intensity_x_ave.max())

    fig, ax = plt.subplots()

    if pixel_size is None:
        z = np.arange(intensity_x_ave.shape[1])
        ax.set_ylabel("Z (pixels)")
    else:
        z = (np.arange(intensity_x_ave.shape[1])[::-1] + 0.5) * pixel_size
        ax.set_ylabel("Z (mm)")

    color_map = ax.pcolormesh(angle, z, intensity_x_ave.T, norm=color_norm, cmap=__color_scheme_choices[color_scheme])

    color_bar = fig.colorbar(color_map, ax=ax, extend=extension)

    ax.set_xlabel("angle (degrees)")
    ax.set_label("Intensity (counts)")

    return fig, ax, color_map, color_bar


"""Internal Functions for loading"""


def __extract_angle_from(file: Path) -> float:
    """
    Get the value of the angle from the name of a file.
    :param filename: The file in question.
    :return: The angle (in degrees) the data was taken at.
    """
    file_name_split = file.name.split("_")
    left_of_decimal = file_name_split[-3]
    angle = float(left_of_decimal)
    right_of_decimal = file_name_split[-2].replace(".tif", "")
    if left_of_decimal[0] == "-":
        angle -= float(right_of_decimal) / 10. ** len(right_of_decimal)
    else:
        angle += float(right_of_decimal) / 10. ** len(right_of_decimal)
    angle = round(angle, 3)
    return angle


def __find_direct_beam_file_index(files: list[str] | tuple[str]) -> int:
    """
    Get the index of the file which contains the direct beam data.
    :param files: list of files.
    :return: index of the direct beam file.
    """
    for ii, tif in enumerate(files):
        if tif.split("/")[-1].replace(".tif", "").lower() in ("db", "direct_beam"):
            return ii
    for ii, tif in enumerate(files):
        if tif.split("\\")[-1].replace(".tif", "").lower() in ("db", "direct_beam"):
            return ii


def __find_direct_beam_pixels(direct_beam_tiff: np.ndarray, threshold: int = 8) -> tuple[tuple[int]]:
    """
    Get list of indices where the direct beam is.
    :param direct_beam_tiff: full tiff
    :param threshold: the count threshold to look for the beam
    :return: list of direct beam locations.
    """
    direct_beam_bool = direct_beam_tiff > threshold     # bool of pixel brightness above a threshold.
    direct_beam_indx_z, direct_beam_indx_x = np.where(direct_beam_bool)  # indices of those bright pixels in z direction

    distances_z = np.diff(direct_beam_indx_z)               # get distances of a bright pixel from the next
    distances_x = np.diff(direct_beam_indx_x)  # get distances of a bright pixel from the next
    longest_repeat_x = direct_beam_indx_x[__get_indices_longest_repeated_element(distances_x)]
    longest_repeat_z = direct_beam_indx_z[__get_indices_longest_repeated_element(distances_z)]

    x_start = longest_repeat_x[0]
    x_end = longest_repeat_x[-1]
    z_start = longest_repeat_z[0]
    z_end = longest_repeat_z[-1]

    return (x_start, x_end), (z_start, z_end)


def __get_indices_longest_repeated_element(array: np.ndarray) -> np.ndarray:
    """
    Find the indices of the longest repeated element in an array.
    :param array: Array to sift through.
    :return: List of indices of locations of the longest repeated element
    """
    return __find_sequence(array, __get_longest_repeated_element(array))


def __get_longest_repeated_element(arr: np.ndarray) -> np.ndarray:
    """
    Return the longest repeated sequence of an element in an array.
    :param arr: array to look through.
    :return: longest sequence of repeated elements
    """
    mask = np.concatenate(([False], np.diff(arr) == 0, [False]))
    idx = np.flatnonzero(mask[1:] != mask[:-1])
    islands = [arr[idx[i]:idx[i+1] + 1] for i in range(0, len(idx), 2)]
    island_lengths = np.array([len(island) for island in islands])
    return islands[np.argmax(island_lengths)]


def __find_sequence(arr: np.ndarray, seq: np.ndarray) -> np.ndarray:
    """
    Search array for a specific sequence
    :param arr: array to search
    :param seq: array to find locations of
    :return: array of indices of where the sequence lies in the array
    """
    n_arr = arr.size
    n_seq = seq.size
    r_seq = np.arange(n_seq)
    match = (arr[np.arange(n_arr - n_seq + 1)[:, None] + r_seq] == seq).all(1)

    # Get the range of those indices as final output
    if match.any() > 0:
        return np.where(np.convolve(match, np.ones(n_seq, dtype=int)) > 0)[0]
    else:
        return np.array([])  # No match found


def __load_image(filename: Path, mask: np.ndarray) -> np.ndarray:
    """
    Load a TIF file as data.
    :param filename: Full path to TIF file.
    :return: Pixel intensity in a 2D numpy array
    """
    image = fabio.open(filename).data
    mask[np.where(image == Detector.MAX_INT)] = 0
    return image * mask
