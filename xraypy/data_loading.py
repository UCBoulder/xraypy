import numpy as np
from pathlib import Path
from scipy import ndimage
from XRDpy import package_params as package
import tomllib
import fabio
import pyFAI.detectors


def get_detector(detector_name: str) -> pyFAI.detectors._common.Detector:
    """
    Load the detector class from a string: i.e. "Eiger1M" -> pyFAI.detectors.Eiger1M
    """
    DetectorClass = getattr(pyFAI.detectors, detector_name, None)
    if DetectorClass is None:
        raise ValueError(f"Detector '{detector_name}' is not a valid detector type. Visit https://pyfai.readthedocs.io/en/stable/api/detectors/index.html#module-pyFAI.detectors")
    return DetectorClass


with open(package.config, "rb") as file:
    config = tomllib.load(file)

detector_config = config["detector"]
stitch_config = config["stitch"]

DetectorFromConfig = get_detector(detector_config["name"])

detector_shape = DetectorFromConfig().shape


class Detector(DetectorFromConfig):

    BAD_PIXELS = (np.array(detector_config['bad_pixels']['rows'], dtype=np.int64),   # rows
                  np.array(detector_config['bad_pixels']['columns'], dtype=np.int64))   # columns
    
    ISSUE_QUADS = (tuple(detector_config['issue_quads']['start_rows']),
                   tuple(detector_config['issue_quads']['start_columns']))
    
    ROWS = detector_shape[0]
    COLS = detector_shape[1]

    MAX_INT = detector_config["max_int"]

    def __init__(self):
        super().__init__(orientation=detector_config["orientation"])
    
    @classmethod
    def calc_mask(cls) -> np.ndarray:
        mask = super().calc_mask()
        mask[cls.BAD_PIXELS] = 1
        for row_start_quad in cls.ISSUE_QUADS[0]:
            for col_start_quad in cls.ISSUE_QUADS[1]:
                mask[row_start_quad:(row_start_quad + 4), col_start_quad:(col_start_quad + 4)] = 1
        return mask
    
    @classmethod
    def calc_mask_dezinger(cls, image: np.ndarray, cut_off: float=None,
                           gaussian_standard_deviation: float=None) -> np.ndarray:
        return np.logical_or(cls.find_zingers(image, cut_off, gaussian_standard_deviation), cls.calc_mask())
    
    @classmethod
    def get_size(cls) -> tuple:
        return (cls.ROWS, cls.COLS)
    
    @classmethod
    def get_rows(cls) -> int:
        return cls.ROWS
    
    @classmethod
    def get_columns(cls) -> int:
        return cls.COLS
    
    @staticmethod
    def find_zingers(image: np.ndarray, cut_off: float=5, gaussian_standard_deviation: float=2) -> np.ndarray:
        """
        Finds abnormally hot pixels by applying a gaussian filter to smooth the image and locate pixels and returns a mask for them
        :params image: image to find zingers in
        :params cut_off: 
        :params gaussian_standard_deviations:
        :return: return mask of zinger locations
        """
        if cut_off is None:
            cut_off = 5.
        if gaussian_standard_deviation is None:
            gaussian_standard_deviation = 2.
        smoothed_img = ndimage.gaussian_filter(image, gaussian_standard_deviation)
        dif_img = image - smoothed_img

        zinger_chart = dif_img / (smoothed_img + 1)
        anomalies = zinger_chart > cut_off
        print(f'Found {np.sum(anomalies)} zingers')
        return anomalies.astype(bool)
        

class Stitcher:

    OVERLAP = stitch_config["image_overlap"]
    MISSING_BAND_PIXELS = stitch_config["missing_band_rows"]
    WEIGHT_PER = stitch_config["exposure_per_image"]
    SEP = stitch_config["separator"]
    LABEL = stitch_config["raw_label"]
    EXT = stitch_config["file_extension"]
    
    def __init__(self, rows: int = 0, columns: int = 0, compensate: bool = stitch_config["compensate"]):
        """
        Define the stitch by the rows and columns and if the detector will be moved to compensate for the missing band.
        If rows=0 and columns=0 is given, it will look for all the raw images present in the directory to stitch.

        :params rows: number of stitching rows.
        :params columns: number of stitching columns.
        :params compensate: Will the detector be moved to compensate for a missing band of pixels?
        """
        self.detector = Detector()
        self.stitch_rows = rows
        self.stitch_columns = columns
        self.compensate = compensate
        if compensate:
            self.format = f"{self.LABEL}*{self.SEP}*{self.SEP}*{self.EXT}"
        else:
            self.format = f"{self.LABEL}*{self.SEP}*{self.EXT}"
        self.shape = (0, 0)
        
    def determine_shape(self):
        """
        Determine the shape of the stitched image.
        """
        rows = self.stitch_rows * (self.detector.get_rows() + Stitcher.MISSING_BAND_PIXELS - Stitcher.OVERLAP) + Stitcher.OVERLAP
        columns = self.stitch_columns * (self.detector.get_columns() - Stitcher.OVERLAP) + Stitcher.OVERLAP
        self.shape = (rows, columns)
        print(f"Stitched image will have shape: ({self.shape[0]}, {self.shape[1]}).")

    def stitch(self, directory: Path, dezinger: bool=False, cut_off: float=None, gaussian_standard_deviation: float=None):
        """
        Load data from directory based on self.rows and self.columns
        """
        print(f"Loading raw images from: {directory.as_posix()}\n")
        if not (self.stitch_rows and self.stitch_columns):  # if either of these are 0
            print("Will look for images to determine stitching size\n")
            # Find stitch size by using all files present
            for file in directory.glob(self.format):
                print(f"found: {file.name}")
                try:
                    name = file.stem.lstrip(self.LABEL)
                    stitch_row, stitch_col, _ = name.split(self.SEP)
                    stitch_row = int(stitch_row)
                    stitch_col = int(stitch_col)
                    print(f" - row: {stitch_row}; col: {stitch_col}")
                    if stitch_row > self.stitch_rows:
                        self.stitch_rows = stitch_row
                    if stitch_col > self.stitch_columns:
                        self.stitch_columns = stitch_col
                except ValueError:
                    print(" - not a valid raw image for stitching")
        stitch_description = f"Stitching {self.stitch_rows} rows and {self.stitch_columns} columns"
        if self.compensate:
            stitch_description += ", with compensation for empty band,"
        stitch_description += " for new image."
        print(stitch_description)
        self.determine_shape()

        data = np.zeros(self.shape, dtype=np.float64)
        flat_field = np.zeros(self.shape, dtype=np.float64)
        base_mask = np.logical_not(self.detector.calc_mask())

        comp_range = int(self.compensate) + 1
        
        for stitch_row in range(1, self.stitch_rows + 1):
            for stitch_col in range(1, self.stitch_columns + 1):
                for stitch_offset in range(1, comp_range + 1):
                    raw_image_file_name = f"{self.LABEL}{stitch_row}{self.SEP}{stitch_col}"
                    if self.compensate:
                        raw_image_file_name += f"{self.SEP}{stitch_offset}"
                    raw_image_file_name += self.EXT
                    file_data = fabio.open(raw_image_file_name).data
                    if dezinger:
                        mask = np.logical_not(self.detector.calc_mask_dezinger(file_data, cut_off, gaussian_standard_deviation))
                    else:
                        mask = base_mask.copy()                
                    mask = np.logical_not(self.detector.calc_mask)
                    mask[file_data == self.detector.MAX_INT] = 0
                    file_data *= mask
                    start_row = (self.stitch_rows - stitch_row) * (self.detector.ROWS - Stitcher.OVERLAP)
                    if self.compensate:
                        start_row += (2 - stitch_offset) * Stitcher.MISSING_BAND_PIXELS
                    start_column = (stitch_col - 1) * (self.detector.COLS - Stitcher.OVERLAP)
                    data[start_row:(start_row + self.detector.ROWS), start_column:(start_column + self.detector.COLS)] += file_data
                    flat_field[start_row:(start_row + self.detector.ROWS), start_column:(start_column + self.detector.COLS)] += mask
        flat_field *= self.WEIGHT_PER

        detector = pyFAI.detectors.Detector(pixel1=self.detector.pixel1,
                                            pixel2=self.detector.pixel2,
                                            max_shape=data.shape,
                                            orientation=int(self.detector.orientation))
        return data, flat_field, detector
