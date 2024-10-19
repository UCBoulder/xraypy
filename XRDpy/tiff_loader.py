import numpy as np
import yaml
from pathlib import Path
import fabio
from pyFAI.detectors import Eiger1M
import matplotlib.pylab as plt


def load_from(directory: Path, exposure_time: int = None):
    exposure_time = handle_yaml(directory, exposure_time)
    data_file_name = "raw-stitched-data.tif"
    weights_file_name = "flat-field.tif"
    if (directory / data_file_name).is_file() and (directory / weights_file_name).is_file():
        print("Found stitched files.")
        data = fabio.open(directory / data_file_name).data
        weight = fabio.open(directory / weights_file_name).data
    else:
        print("Stitching raw images.")
        stitcher = Stitcher()
        data, weight = stitcher.load_data(directory)
        data_im = fabio.tifimage.tifimage(data)
        weight_im = fabio.tifimage.tifimage(weight)
        data_im.write(directory / data_file_name)
        weight_im.write(directory / weights_file_name)
    return data, weight


def handle_yaml(directory: Path, exposure_time: int):
    yaml_file = directory / "params.yaml"
    if (directory / yaml_file).is_file():
        with open(yaml_file, 'r') as f:
            params = yaml.safe_load(f)
        if exposure_time is None:
            try:
                exposure_time = params["time"]
            except KeyError:
                raise KeyError("params.yaml does not contain 'time' key. You must give an exposure time (s) to create one.")
        else:
            try:
                params["time"] = exposure_time
            except TypeError:
                params = {"time": exposure_time}
            with open(yaml_file, 'w') as f:
                yaml.dump(params, f)
    else:
        if exposure_time is None:
            raise FileNotFoundError("params.yaml does not exist. You must give an exposure time (s) to create one.")
        else:
            with open(yaml_file, 'w') as f:
                yaml.dump({"time": exposure_time}, f)
    return exposure_time


class Detector(Eiger1M):

    DEAD_PIXELS = (np.array([681, 904, 398], dtype=np.int64),   # rows
                   np.array([646, 561, 579], dtype=np.int64))   # columns

    ISSUE_QUADS = ((255, 808),       # start rows
                   (255, 513, 771))  # start columns

    ROWS = 1065
    COLS = 1030

    MAX_INT = (1 << 32) - 1

    def __init__(self):
        super().__init__()
    
    def calc_mask(self) -> np.ndarray:
        mask = super().calc_mask()
        mask[self.DEAD_PIXELS] = 1
        for row_start_quad in self.ISSUE_QUADS[0]:
            for col_start_quad in self.ISSUE_QUADS[1]:
                mask[row_start_quad:(row_start_quad + 4), col_start_quad:(col_start_quad + 4)] = 1
        return mask
    
    @classmethod
    def get_size(cls) -> tuple:
        return (cls.ROWS, cls.COLS)
    
    @classmethod
    def get_rows(cls) -> int:
        return cls.ROWS
    
    @classmethod
    def get_columns(cls) -> int:
        return cls.COLS
        

class Stitcher:

    OVERLAP = 10
    DEAD_BAND_PIXELS = 37
    
    def __init__(self, rows: int = 0, columns: int = 0):
        self.weight_per = 0.5
        self.detector = Detector()
        self.stitch_rows = rows
        self.stitch_columns = columns
        self.size = (0, 0)
        
    def _determine_size(self):
        rows = self.stitch_rows * (self.detector.get_rows() + Stitcher.DEAD_BAND_PIXELS - Stitcher.OVERLAP) + Stitcher.OVERLAP
        columns = self.stitch_columns * (self.detector.get_columns() - Stitcher.OVERLAP) + Stitcher.OVERLAP
        self.size = (rows, columns)
        print(self.size)

    def load_data(self, directory: Path):
        if self.size == (0, 0):
            for file in directory.glob("*_*_*"):
                try:
                    stitch_row, stitch_col, _ = file.name.split("_")
                    stitch_row = int(stitch_row)
                    stitch_col = int(stitch_col)
                    if stitch_row > self.stitch_rows:
                        self.stitch_rows = stitch_row
                    if stitch_col > self.stitch_columns:
                        self.stitch_columns = stitch_col
                except ValueError:
                    pass
        print(f"Stitching an eiger_stitch {self.stitch_rows} {self.stitch_columns}")
        self._determine_size()

        data = np.zeros(self.size, dtype=np.float64)
        weight = np.zeros(self.size, dtype=np.float64)
        base_mask = np.logical_not(self.detector.calc_mask())
        
        for stitch_row in range(1, self.stitch_rows + 1):
            for stitch_col in range(1, self.stitch_columns + 1):
                for stitch_offset in range(1, 2 + 1):
                    
                    file_data = fabio.open(f"{stitch_row}_{stitch_col}_{stitch_offset}").data
                    mask = base_mask.copy()
                    mask[np.where(file_data == self.detector.MAX_INT)] = 0
                    file_data *= mask
                    start_row = (2 - stitch_offset) * Stitcher.DEAD_BAND_PIXELS + (self.stitch_rows - stitch_row) * (self.detector.ROWS - Stitcher.OVERLAP)
                    start_column = (stitch_col - 1) * (self.detector.COLS - Stitcher.OVERLAP)
                    data[start_row:(start_row + self.detector.ROWS), start_column:(start_column + self.detector.COLS)] += file_data
                    weight[start_row:(start_row + self.detector.ROWS), start_column:(start_column + self.detector.COLS)] += mask
        weight *= self.weight_per
        return data, weight
