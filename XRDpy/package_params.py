from pathlib import Path

name = "XRDpy"
directory = Path.home() / "Documents" / name
data_path = Path.home() / "DATA"
config_name = "config.txt"
default_config = {"data-file": "raw-stitched-data",
                  "weight-file": "stitched-exposure-time",
                  "process": "default",
                  "incident": None,
                  "tilt": None,
                  "exposure": None,
                  "plot": True}