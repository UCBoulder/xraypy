from pathlib import Path

name = "XRDpy"

class Directory:
    DATA = (Path.home() / "DATA").resolve()
    home = (Path.home() / "Documents" / name).resolve()
    notebooks = home / "Jupyter-Notebooks"
    user = home / "User-Files"
    macros = user / "macros"


config = Directory.home / "config.toml"
