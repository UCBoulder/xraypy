from setuptools import setup, find_packages
from pathlib import Path
import shutil
import XRDpy.package_params as package

detector_dir = package.directory / "Detectors"
notebooks_dir = package.directory / "Jupyter-Notebooks"
user_dir = package.directory / "User-Files"
if package.directory.exists():
    if detector_dir.exists():
        shutil.rmtree(detector_dir)
    if notebooks_dir.exists():
        shutil.rmtree(notebooks_dir)
else:
    package.directory.mkdir(parents=True, exist_ok=True)
detector_dir.mkdir(parents=True)
notebooks_dir.mkdir(parents=True)
if not user_dir.exists():
    user_dir.mkdir(parents=True)

for det_file in (Path("files") / "Detectors").glob("*.h5"):
    shutil.copyfile(det_file, detector_dir / det_file.name)
for notebook_file in (Path("files") / "Jupyter-Notebooks").glob("*.ipynb"):
    shutil.copyfile(notebook_file, notebooks_dir / notebook_file.name)

with open(package.directory / "config.txt", "w") as f:
    f.write(package.data_path.as_posix())

setup(
    packages=find_packages(include=['XRDpy', 'XRDpy.*']),
    scripts=["XRDpy/main.py",],
    py_modules=["XRDpy.transform"],
)