from setuptools import setup, find_packages
from pathlib import Path
import shutil
import sys

sys.path.append(str(Path(__file__).parent / "XRDpy"))
import package_params as package

with open ("README.md", "r") as f:
    long_description = f.read()

if package.Directory.home.exists():
    if package.Directory.notebooks.exists():
        shutil.rmtree(package.Directory.notebooks)
else:
    package.Directory.home.mkdir(parents=True, exist_ok=True)
package.Directory.notebooks.mkdir(parents=True)

if not package.Directory.user.exists():
    package.Directory.user.mkdir(parents=True)
for notebook_file in Path("files/Jupyter-Notebooks").glob("*.ipynb"):
    shutil.copyfile(notebook_file, package.Directory.notebooks / notebook_file.name)

shutil.copyfile(Path("files/config.toml"), package.config)

# See pyproject.toml for metadata
setup(
    long_description=long_description,
    long_description_content_type="text/markdown",
)