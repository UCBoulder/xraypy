from setuptools import setup, find_packages
from pathlib import Path
import shutil
import tomllib

this_path = Path(__file__).parent.resolve()

def create_directories(config):
    home_path = Path.home() / "Documents" / config["home"]
    notebooks_path = home_path / config["notebooks"]
    user_path = home_path / config["user"]
    macros_path = user_path / config["macros"]

    if not home_path.exists():
        home_path.mkdir(parents=True, exist_ok=True)
    if notebooks_path.exists():
        shutil.rmtree(notebooks_path)
    notebooks_path.mkdir(parents=True)

    if not user_path.exists():
        user_path.mkdir(parents=True)
    for notebook_file in Path(this_path / "xraypy/Jupyter-Notebooks").glob("*.ipynb"):
        shutil.copyfile(notebook_file, notebooks_path / notebook_file.name)
    if not macros_path.exists():
        macros_path.mkdir(parents=True)
    shutil.copyfile(this_path / "xraypy/config.toml", home_path / "config.toml")

with open (this_path / "README.md", "r") as f:
    long_description = f.read()

with open(this_path / "xraypy/config.toml", "rb") as f:
    paths = tomllib.load(f)["paths"]

create_directories(paths)

# See pyproject.toml for metadata
setup(
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
)