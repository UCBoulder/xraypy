from pathlib import Path
import tomllib
import sys

with open(Path(sys.prefix) / "Lib/site-packages/xraypy/config.toml", "rb") as file:
    config = tomllib.load(file)["paths"]

class Directory:
    DATA = (Path.home() / config["data"]).resolve()
    home = (Path.home() / "Documents" / config["home"]).resolve()
    notebooks = home / config["notebooks"]
    user = home / config["user"]
    macros = user / config["macros"]

config = Directory.home / "config.toml"
