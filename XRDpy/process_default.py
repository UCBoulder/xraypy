import sys
import yaml
from pathlib import Path
import cuxrd.package_params as package

class ConfigHandler:
    def __init__(self, args) -> None:
        self.user_config_dir = Path.home() / package.name
        self.config = ConfigHandler.load()
        self.iterator = iter(args)
        self.flag = None
        self.item = None
        self.start()
    
    def next(self):
        self.flag = next(self.iterator).lower()
        self.handle_flag()

    def handle_flag(self):
        if "=" in self.flag:
            flag, param = self.flag.split("=")
            if flag in ["--alpha", "--incident", "--in", "-i"]:
                incident = float(param)
                self.config["incident"] = incident
                print("Changing default incident angle to {} degrees".format(incident))
            elif flag in ["--tilt", "-t"]:
                tilt = float(param)
                self.config["tilt"] = tilt
                print("Changing default tilt to {} degrees".format(tilt))
            elif flag in ["--expose", "--exposure", "--exposuretime", "--time", '-e']:
                exposure_time = int(param)
                self.config["exposure"] = exposure_time
                print("Changing default exposure time to {} s".format(exposure_time))
            elif flag in ["--process", "-p"]:
                self.config["process"] = param
                print("Changing default process to {}".format(param))
        elif self.flag in ["--help", "-h"]:
            with open("help_default.txt", "r") as f:
                print(f.read())
        elif self.flag in ["--reset", "--remove", "-r"]:
            print("Resetting config file to program defaults.")
            self.config = ConfigHandler.load_default()

    def save_user_config(self):
        print("Saving new config file: {}".format(self.user_config_dir / package.name))
        if not self.user_config_dir.exists():
            self.user_config_dir.mkdir()
        with open(self.user_config_dir / package.config_name, 'w') as f:
            yaml.dump(self.config, f)

    @staticmethod  
    def load():
        user_config = Path.home() / package.name / package.config_name
        if user_config.is_file():
            with open(user_config, 'r') as f:
                config = yaml.safe_load(f)
        else:
            config = ConfigHandler.load_default()
        return config
    
    @staticmethod
    def load_default():
        config = package.default_config
        return config
    
    def start(self):
        while True:
            try:
                self.next()
            except StopIteration:
                break
        self.save_user_config()
    

def main(args:list=None):
    if args is None:
        args = sys.argv[1:]
    handler = ConfigHandler(args)
    