# functions used by multiple scripts throughout the project

import os
import sys
import subprocess
import logging
import yaml
# project-specific imports
import utilities.parse_config as pc

# check if docker environment
# note: this method may not always be reliable
def is_docker_env():
    try:
        with open('/proc/1/cgroup', 'rt') as f:
            return 'docker' in f.read()
    except FileNotFoundError:
        return False

def is_gzipped(file_path):
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def pigz_available(quiet=False):
    try:
        result = subprocess.run(["pigz", "--version"], capture_output=True, text=True)
        if result.returncode == 0:
            version_line = result.stdout.split('\n')[0]
            version = version_line.split()[1]
            if not quiet:
                logging.info(f"pigz version: {version}")
            return True
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        if not quiet:
            logging.warning("pigz is not available. Falling back to gzip.")
        return False

# generic function to run a shell command
# returns stdout
def run_command(command):
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(result.stderr)
        sys.exit(1)
    return result.stdout

# load config as either yaml or sh
def load_config(quiet=False,config_file=''):
    if not config_file:
        if os.path.exists("config.yml") or os.path.exists("config.yaml"):
            config_file = "config.yml" if os.path.exists("config.yml") else "config.yaml"
        else:
            config_file = "config.sh"
    config = pc.parse_config(config_file, quiet)
    config["config_file"] = config_file
    return config

