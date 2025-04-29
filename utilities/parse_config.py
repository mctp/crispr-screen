import logging
import sys
import argparse
import os
import yaml

# config file is a shell script (.sh) or YAML (.yaml/.yml) 
# if shell script
#     variables are defined as key=value pairs
#     lines beginning with # are ignored as well as text after # on a given line
#     variables can be defined as lists by enclosing values in parentheses (arrays)
#     boolean values are defined as TRUE or FALSE but case is technically ignored
# if YAML
#     variables are defined as key: value pairs
#     lists are defined as key:
#       - value1
#       - value2
#     boolean values are defined as true or false

def parse_config(file_path, quiet=False):
    config = {}
    file_extension = os.path.splitext(file_path)[1].lower()

    if file_extension == '.sh':
        with open(file_path, 'r') as file:
            for line in file:
                if line.strip() and not line.startswith('#'):
                    key_value = line.split('=', 1)
                    key = key_value[0].strip()
                    value = key_value[1].split('#', 1)[0].strip() if len(key_value) > 1 else ''
                    if value.startswith('(') and value.endswith(')'):
                        value = value[1:-1].split()
                    elif value.upper() == "TRUE":
                        value = True
                    elif value.upper() == "FALSE":
                        value = False
                    config[key] = value
    elif file_extension == '.yaml' or file_extension == '.yml':
        with open(file_path, 'r') as file:
            config = yaml.safe_load(file)
    else:
        raise ValueError("Unsupported file type. Only .sh and .yaml/.yml files are supported.")

    if not quiet:
        logging.info(f"Configuration file parsed: {file_path}")
    return config

if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format="[{asctime}] [{levelname}] {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler("output/parse_config.log")
        ]
    )
    parser = argparse.ArgumentParser(description="Parse configuration and set variables.")
    parser.add_argument("--config", help="Path to the configuration file.")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress logging output.")
    args = parser.parse_args()

    parse_config(args.config, args.quiet)
