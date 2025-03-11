import logging
import sys
import argparse

logging.basicConfig(
    level=logging.INFO,
    format="[{asctime}] [{levelname}] {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("output/reorient_fastq_parallel.log")
    ]
)

# config file is a shell script with variables defined as key=value pairs
# lines beginning with # are ignored as well as text after # on a given line
# variables can be defined as lists by enclosing values in parentheses (arrays)
def parse_config(file_path):
    config = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith('#'):
                key_value = line.split('=', 1)
                key = key_value[0].strip()
                value = key_value[1].split('#', 1)[0].strip() if len(key_value) > 1 else ''
                if value.startswith('(') and value.endswith(')'):
                    value = value[1:-1].split()
                config[key] = value
    logging.info(f"Configuration file parsed: {config}")
    return config

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse configuration and set variables.")
    parser.add_argument("--config", help="Path to the configuration file.")
    args = parser.parse_args()
