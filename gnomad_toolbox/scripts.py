"""Script to copy Jupyter notebooks and configurations to a user-specified directory."""

import json
import logging
import os
import shutil
import subprocess
import sys
from typing import Union

import click

CONFIG_FILE = os.path.expanduser("~/.gnomad_toolbox_config.json")
CONFIGS_DIR = "configs"
NOTEBOOKS_DIR = "notebooks"

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("gnomad_toolbox")
logger.setLevel(logging.INFO)


def load_config(config_file: str = CONFIG_FILE) -> dict:
    """
    Load configuration from a JSON file.

    :param config_file: Path to the configuration file.
    :return: The configuration dictionary.
    """
    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            return json.load(f)
    return {}


def save_config(config: dict, config_file: str = CONFIG_FILE) -> None:
    """
    Save configuration to a JSON file.

    :param config: The configuration dictionary.
    :param config_file: Path to the configuration file.
    :return: None.
    """
    with open(config_file, "w") as f:
        json.dump(config, f, indent=4)


def set_config(
    key: str, value: Union[str, dict], config_file: str = CONFIG_FILE
) -> None:
    """
    Add a key-value pair to the configuration file, supporting nested keys.

    :param key: The key to save, with nesting indicated by periods.
    :param value: The value to save.
    :param config_file: Path to the configuration file.
    :return: None.
    """
    config = load_config(config_file)
    keys = key.split(".")
    current = config

    # Traverse or create nested dictionaries.
    for k in keys[:-1]:
        current = current.setdefault(k, {})

    current[keys[-1]] = value
    save_config(config, config_file)


def get_config(key: str, config_file: str = CONFIG_FILE) -> Union[str, None]:
    """
    Retrieve a value from the configuration file by key, supporting nested keys.

    :param key: The key to retrieve, with nesting indicated by periods.
    :param config_file: Path to the configuration file.
    :return: The value associated with the key, or None if the key doesn't exist.
    """
    config = load_config(config_file)
    keys = key.split(".")
    current = config

    for k in keys:
        if k in current:
            current = current[k]
        else:
            return None

    return current


def copy_directory(src: str, dest: str, overwrite: bool = False) -> None:
    """
    Copy a directory to a destination.

    :param src: Source directory.
    :param dest: Destination directory.
    :param overwrite: Whether to overwrite if the destination exists.
    :return: None.
    """
    if os.path.exists(dest):
        if overwrite:
            shutil.rmtree(dest)
            logger.info("Overwriting existing directory: %s", dest)
        else:
            raise FileExistsError(f"Directory '{dest}' already exists.")
    shutil.copytree(src, dest)
    logger.info("Copied %s to %s", src, dest)


def copy_notebooks(destination: str, overwrite: bool = False) -> None:
    """
    Copy Jupyter notebooks and configurations to a user-specified directory.

    :param destination: The target directory.
    :param overwrite: Whether to overwrite existing files/directories.
    :return: None.
    """
    pkg_dir = os.path.dirname(__file__)
    notebook_dir = os.path.join(pkg_dir, NOTEBOOKS_DIR)
    config_dir = os.path.join(pkg_dir, CONFIGS_DIR)

    # Validate source directories.
    if not os.path.exists(notebook_dir):
        raise FileNotFoundError(f"No notebooks directory found at {notebook_dir}")
    if not os.path.exists(config_dir):
        raise FileNotFoundError(f"No configs directory found at {config_dir}")

    # Copy example jupyter notebooks.
    copy_directory(notebook_dir, destination, overwrite)

    # Copy jupyter configs.
    config_dest = os.path.join(destination, "jupyter_configs")
    copy_directory(config_dir, config_dest, overwrite)

    # Update gnomAD Toolbox configuration.
    set_config("notebook_dir", destination)
    logger.info("Default notebook directory set to: %s", destination)

    # Modify the Jupyter config file to set the notebook directory.
    jupyter_config = os.path.join(
        destination, "jupyter_configs/jupyter_notebook_config.json"
    )
    set_config("NotebookApp.notebook_dir", destination, jupyter_config)


@click.command()
@click.argument("destination", type=click.Path())
@click.option(
    "--overwrite", is_flag=True, help="Overwrite existing files if necessary."
)
def copy_notebooks_cli(destination: str, overwrite: bool) -> None:
    """
    CLI command to copy Jupyter notebooks and configurations.

    :param destination: Target directory for the notebooks and configs.
    :param overwrite: Whether to overwrite existing files.
    :return: None.
    """
    try:
        copy_notebooks(destination, overwrite)
        logger.info("Notebooks successfully copied to %s.", destination)
    except Exception as e:
        logger.error("Error during notebook copy: %s", e)


def run_jupyter_cli() -> None:
    """
    Launch Jupyter Lab or Notebook using the configured directory.

    :return: None.
    """
    notebook_dir = get_config("notebook_dir")
    if not notebook_dir:
        logger.error("No notebook directory configured. Run `copy-notebooks` first.")
        return

    if not os.path.exists(notebook_dir):
        logger.error("Configured notebook directory does not exist: %s", notebook_dir)
        return

    # Set Jupyter configuration directory.
    jupyter_config_dir = os.path.join(notebook_dir, "jupyter_configs")
    os.environ["JUPYTER_CONFIG_DIR"] = jupyter_config_dir
    logger.info("Launching Jupyter with config directory: %s", jupyter_config_dir)

    # Launch Jupyter.
    command = sys.argv[1] if len(sys.argv) > 1 else "lab"
    subprocess.run(["jupyter", command])
