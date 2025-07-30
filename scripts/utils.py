# Title: Utils
# Last Updated: 2025-07-30
# Description: This script contains utils functions for searching PubMed using the Entrez Programming Utilities (E-utilities). 
# Source: https://github.com/SoftwareImpacts/SIMPAC-2024-316

#=============== Libraries ==============
import yaml
import os
import time
import pandas as pd
import logging
from Bio import Entrez
from pathlib import Path

#=============== Global Var ==============
CONFIG = None
FILE_PATHS = None
PATH_ROOT = None

#=============== Functions ==============
def load_config(config_path):
    """
    Load configuration from a YAML file.
    
    :param config_path: Path to the YAML configuration file.
    :return: Dictionary containing configuration data.
    """
    try:
        with open(config_path, "r") as config_file:
            return yaml.safe_load(config_file)
    
    except FileNotFoundError:
        # Log error if configuration file is not found
        raise FileNotFoundError(f"Configuration file not found at {config_path}. Please provide a valid path.")
    
    except yaml.YAMLError as e:
        # Log error if there is an issue parsing the YAML file
        raise ValueError(f"Error parsing YAML file: {e}")



def get_output_file_paths(file_config):
    """
    Get file paths for outputs based on configuration.
    
    :param file_config: Dictionary with file paths.
    :return: Dictionary of file paths.
    """
    try:
        return {
            key: Path(file_config[key]) for key in file_config
        }
    except KeyError as e:
        logging.error(f"Missing file path in configuration: {e}")
        raise KeyError(f"Missing file path in configuration: {e}")
    


def save_data_to_file(df, file_path):
    """
    Save the data from a DataFrame to a CSV file. 
    If the file does not exist, it will be created with headers.
    If the file already exists, the data will be appended without headers.

    Parameters:
        df (pd.DataFrame): The DataFrame containing the data to be saved.
        file_path (str): The full path of the CSV file.

    Error Handling:
        - Logs any exception that occurs while attempting to save the file.
    
    Example:
        save_data_to_file(my_dataframe, "data.csv")
    """
    try:
        if not os.path.exists(file_path):
            df.to_csv(file_path, sep='|', index=False)
        else:
            df.to_csv(file_path, sep='|', mode='a', header=False, index=False)
    
    except Exception as e:
        logging.error(f"Error saving file: {e}")




def configure_logging(path_root, logging_config):
    """
    Configures logging based on the provided dictionary.

    Args:
        logging_config (dict): A dictionary containing logging configuration.

    Raises:
        ValueError: If logging configuration fails.
    """
    try:
        # Validate the logging configuration dictionary
        required_keys = ["log_file", "level", "format", "date_format"]
        for key in required_keys:
            if key not in logging_config:
                raise KeyError(f"Missing required logging configuration key: '{key}'")

        # Extract configuration values
        log_file = logging_config["log_file"]
        level = logging_config["level"].upper()

        log_dir1 = os.path.join(path_root, log_file)
        log_dir2 = os.path.normpath(log_dir1)
                
        # Configure logging
        logging.basicConfig(
            filename=log_dir2,
            level=getattr(logging, level, logging.INFO),
            format=logging_config["format"],
            datefmt=logging_config["date_format"]
        )
        # Test logging setup
        logging.info("Logging configuration successfully loaded.")

    except Exception as e:
        logging.error(f"Error configuring logging: {e}")
        raise ValueError(f"Error configuring logging: {e}")