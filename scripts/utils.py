# Title: Utils
# Last Updated: 2025-07-30
# Description: This script contains utils functions for searching PubMed using the Entrez Programming Utilities (E-utilities). 
# Source: https://github.com/SoftwareImpacts/SIMPAC-2024-316

#=============== Libraries ==============
import yaml
import os
import sys
import time
import pandas as pd
import logging
from Bio import Entrez
from pathlib import Path

#=============== Global Var ==============
CONFIG = None
FILE_PATHS = None
PATH_ROOT = None
OUTPUT_PATH = None

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
            df.to_csv(file_path, index=False)
        else:
            df.to_csv(file_path, mode='a', header=False, index=False)
    
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






def setup_directories(path_root, directories):
    """
    Create necessary directories if they don't exist.
    
    :param directories: List of directory paths to create.
    """
    try:
        # Iterate through the directory configurations
        for dir_type, dir_path in directories.items():
            
            dir_path = os.path.join(path_root, dir_path)
            dir_path = os.path.normpath(dir_path)

            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
                logging.info(f"Created directory: {dir_path}")
            else:
                logging.info(f"Directory already exists: {dir_path}")

    except Exception as e:
        logging.error(f"Error while processing YAML: {e}")
        print(f"Error while processing YAML: {e}")




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
    




def initialize_environment():
    """
    Initialize the application environment by loading configuration, 
    setting up directories, and configuring logging.
    """
    global CONFIG, FILE_PATHS, PATH_ROOT, OUTPUT_PATH

    # Root path for the application
    PATH_ROOT = os.path.abspath(os.path.join(os.getcwd()))
    sys.path.append(PATH_ROOT)

    config_path = os.path.join(PATH_ROOT, 'config.yaml')
    CONFIG = load_config(config_path)
    
    configure_logging(PATH_ROOT, CONFIG["logging"])

    # Define email and API key for Entrez
    Entrez.email = CONFIG['pubmed_config']["api_email"]
    Entrez.api_key = CONFIG['pubmed_config']["api_key"]

    setup_directories(PATH_ROOT, CONFIG["directories"])

    FILE_PATHS = get_output_file_paths(CONFIG["files"])
    logging.info(f"Environment initialized successfully.{PATH_ROOT}")

    #path input process
    #DATA_INPUT_PATH_TMP = CONFIG["directories"]['input']
    #INPUT_PATH = os.path.normpath(os.path.join(PATH_ROOT, DATA_INPUT_PATH_TMP))

    #path output process
    RESEARCH_OUTPUT_PATH_TMP = CONFIG["directories"]['output'] 
    OUTPUT_PATH     = os.path.normpath(os.path.join(PATH_ROOT, RESEARCH_OUTPUT_PATH_TMP))
