# Title: Configs
# Last Updated: 2025-07-30
# Description: This script contains configuration settings for the PubMed API client.
# Source: https://github.com/SoftwareImpacts/SIMPAC-2024-316

#=============== Libraries ==============
import os
import sys
import time
import logging
from pathlib import Path
import yaml
import pandas as pd
from Bio import Entrez



#=============== Importing modules internal to the project ==============
from scripts.utils import initialize_environment, save_data_to_file


#=============== List of elements to export ==============================
__all__ = [
    "os",
    "sys",
    "time",
    "logging",
    "Path",
    "yaml",
    "pd",
    "Entrez",
    "initialize_environment",
    "save_data_to_file"
]

#=============== Initialize the environment ==============================
initialize_environment()