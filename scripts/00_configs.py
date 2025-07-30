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
from scripts.utils import save_data_to_file


#=============== List of elements to export ==============================
# Lista de elementos exportados
__all__ = [
    "os",
    "sys",
    "time",
    "logging",
    "Path",
    "yaml",
    "pd",
    "Entrez",
    "save_data_to_file"
]