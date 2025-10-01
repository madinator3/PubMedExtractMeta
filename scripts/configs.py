# Title: Configs
# Last Updated: 2025-07-30
# Description: This script contains configuration settings for the PubMed API client.
# Source: https://github.com/SoftwareImpacts/SIMPAC-2024-316

#=============== Libraries ==============
# Global Libraries
import os
import sys
import time
import logging
from pathlib import Path
import yaml
import pandas as pd

# PubMed Search Libraries
from Bio import Entrez

# Arxiv Search Libraries

# OpenAI Libraries



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