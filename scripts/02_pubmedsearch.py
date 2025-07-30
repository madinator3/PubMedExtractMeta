# Title: Searching PubMed
# Last Updated: 2025-07-30
# Description: This script searches PubMed for articles related to a specific topic using the Entrez Programming Utilities (E-utilities).
# Source: https://github.com/TLDWTutorials/PubmedAPI/blob/main/pubmed_api_in_python_2024.py

#=============== Libraries ==============
import os
import sys
import pandas as pd
import json
from Bio import Entrez


#================ Project setup ==============
# Define the project root directory and add it to the system path
project_root = os.path.abspath(os.path.join(os.getcwd(), ".."))
sys.path.append(project_root)
print(project_root)

from scripts.configs import *
from scripts.utils import  CONFIG, FILE_PATHS, PATH_ROOT, OUTPUT_PATH

config.load_config(os.path.join(PATH_ROOT, "config.yaml"))

# Acessa as variáveis globais
print("Configuration:", CONFIG)
print("Path:", FILE_PATHS)
print("Source:", PATH_ROOT)
print("Output:", OUTPUT_PATH)
print("File", CONFIG["files"]["file_data_pubmed"])