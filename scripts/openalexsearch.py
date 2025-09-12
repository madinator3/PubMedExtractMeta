# Title: Searching PubMed
# Last Updated: 2025-07-30
# Description: This script searches PubMed for articles related to a specific topic using the Entrez Programming Utilities (E-utilities).
# Sources: 
# 


#=============== Libraries ===================================================================
import os
import sys
import pandas as pd
import json
import datetime as dt
import requests

#================ Project setup =============================================================
# Define the project root directory and add it to the system path
project_root = os.path.abspath(os.path.join(os.getcwd()))
sys.path.append(project_root)
print(project_root)

from scripts.configs import *
from scripts.utils import CONFIG, FILE_PATHS, PATH_ROOT, OUTPUT_PATH


# Checking if the necessary directories exist
print("Configuration:", CONFIG)
print("Path:", FILE_PATHS)
print("Source:", PATH_ROOT)
print("Output:", OUTPUT_PATH)
print("File", CONFIG["files"]["file_data_pubmed"])

# Get system date and time
current_datetime = dt.datetime.now().strftime("%Y-%m-%d_%H-%M")
print("Current Date_Time:", current_datetime)

#================ Import query =============================================================

from scripts.querysetup import full_query

print("Full Query:", full_query)


#================ Search PubMed for relevant records ========================================