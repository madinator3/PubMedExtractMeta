# Title: Searching Open Alex
# Last Updated: 2025-09-12
# Description: This script searches Open Alexfor articles related to a specific topics.
# Sources: 
# https://docs.openalex.org/how-to-use-the-api/api-overview
# https://github.com/J535D165/pyalex
# https://ua-libraries-research-data-services.github.io/UALIB_ScholarlyAPI_Cookbook/src/python/openalex.html



#=============== Libraries ===================================================================
import os
import sys
import pandas as pd
import json
import datetime as dt
import requests
from time import sleep

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

from scripts.querysetup import oaq

print("Full Query:", oaq)


#================ Search OpenAlex for relevant records ========================================
# Parameters for obtaining UA publications
BASE_URL = 'https://api.openalex.org/'
endpoint = 'works'
filters = oaq

# The cursor to the next page is given in the metadata of the response
# To get the first page, we set the cursor parameter to '*'
cursor = '*'

ua_articles = []

# When the end of the results has been reached, the 'next_cursor' variable will be null
# This while loop iterates through all pages of the results
while cursor is not None:
    params = {
        'filter': filters,
        'cursor': cursor,
        'per-page': CONFIG['openalex_config']['results_per_iteration'],
        'mailto': CONFIG['openalex_config']['api_email']
    }
    page_data = requests.get(BASE_URL + endpoint, params=params).json()

    # Set cursor to the next page
    cursor = page_data['meta']['next_cursor']
    
    # Add results to the ua_articles list
    ua_articles.extend(page_data['results'])
    
    # Wait 3 seconds between requests to follow OpenAlex's rate limit
    sleep(CONFIG['openalex_config']['time_sleep'])

# Display number of UA articles in database
print(len(ua_articles))

# Convert the list of articles to a DataFrame
df = pd.json_normalize(ua_articles) 


#TO DO: annotate function, move citation information function into utils and its configs in configs.yaml
#================ Export results =============================================================

from scripts.utils import save_data_to_file

# Save the DataFrame to a CSV file
save_data_to_file(df,  OUTPUT_PATH + "/" + current_datetime + "_" + "search_results_openalex.csv")
