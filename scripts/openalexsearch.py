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
# df = pd.json_normalize(ua_articles) 

# DataFrame to store the extracted data
df = pd.DataFrame(columns=['OAID', 'Title', 'Abstract', 'Authors','Author_ORCIDs','Affiliations',
                           'Journal_Title', 'Journal_ISO', 'Keywords', 'Date_Pub_Year', 
                           'Date_Pub_Month', 'Date_Pub_Day', 'Date_Revise_Year', 
                           'Date_Revise_Month', 'Date_Revise_Day', 'Grants', 'URL', 'DOI'])

# Fetch information for each record in the ua_articles list
for record in ua_articles:
    # Print the record in a formatted JSON style
    print(json.dumps(record, indent=4, default=str))  # default=str handles types JSON can't serialize like datetime

    # Extract title information
    title = record['title']

    # Extract open alex id and url
    url = record['id']
    oaid = record['id']

    # Extract doi
    if 'doi' in record and record['doi'] is not None:
        doi_url = record['doi']
        doi = doi_url.replace("https://doi.org/", "")
    else:
        doi = "NONE"

    # Extract author information
    authors = []
    authors_orcid = []
    affiliations = []
    affiliations_ror = []
    affiliations_oaid = []
    affiliations_country = []
    if 'authorships' in record and record['authorships'] is not None: 
        for author in record['authorships']['author']:
            # Extract author names
            if 'display_name' in author and author['display_name'] is not None:
                name = author['display_name']
            else: 
                name = "NONE"
            authors.append(name)

            # Extract author ORCIDs
            if 'orcid' in author and author['orcid'] is not None:
                orcid = author['orcid'] 
            else: 
                orcid = "NONE"
            authors_orcid.append(orcid)
        authors = '; '.join(authors)
        authors_orcid = '; '.join(authors_orcid)

    else: 
        authors = "NONE"
            
    
 # Create a new row with the extracted data
    new_row = pd.DataFrame({
        'OAID': [oaid],
        'Title': [title],
            # 'Abstract': [abstract],
        'Authors': [authors],
        'Author_ORCIDs': [authors],
            # 'Affiliations': [affiliations],
            # 'Journal_Title': [journal_title],
            # 'Journal_ISO': [journal_iso],
            # 'Keywords': [keywords],
            # 'Date_Pub_Year': [date_pub_yr],
            # 'Date_Pub_Month': [date_pub_month],
            # 'Date_Pub_Day': [date_pub_day],
            # 'Date_Revise_Year': [date_revise_yr],
            # 'Date_Revise_Month': [date_revise_month],
            # 'Date_Revise_Day': [date_revise_day],
            # 'Grants': [grants],
        'URL': [url],
        'DOI': [doi]
        })

    df = pd.concat([df, new_row], ignore_index=True)

print(df)

# #TO DO: parse the DataFrame to only include relevant columns; Add authors to query if possible
# #================ Export results =============================================================

from scripts.utils import save_data_to_file

# Save the DataFrame to a CSV file
save_data_to_file(df,  OUTPUT_PATH + "/" + current_datetime + "_" + "search_results_openalex.csv")
