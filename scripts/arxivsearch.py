# Title: Searching PubMed
# Last Updated: 2025-07-30
# Description: This script searches PubMed for articles related to a specific topic using the Entrez Programming Utilities (E-utilities).
# Sources: 
# https://pypi.org/project/arxiv/
# https://info.arxiv.org/help/api/user-manual.html#query_details
# https://info.arxiv.org/help/api/examples/python_arXiv_parsing_example.txt
# https://mukhlisraza.medium.com/how-to-export-arxiv-papers-with-pythons-atom-api-e084c2970484





#=============== Libraries ===================================================================
import os
import sys
import pandas as pd
import feedparser
import urllib
import xml.etree.ElementTree as ET
import datetime as dt
import csv

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
print("File:", CONFIG["files"]["file_data_arxiv"])

# Get system date and time
current_datetime = dt.datetime.now().strftime("%Y-%m-%d_%H-%M")
print("Current Date_Time:", current_datetime)

#================ Import query =============================================================

from scripts.querysetup import arxivq

print("Full Query:", arxivq)


#================ Search Arxiv for relevant records ========================================
arxiv_api_url = "http://export.arxiv.org/api/query?search_query="

ns   = {'atom': 'http://www.w3.org/2005/Atom',
        'arxiv': 'http://arxiv.org/schemas/atom'}

# Set up Arxiv parameters
start = 0 # Start at the first record
max_tries = CONFIG['arxiv_config']['max_tries']  # Set the maximum number of attempts to fetch data from PubMed in case of failure.
sleep_between_tries = CONFIG['arxiv_config']['time_sleep']  # Set the time to sleep between queries in seconds
total_results = CONFIG['arxiv_config']['retmax']  # Set the maximum number of records to retrieve from PubMed in a single query. Max is 30000.
results_per_iteration = CONFIG['arxiv_config']['results_per_iteration']  # Number of results to fetch per iteration

# DataFrame to store the extracted data
df = pd.DataFrame(columns=['ARXIVID', 'Title', 'Abstract', 'Authors', 'Affiliations', 
                           'Date_published', 'URL'])

for i in range(start,total_results,results_per_iteration):
    
    print("Results %i - %i" % (i,i+results_per_iteration))

    params = '%s&start=%i&max_results=%i' % (arxivq,
                                            i,
                                            results_per_iteration)

    url = arxiv_api_url + params

    print("Query URL:", url)

    request = urllib.request.Request(url)

    with urllib.request.urlopen(request) as response:
        xml = response.read().decode('utf-8')

    root = ET.fromstring(xml)
    entries = root.findall('atom:entry', ns)

    # Fetch information for each record in the id_list
    for e in entries:
            eid       = e.find('atom:id', ns).text.strip().replace('http://arxiv.org/abs/', '')
            
            title     = e.find('atom:title', ns).text.strip().replace('\n', ' ')
            
            abstract     = e.find('atom:summary', ns).text.strip().replace('\n', ' ')
            
            authors = []
            for author in e.findall('atom:author', ns):
                name = author.find('atom:name', ns).text
                authors.append(name)
            authors = '; '.join(authors)

            affiliations = []
            for author in e.findall('atom:author', ns):
                aff = author.find('arxiv:affiliation', ns)
                if aff is not None and aff.text:
                    affiliations.append(aff.text.strip())
            affiliations = '; '.join(affiliations)

            published = e.find('atom:published', ns).text

            link      = e.find("atom:link[@type='text/html']", ns).attrib['href']

            # Create a new row with the extracted data
            new_row = pd.DataFrame({
                'ARXIVID': [eid],
                'Title': [title],
                'Abstract': [abstract],
                'Authors': [authors],
                'Affiliations': [affiliations],
                'Date_published': [published],
                'URL': [link],
            })

            df = pd.concat([df, new_row], ignore_index=True)
    
    print(f"Found {len(entries)} records from arXiv")

    # Rest before hitting the api again
    print('Sleeping for %i seconds' % sleep_between_tries)
    time.sleep(sleep_between_tries)

print(df)

#TO DO: EXTRACT DOI IF POSSIBLE AND SCRAPE AFFILIATIONS FROM PDFs
#================ Export results =============================================================

from scripts.utils import save_data_to_file

# Save the DataFrame to a CSV file
save_data_to_file(df,  OUTPUT_PATH + "/" + current_datetime + "_" + "search_results_arxiv.csv")
