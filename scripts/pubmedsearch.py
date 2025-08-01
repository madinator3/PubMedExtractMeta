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





#================ Import query ==============

from scripts.querysetup import full_query

print("Full Query:", full_query)





#================ Search PubMed for relevant records
Entrez.api_key = CONFIG['pubmed_api_key'] # Set your API key for NCBI Entrez. If not set, only 3 queries per second are allowed. 10 queries per seconds otherwise with a valid API key.
Entrez.email = CONFIG['email']  # Set your email to comply with NCBI's policy 
Entrez.max_tries = CONFIG['max_tries']  # Set the maximum number of attempts to fetch data from PubMed in case of failure.
Entrez.sleep_between_tries = CONFIG['time_sleep']  # Set the time to sleep between queries in seconds
retmax = CONFIG['retmax']  # Set the maximum number of records to retrieve from PubMed in a single query. Max is 10000.


handle = Entrez.esearch(db='pubmed', retmax=retmax, term=full_query, ma)
record = Entrez.read(handle)
id_list = record['IdList']

# DataFrame to store the extracted data
df = pd.DataFrame(columns=['PMID', 'Title', 'Abstract', 'Authors', 'Journal', 'Keywords', 'URL', 'Affiliations'])

# Fetch information for each record in the id_list
for pmid in id_list:
    handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
    records = Entrez.read(handle)

    # Process each PubMed article in the response
    for record in records['PubmedArticle']:
        # Print the record in a formatted JSON style
        print(json.dumps(record, indent=4, default=str))  # default=str handles types JSON can't serialize like datetime
        title = record['MedlineCitation']['Article']['ArticleTitle']
        abstract = ' '.join(record['MedlineCitation']['Article']['Abstract']['AbstractText']) if 'Abstract' in record['MedlineCitation']['Article'] and 'AbstractText' in record['MedlineCitation']['Article']['Abstract'] else ''
        authors = ', '.join(author.get('LastName', '') + ' ' + author.get('ForeName', '') for author in record['MedlineCitation']['Article']['AuthorList'])
        
        affiliations = []
        for author in record['MedlineCitation']['Article']['AuthorList']:
            if 'AffiliationInfo' in author and author['AffiliationInfo']:
                affiliations.append(author['AffiliationInfo'][0]['Affiliation'])
        affiliations = '; '.join(set(affiliations))

        journal = record['MedlineCitation']['Article']['Journal']['Title']
        keywords = ', '.join(keyword['DescriptorName'] for keyword in record['MedlineCitation']['MeshHeadingList']) if 'MeshHeadingList' in record['MedlineCitation'] else ''
        url = f"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}"

        new_row = pd.DataFrame({
            'PMID': [pmid],
            'Title': [title],
            'Abstract': [abstract],
            'Authors': [authors],
            'Journal': [journal],
            'Keywords': [keywords],
            'URL': [url],
            'Affiliations': [affiliations]
        })

        df = pd.concat([df, new_row], ignore_index=True)

