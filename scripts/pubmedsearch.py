# Title: Searching PubMed
# Last Updated: 2025-07-30
# Description: This script searches PubMed for articles related to a specific topic using the Entrez Programming Utilities (E-utilities).
# Sources: 
# https://github.com/TLDWTutorials/PubmedAPI/blob/main/pubmed_api_in_python_2024.py
# https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
# https://biopython.org/docs/1.76/api/Bio.Entrez.html


#=============== Libraries ===================================================================
import os
import sys
import pandas as pd
import json
from Bio import Entrez
import datetime as dt

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
current_date = dt.date.today().strftime("%Y-%m-%d")
print("Current Date:", current_date)

#================ Import query =============================================================

from scripts.querysetup import full_query

print("Full Query:", full_query)


#================ Search PubMed for relevant records ========================================
# Set up Entrez API parameters
Entrez.api_key = CONFIG['pubmed_api_key'] # Set your API key for NCBI Entrez. If not set, only 3 queries per second are allowed. 10 queries per seconds otherwise with a valid API key.
Entrez.email = CONFIG['pubmed_api_email']  # Set your email to comply with NCBI's policy 
Entrez.max_tries = CONFIG['config']['max_tries']  # Set the maximum number of attempts to fetch data from PubMed in case of failure.
Entrez.sleep_between_tries = CONFIG['config']['time_sleep']  # Set the time to sleep between queries in seconds
retmax = CONFIG['config']['retmax']  # Set the maximum number of records to retrieve from PubMed in a single query. Max is 10000.


handle = Entrez.esearch(db='pubmed', retmax=retmax, term=full_query)
record = Entrez.read(handle)
id_list = record['IdList']
print(id_list)

# DataFrame to store the extracted data
df = pd.DataFrame(columns=['PMID', 'Title', 'Abstract', 'Authors','Affiliations',
                           'Journal_Title', 'Journal_ISO', 'Keywords', 'Date_Pub_Year', 
                           'Date_Pub_Month', 'Date_Pub_Day', 'Date_Revise_Year', 
                           'Date_Revise_Month', 'Date_Revise_Day', 'URL', 'DOI'])

# Fetch information for each record in the id_list
for pmid in id_list:
    handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
    records = Entrez.read(handle)

    # Process each PubMed article in the response
    for record in records['PubmedArticle']:
        # Print the record in a formatted JSON style
        print(json.dumps(record, indent=4, default=str))  # default=str handles types JSON can't serialize like datetime

        # Extract title information
        title = record['MedlineCitation']['Article']['ArticleTitle']

        # Extract abstract
        abstract = ' '.join(record['MedlineCitation']['Article']['Abstract']['AbstractText']) if 'Abstract' in record['MedlineCitation']['Article'] and 'AbstractText' in record['MedlineCitation']['Article']['Abstract'] else ''
        
        # Extract authors
        authors = ', '.join(author.get('LastName', '') + ' ' + author.get('ForeName', '') for author in record['MedlineCitation']['Article']['AuthorList'])
        
        # Extract author affiliations
        affiliations = []
        for author in record['MedlineCitation']['Article']['AuthorList']:
            if 'AffiliationInfo' in author and author['AffiliationInfo']:
                affiliations.append(author['AffiliationInfo'][0]['Affiliation'])
        affiliations = '; '.join(set(affiliations))

        # Extract journal information
        journal_title = record['MedlineCitation']['Article']['Journal']['Title']
        journal_iso = record['MedlineCitation']['Article']['Journal']['ISOAbbreviation'] if 'ISOAbbreviation' in record['MedlineCitation']['Article']['Journal'] else ''

        # Extract keywords
        keywords = ', '.join(keyword['DescriptorName'] for keyword in record['MedlineCitation']['MeshHeadingList']) if 'MeshHeadingList' in record['MedlineCitation'] else ''

        # Extract URL
        url = f"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}"

        # Extract DOI if available
        doi = record['MedlineCitation']['Article']['ELocationID'][1] if 'ELocationID' in record['MedlineCitation']['Article']['ELocationID'] else ''

        # Extract Date published
        date_pub_yr = record['MedlineCitation']['Article']['ArticleDate']['Year'] if 'ArticleDate' in record['MedlineCitation']['Article'] and 'Year' in record['MedlineCitation']['Article']['ArticleDate'] else ''
        date_pub_month = record['MedlineCitation']['Article']['ArticleDate']['Month'] if 'ArticleDate' in record['MedlineCitation']['Article'] and 'Month' in record['MedlineCitation']['Article']['ArticleDate'] else ''
        date_pub_day = record['MedlineCitation']['Article']['ArticleDate']['Day'] if 'ArticleDate' in record['MedlineCitation']['Article'] and 'Day' in record['MedlineCitation']['Article']['ArticleDate'] else ''

        # Extract last revised date
        date_revise_yr = record['MedlineCitation']['DateRevised']['Year'] if 'DateRevised' in record['MedlineCitation'] and 'Year' in record['MedlineCitation']['DateRevised'] else ''
        date_revise_month = record['MedlineCitation']['DateRevised']['Month'] if 'DateRevised' in record['MedlineCitation'] and 'Month' in record['MedlineCitation']['DateRevised'] else ''
        date_revise_day = record['MedlineCitation']['DateRevised']['Day'] if 'DateRevised' in record['MedlineCitation'] and 'Day' in record['MedlineCitation']['DateRevised'] else ''
        # Extract number of citations if available


        # Create a new row with the extracted data
        new_row = pd.DataFrame({
            'PMID': [pmid],
            'Title': [title],
            'Abstract': [abstract],
            'Authors': [authors],
            'Affiliations': [affiliations],
            'Journal_Title': [journal_title],
            'Journal_ISO': [journal_iso],
            'Keywords': [keywords],
            'Date_Pub_Year': [date_pub_yr],
            'Date_Pub_Month': [date_pub_month],
            'Date_Pub_Day': [date_pub_day],
            'Date_Revise_Year': [date_revise_yr],
            'Date_Revise_Month': [date_revise_month],
            'Date_Revise_Day': [date_revise_day],
            'URL': [url],
            'DOI': [doi]
        })

        df = pd.concat([df, new_row], ignore_index=True)

print(df)

#TO DO: Change for loop to function, annotate function, add in DOI extraction, and number of citations extraction, and error handling

#================ Export results =============================================================

from scripts.utils import save_data_to_file

# Save the DataFrame to a CSV file
save_data_to_file(df,  OUTPUT_PATH + "/" + current_date + "_" + "search_results_pubmed.csv")