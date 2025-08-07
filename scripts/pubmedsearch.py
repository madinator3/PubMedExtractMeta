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
current_datetime = dt.datetime.now().strftime("%Y-%m-%d_%H-%M")
print("Current Date_Time:", current_datetime)

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


# DataFrame to store the extracted data
df = pd.DataFrame(columns=['PMID', 'Title', 'Abstract', 'Authors','Affiliations',
                           'Journal_Title', 'Journal_ISO', 'Keywords', 'Date_Pub_Year', 
                           'Date_Pub_Month', 'Date_Pub_Day', 'Date_Revise_Year', 
                           'Date_Revise_Month', 'Date_Revise_Day', 'Grants', 'URL', 'DOI'])

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
                for affiliation in author['AffiliationInfo']:
                    if 'Affiliation' in affiliation:
                        affiliations.append(affiliation['Affiliation'])
        affiliations = '; '.join(set(affiliations))

        # Extract journal information
        journal_title = record['MedlineCitation']['Article']['Journal']['Title']
        journal_iso = record['MedlineCitation']['Article']['Journal']['ISOAbbreviation'] if 'ISOAbbreviation' in record['MedlineCitation']['Article']['Journal'] else ''

        # Extract keywords
        keywords = ', '.join(keyword['DescriptorName'] for keyword in record['MedlineCitation']['MeshHeadingList']) if 'MeshHeadingList' in record['MedlineCitation'] else ''

        # Extract URL
        url = f"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}"

        # Extract DOI if available
        eid_list = record['MedlineCitation']['Article']['ELocationID']
        prefix = '10.'
        doi_list = []

        for doi in eid_list:
            if doi.startswith(prefix):
                doi_list.append(doi)
        
        doi = doi_list[0] if doi_list else ''  # Use the first DOI if available, otherwise an empty string

        # Extract Date published
        pubdate_list = record['MedlineCitation']['Article']['ArticleDate']
        for dmy in pubdate_list:
            if 'Year' in dmy and 'Month' in dmy and 'Day' in dmy:
                date_pub_yr = dmy['Year']
                date_pub_month = dmy['Month']
                date_pub_day = dmy['Day']
                break
            else:
                date_pub_yr = ''
                date_pub_month = ''
                date_pub_day = ''

        # Extract last revised date
        date_revise_yr = record['MedlineCitation']['DateRevised']['Year'] if 'DateRevised' in record['MedlineCitation'] and 'Year' in record['MedlineCitation']['DateRevised'] else ''
        date_revise_month = record['MedlineCitation']['DateRevised']['Month'] if 'DateRevised' in record['MedlineCitation'] and 'Month' in record['MedlineCitation']['DateRevised'] else ''
        date_revise_day = record['MedlineCitation']['DateRevised']['Day'] if 'DateRevised' in record['MedlineCitation'] and 'Day' in record['MedlineCitation']['DateRevised'] else ''
        
        # Extract Grants
        grants = []
        for gt in record['MedlineCitation']['Article'].get('GrantList', []):
            # Extracting GrantID, Acronym, and Agency
            if 'GrantID' in gt:
                gID = gt.get('GrantID',"NAN")  
            if 'Acronym' in gt:
                acronym = gt.get('Acronym',"NAN")  
            if 'Agency' in gt:
                agency = gt.get('Agency',"NAN")  
            if 'Country' in gt:
                country = gt.get('Country',"NAN")
            
            grant_info = f"{gID}_{acronym}_{agency}_{country}" 
            grants.append(grant_info)
  
        grants = '; '.join(set(grants))

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
            'Grants': [grants],
            'URL': [url],
            'DOI': [doi]
        })

        df = pd.concat([df, new_row], ignore_index=True)

print(df)

import requests

BASE_URL = "https://icite.od.nih.gov/api"
parameter_list = ['pmid', 'doi', 'citation_count', 'cited_by_clin', 'cited_by']
MAX_PMIDS_PER_REQUEST = 1000
sleep = 0.05
api_url = f'{BASE_URL}/pubs'


def get_icites(pmid_list=[], field_list=[], timeout=500):
    """
    Fetches publication data in bulk. If pmid_list exceeds 1000, it splits the requests and merges the results.
     Parameters:
            pmid_list: List of PMIDs,
            field_list: List of fields to be fetched,
            timeout: Timeout for the request, default is 500 seconds.
    Returns:
    A JSON object containing merged results of the requests.
    """
    BASE_URL = "https://icite.od.nih.gov/api"
    MAX_PMIDS_PER_REQUEST = 1000
    sleep = 0.05
    api_url = f'{BASE_URL}/pubs'

    responses_data = []  # To store response data from all batches

    # Split pmid_list into sublists of up to 1000 PMIDs each
    for i in range(0, len(pmid_list), MAX_PMIDS_PER_REQUEST):
        time.sleep(sleep)
        sub_pmid_list = pmid_list[i:i + MAX_PMIDS_PER_REQUEST]
        payload = {
            'pmids': ','.join(map(str, sub_pmid_list)),
            'fl': ','.join(field_list)
        }

        try:
            response = requests.get(api_url, params=payload, timeout=timeout)
            if response.status_code == 200:
                # Add this batch's results to the list of response data
                responses_data.extend(response.json().get('data', []))
            else:
                response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"Request failed: {e}")
            return None

    # Merge results from all batches and return
        
    return {
        'meta': {
            'pmids': ','.join(map(str, pmid_list)),
            'fl': ','.join(field_list)
        },
        'data': responses_data
    }

citations = get_icites(pmid_list=id_list, field_list=parameter_list)

df_cite = pd.DataFrame(columns=['PMID', 'DOI', 'Citation_Count', 'Cited_By_Clin', 'Cited_By'])

# Process each PubMed article in the response
for record in citations['data']:
    # Print the record in a formatted JSON style
    print(json.dumps(record, indent=4, default=str))  # default=str handles types JSON can't serialize like datetime

    # Extract PMID information
    pmid = record['pmid']

    # Extract DOI information
    doi = record['doi']

    # Extract DOI information
    citation_count = record['citation_count']

    # Extract Cited By Clin List
    cited_by_clin = record['cited_by_clin']

    # Extract Cited By List
    cited_by = record['cited_by']

    # Create a new row with the extracted data
    new_row = pd.DataFrame({
        'PMID': [pmid],
        'DOI': [doi],
        'Citation_Count': [citation_count],
        'Cited_By_Clin': [cited_by_clin],
        'Cited_By': [cited_by]
        })

    df_cite = pd.concat([df_cite, new_row], ignore_index=True)


df_final = pd.merge(df, df_cite, on='DOI', how='outer')

#TO DO: Change for loop to function, annotate function,  and number of citations extraction
#================ Export results =============================================================

from scripts.utils import save_data_to_file

# Save the DataFrame to a CSV file
save_data_to_file(df_final,  OUTPUT_PATH + "/" + current_datetime + "_" + "search_results_pubmed.csv")