# Title: Searching Open Alex
# Last Updated: 2025-09-12
# Description: This script searches Open Alexfor articles related to a specific topics.
# Sources: 
# https://docs.openalex.org/how-to-use-the-api/api-overview
# https://github.com/J535D165/pyalex
# https://ua-libraries-research-data-services.github.io/UALIB_ScholarlyAPI_Cookbook/src/python/openalex.html
# https://api.openalex.org/works?per-page=50&page=2



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
df = pd.DataFrame(columns=['OAID', 'Title', 'Abstract', 'Authors','Author_ORCIDs',
                           'Affiliations','Affiliations_ROR','Affiliations_OAID','Affiliations_Country',
                           'Journal_Title', 'Journal_ISO', 'Journal_Host', 'Keywords', 'Date_Published', 'Year_Published',
                           'Work_Type','OpenAccess_Status','Publication_Status','Retraction_Status', 'Cited_By', 'URL', 'DOI'])

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

    if 'authorships' in record and record['authorships']: 
        for authorship in record['authorships']:
            # Define authorship tables in a variable for easier access
            author = authorship['author']
            institution = authorship['institutions']

            # Extract author names
            author = authorship['author']
            if 'display_name' in author and author['display_name']:
                name = author['display_name']
            else: 
                name = "NONE"
            authors.append(name)

            # Extract author ORCIDs
            if 'orcid' in author and author['orcid']:
                orcid = author['orcid'] 
            else: 
                orcid = "NONE"
            authors_orcid.append(orcid)


        authors = '; '.join(authors)
        authors_orcid = '; '.join(authors_orcid)

    else: 
        authors = "NONE"
        authors_orcid = "NONE"
            
    affiliations = []
    affiliations_ror = []
    affiliations_oaid = []
    affiliations_country = []
    if 'authorships' in record and record['authorships'] is not None: 
        for places in record['authorships']:
            for institution in places['institutions']:
                # Extract affiliation names
                if 'display_name' in institution and institution['display_name'] is not None:
                    aff = institution['display_name'] 
                else:
                    aff = "NONE"
                affiliations.append(aff)
                
                # Extract ror
                if 'ror' in institution and institution['ror']:
                    ror = institution['ror'] 
                else:
                    ror = "NONE"
                affiliations_ror.append(ror)

                # Extract openalex institution id
                if 'id' in institution and institution['id']:
                    int_oaid = institution['id'] 
                else:
                    int_oaid = "NONE"
                affiliations_oaid.append(int_oaid)

                # Extract institution country
                if 'country_code' in institution and institution['country_code']:
                    country = institution['country_code']         
                else:
                    country = "NONE"    
                affiliations_country.append(country)

        affiliations = '; '.join(affiliations)
        affiliations_ror = '; '.join(affiliations_ror)  
        affiliations_oaid = '; '.join(affiliations_oaid)
        affiliations_country = '; '.join(affiliations_country)

    else: 
        affiliations = "NONE"
        affiliations_ror = "NONE"
        affiliations_oaid = "NONE"
        affiliations_country = "NONE"   

    if affiliations == '':
        affiliations = "NONE"
    if affiliations_ror == '':
        affiliations_ror = "NONE"
    if affiliations_oaid == '':
        affiliations_oaid = "NONE"
    if affiliations_country == '':
        affiliations_country = "NONE"
    
    # Extract year published
    pub_yr = record['publication_year']

    # Extract date published
    date_pub = record['publication_date']

    # Extract Journal title, ISSN, and host
    if 'primary_location' in record and record['primary_location'] is not None:
        if 'source' in record['primary_location'] and record['primary_location']['source'] is not None:
            # Extract jounral title
            if 'display_name' in record['primary_location']['source'] and record['primary_location']['source']['display_name'] is not None:
                journal_title = record['primary_location']['source']['display_name']
            else:
                journal_title = "NONE"
            # Extract journal ISSN
            if 'issn_l' in record['primary_location']['source'] and record['primary_location']['source']['issn_l'] is not None:
                journal_iso = record['primary_location']['source']['issn_l']
            else:
                journal_iso = "NONE"
            # Extract journal host
            if 'host_organization' in record['primary_location']['source'] and record['primary_location']['source']['host_organization'] is not None:
                journal_host = record['primary_location']['source']['host_organization_name']
            else:
                journal_host = "NONE"
        else:
            journal_title = "NONE"
            journal_iso = "NONE"
            journal_host = "NONE"
    else:
        journal_title = "NONE"
        journal_iso = "NONE"
        journal_host = "NONE"

    # Extract work type
    if 'type' in record and record['type'] is not None:
        work_type = record['type']  
    else: 
        work_type = "NONE"

    # Extract open access status
    if 'open_access' in record and record['open_access'] is not None:
        if 'is_oa' in record['open_access'] and record['open_access']['is_oa'] is not None:
            is_oa = record['open_access']['is_oa'] 
        else:
            is_oa = "NONE"
    else:
        is_oa = "NONE"

    # Extract Cited by count
    if 'cited_by_count' in record and record['cited_by_count'] is not None:
        cited_by = record['cited_by_count']
    else:
        cited_by = 0

    # Extract Publication status
    if 'primary_location' in record and record['primary_location'] is not None:
        if 'is_published' in record['primary_location'] and record['primary_location']['is_published'] is not None:
            pub_status = record['primary_location']['is_published']
        else:
            pub_status = "NONE"
    else:
        pub_status = "NONE"

    # Extract Retraction status
    if 'is_retracted' in record and record['is_retracted'] is not None:
        retraction_status = record['is_retracted']
    else:
        retraction_status = "NONE"

    # Extract keywords
    if 'concepts' in record and record['concepts'] is not None:
        keywords = []
        for concept in record['concepts']:
            if 'display_name' in concept and concept['display_name'] is not None:
                keywords.append(concept['display_name'])
            else:
                pub_status = "NONE"
    else:
        pub_status = "NONE"

    keywords = '; '.join(keywords)

    # Extract abstract 
    def invert_abstract(inv_index):
        """Invert OpenAlex abstract index. Function from pyAlex.
        https://github.com/J535D165/pyalex/blob/main/pyalex/api.py

        Parameters
        ----------
        inv_index : dict
            Inverted index of the abstract.

        Returns
        -------
        str
            Inverted abstract.
        """

        if inv_index is not None:
            l_inv = [(w, p) for w, pos in inv_index.items() for p in pos]
            return " ".join(map(lambda x: x[0], sorted(l_inv, key=lambda x: x[1])))
        else:
            return "NONE"
    
    abstract = invert_abstract(record['abstract_inverted_index'])

    # Create a new row with the extracted data
    new_row = pd.DataFrame({
        'OAID': [oaid],
        'Title': [title],
        'Authors': [authors],
        'Author_ORCIDs': [authors_orcid],
        'Affiliations': [affiliations],
        'Affiliations_ROR': [affiliations_ror],
        'Affiliations_OAID': [affiliations_oaid],
        'Affiliations_Country': [affiliations_country],
        'Date_Published': [date_pub],
        'Year_Published': [pub_yr],
        'Journal_Title': [journal_title],
        'Journal_ISO': [journal_iso],
        'Journal_Host': [journal_host],
        'Work_Type': [work_type],
        'OpenAccess_Status': [is_oa],
        'Cited_By': [cited_by],
        'Publication_Status': [pub_status],
        'Retraction_Status': [retraction_status],
        'Keywords': [keywords],
        'Abstract': [abstract],
        'URL': [url],
        'DOI': [doi]
        })

    df = pd.concat([df, new_row], ignore_index=True)

print(df)

# #TO DO: parse the DataFrame to extract Grants, Keywords, Abstract (look at pyAlex for this one); Add authors to query if possible; Sanity/Verify check export

# #================ Export results =============================================================

from scripts.utils import save_data_to_file

# Save the DataFrame to a CSV file
save_data_to_file(df,  OUTPUT_PATH + "/" + current_datetime + "_" + "search_results_openalex.csv")
