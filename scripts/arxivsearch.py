# Title: Searching PubMed
# Last Updated: 2025-07-30
# Description: This script searches PubMed for articles related to a specific topic using the Entrez Programming Utilities (E-utilities).
# Sources: 
# https://pypi.org/project/arxiv/
# https://info.arxiv.org/help/api/user-manual.html#query_details
# https://info.arxiv.org/help/api/examples/python_arXiv_parsing_example.txt




#=============== Libraries ===================================================================
import os
import sys
import pandas as pd
import feedparser
import urllib
import json
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
print("File:", CONFIG["files"]["file_data_arxiv"])

# Get system date and time
current_datetime = dt.datetime.now().strftime("%Y-%m-%d_%H-%M")
print("Current Date_Time:", current_datetime)

#================ Import query =============================================================

from scripts.querysetup import arxivq

print("Full Query:", arxivq)


#================ Search Arxiv for relevant records ========================================
arxiv_api_url = "http://export.arxiv.org/api/query?search_query="

# Set up Arxiv parameters
start = 0
max_tries = CONFIG['arxiv_config']['max_tries']  # Set the maximum number of attempts to fetch data from PubMed in case of failure.
sleep_between_tries = CONFIG['arxiv_config']['time_sleep']  # Set the time to sleep between queries in seconds
retmax = CONFIG['arxiv_config']['retmax']  # Set the maximum number of records to retrieve from PubMed in a single query. Max is 30000.
results_per_iteration = CONFIG['arxiv_config']['results_per_iteration']  # Number of results to fetch per iteration


query = '%s&start=%i&max_results=%i' % (arxivq,
                                        start,
                                        results_per_iteration)

# This is a hack to expose both of these namespaces in
# feedparser v4.1
#feedparser._FeedParserMixin.namespaces['http://a9.com/-/spec/opensearch/1.1/'] = 'opensearch'
#feedparser._FeedParserMixin.namespaces['http://arxiv.org/schemas/atom'] = 'arxiv'

# perform a GET request using the base_url and query
response = urllib.request.urlopen(arxiv_api_url+query).read()

# parse the response using feedparser
feed = feedparser.parse(response)

# print out feed information
print('Feed title: %s' % feed.feed.title)
print('Feed last updated: %s' % feed.feed.updated)

# print opensearch metadata
print('totalResults for this query: %s' % feed.feed.opensearch_totalresults)
print('itemsPerPage for this query: %s' % feed.feed.opensearch_itemsperpage)
print('startIndex for this query: %s'   % feed.feed.opensearch_startindex)