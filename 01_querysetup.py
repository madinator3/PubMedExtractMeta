# Title: Defining query for PubMed Metadata Extraction
# Last Updated: 2025-07-23
# Description: This script imports necessary libraries for extracting PubMed Metadata 

# Define lists of authors and topics
authors = ['Bryan Holland', 'Mehmet Oz', 'Anthony Fauci']  # Example authors, adjust as needed
topics = ['RNA', 'cardiovascular']  # Example topics, adjust as needed

# Define date range
date_range = '("2012/03/01"[Date - Create] : "2022/12/31"[Date - Create])'

# Build the query dynamically based on the available authors and topics
queries = []

if authors:
    author_queries = ['{}[Author]'.format(author) for author in authors]
    queries.append('(' + ' OR '.join(author_queries) + ')')

if topics:
    topic_queries = ['{}[Title/Abstract]'.format(topic) for topic in topics]
    queries.append('(' + ' OR '.join(topic_queries) + ')')

full_query = ' AND '.join(queries) + ' AND ' + date_range

