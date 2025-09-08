# Title: Defining query for PubMed Metadata Extraction
# Last Updated: 2025-09-08
# Description: This script imports necessary libraries for extracting PubMed Metadata 
# Source: https://github.com/TLDWTutorials/PubmedAPI/blob/main/pubmed_api_in_python_2024.py

#================ QUERY PARAMETERS ========================================
# Define lists of authors
authors = ['Bryan Holland', 'Mehmet Oz', 'Anthony Fauci']  # Example authors, adjust as needed

# Define list of topics
topics = ['RNA', 'cardiovascular']  # Example topics, adjust as needed

# Define date range
date_start = '2012/03/01'
date_end = '2024/12/31'


#================ QUERY FORMATTING ========================================
# Build the query dynamically based on the available authors and topics.
def pubmed_query(authors, topics, date_start, date_end):
    '''
    Constructs a PubMed query string based on authors, topics, and date range.

    authors: List of author names. E.g., ['Author One', 'Author Two']
    topics: List of topics/keywords. E.g., ['cancer', 'diabetes']
    date_start: Start date in "YYYY/MM/DD" format. E.g., "2010/01/01"
    date_end: End date in "YYYY/MM/DD" format. E.g., "2020/12/31"

    Returns: A formatted PubMed query string.
    '''
    queries = []

    if authors:
        author_queries = ['{}[Author]'.format(author) for author in authors]
        queries.append('(' + ' OR '.join(author_queries) + ')')

    if topics:
        topic_queries = ['{}[Title/Abstract]'.format(topic) for topic in topics]
        queries.append('(' + ' OR '.join(topic_queries) + ')')

    if date_start and date_end:
        date_range = '("' + date_start + '"' + "[Date - Create]" + ' : "' + date_end + '"' + "[Date - Create]" + ')'

    full_query = ' AND '.join(queries) + ' AND ' + date_range

    return full_query


def arxiv_query(authors, topics, date_start, date_end):
    '''
    Constructs a PubMed query string based on authors, topics, and date range.

    authors: List of author names. E.g., ['Author One', 'Author Two']
    topics: List of topics/keywords. E.g., ['cancer', 'diabetes']
    date_start: Start date in "YYYY/MM/DD" format. E.g., "2010/01/01"
    date_end: End date in "YYYY/MM/DD" format. E.g., "2020/12/31"
    
    Returns: A formatted PubMed query string.
    '''
    queries = []

    if authors:
        author_queries = ['au:{}'.format(author) for author in authors]
        author_queries = map(lambda authors: authors.replace(' ', '_'), authors)
        queries.append('(' + ' OR '.join(author_queries) + ')')

    if topics:
        topic_queries = ['all:{}'.format(topic) for topic in topics]
        queries.append('(' + ' OR '.join(topic_queries) + ')')

    if date_start and date_end:
        date_start_arxiv = date_start.replace('/', '')
        date_end_arxiv = date_end.replace('/', '')
        date_range = 'submittedDate:' + '[' + date_start_arxiv + '+TO+' + date_end_arxiv + ']' 

    full_query = ' AND '.join(queries) + ' AND ' + date_range

    return full_query



pmq = pubmed_query(authors, topics, date_start, date_end)
arxivq = arxiv_query(authors, topics, date_start, date_end)

print(pmq)
print(arxivq)