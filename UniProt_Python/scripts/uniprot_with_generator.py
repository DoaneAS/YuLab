import requests
import sys
import argparse

url = 'http://www.uniprot.org/'
tag = None
props = {}

def pull_uni(ids, format='txt'):
    """ request entries by uniprot acc using batch retrieval
    Args:
        ids: list of ids to retrieve
        format: txt by default
        possible formats:
        txt, xml, rdf, fasta, gff"""

    tool = 'batch/'
    if type(ids) is not list:
        ids = [ids]
    query = ' '.join(ids)
    query = list(set(query.split()))
    queries = [query[i:i + 100] for i in xrange(0, len(query), 100)]
    data = {'format': format}
    responses = [
        requests.post(url + tool, data=data, files={'file': ' '.join(query)}) for query in queries]
    uni_txt = ''.join([response.text for response in responses])
    # return uni_txt
    res = get_items(uni_txt)
    return res

def testfn(l):
    l != None

def get_items(uni_txt, testfn=None):
    """Return all the items in src; if testfn then include only those
    items for which testfn is true"""
    return [line for line in item_gen(uni_txt)
        if not testfn or testfn(line)]

def item_gen(uni_txt, tag=None):
    tag = None
    uniprot_id = None
    #props = {}
    line = uni_txt.splitlines()
    while line:
        ttag = line[:5].strip()
        if ttag and ttag != tag:
            tag = ttag
        #words = (line[5:].strip()).split()
        #uniprot_id = words[0]
        # need to keep line to feed back to read_words
        #yield line
        #yield [line[:5].strip(), line[5:].strip()]
            feature = read_words(uni_txt, line, tag)
        # need to keep line to feed back to read_words
            yield feature

def is_feature_start(line):
    return line and line[5] != ' '

def read_words(uni_txt, line, tag):
    #uniprot_id = None
    #props = {}
    words = (line[5:].strip()).split()

    if tag == "ID":
            uniprot_id = words[0]
            is_reviewed = words[1].startswith('Reviewed')
            length = int(words[2])
            props[uniprot_id] = {
                'id': uniprot_id,
                'is_reviewed': is_reviewed,
                'length': length,
                #'sequence': '',
                #'accs': [],
                #'gene': [],
                #'refseq': []
            }


            #props[uniprot_id]['accs'].extend(accs)


            return props

if __name__ == '__main__':
    ids = ['P01116', 'P15056']
    res = pull_uni(ids)
    #print res
    import pprint
    pprint.pprint(res)