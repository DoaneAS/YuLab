from pandas import Series, DataFrame
import pandas as pd

filename= "HumanBinary_All.txt"

def get_items_from_file(filename):
    """Return all the items in the file named filename; if testfn
    then include only those items for which testfn is true"""
    with open(filename) as file:
        #return get_pairs(file)
        p= get_Ps(file)
        pp = p[0]
        print pp
        ppid = get_ppi(pp)
        return ppid


#def get_items

# def get_pairs(src):
#     pairs = [[line.split()[0], line.split()[1]] for line in src if line[0] != '']
#     for pp in pairs:
#         k = (pp[0], pp[1])
#     return k
def get_Ps(src):
    Ps = [[line.split()[0:5] for line in src]]
    return Ps

def get_pairs(src):
    pairs = ([line.split()[0], line.split()[1]] for line in src if line[0] != '')
    return pairs







def get_ppi(src):
    ppi_dict ={}
    ppi = None
    for i in src:
        symbol = (i[2], i[3])
        ppi = tuple([i[0], i[1]])
        ppi_dict[ppi] = {
        "gene_id": ppi,
        'symbol': symbol,
        }
        #ppi[tuple([i[0], i[1]])] = val
    return ppi_dict




p = get_items_from_file(filename)



import requests
import sys
import argparse

url = 'http://www.uniprot.org/'


def get_parts(uni_txt):
    for line in uni_txt.splitlines():
        k, v = [line[:5].strip(), line[5:].strip()]
        return k, v.split()


def parse_uniprot(uni_txt):
    """Parse metadata text from uniprot.org

     Returns a dictionary with the UNIPROT ACC as keys.
     """
    uniprot_id = None
    metadata_by_acc = {}
    tag = None
    # with get_parts
    for l in uni_txt.splitlines():
        ttag = l[:5].strip()
        if ttag and ttag != tag:
            tag = ttag
        words = (l[5:].strip()).split()
        if tag == "ID":
            uniprot_id = words[0]
            is_reviewed = words[1].startswith('Reviewed')
            length = int(words[2])
            metadata_by_acc[uniprot_id] = {
                'id': uniprot_id,
                'is_reviewed': is_reviewed,
                'length': length,
                'sequence': '',
                'accs': [],
                'gene': [],
                'refseq': [],
                'geneid':[]
            }
            entry = metadata_by_acc[uniprot_id]

        if tag == "SQ":
            if words[0] != "SEQUENCE":
                entry['sequence'] += ''.join(words)

        if tag == "AC":
            accs = [w.replace(";", "") for w in words]
            entry['accs'].extend(accs)

        if tag == "GN":
            if 'gene' not in entry:
                entry['gene'] = []
            parts = words[0].split("=")
            entry['gene'].append(parts[0:])

        if tag == "DE":
            # if "RecName" in words[0]:
            if 'full_name' not in entry:
                fid = [w.replace("Full=", " ") for w in words[1:]]
            entry['full_name'] = [" ".join(fid[0:])]

        if tag == "DR":
            if 'PDB' in words[0]:
                if 'pdb' not in entry:
                    entry['pdb'] = [words[1][:-1]]
                if 'pdbs' not in entry:
                    entry['pdbs'] = []
                entry['pdbs'].append(words[1][:-1])
            if 'RefSeq' in words[0]:
                if 'refseq' not in entry:
                    entry['refseq'] = []
                ids = [w[:-1] for w in words[1:]]
                entry['refseq'].extend(ids)
            if 'GeneID' in words[0]:
                if 'geneid' not in entry:
                    entry['geneid'] = []
                #ids = [w[:-1] for w in words[1:]]
                entry['geneid'].append(words[1][:-1])



    return metadata_by_acc


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
    res = parse_uniprot(uni_txt)
    resl = []
    resit = res.iteritems()
    for i in resit:
        resl.append(i)
    return resl


import sys
import pprint




filename= "CA_HU.txt"
#filename = "C.test.txt"

def ortho_parse(filename):
    src = open(filename)
    reslist = []
    Can_Hu = {}
    for line in open(filename):
        if line[24:30] == "100.00":
            key = line[0:6].strip()
            huval = line[35:43].strip()
            Can_Hu[key] = huval
    return Can_Hu

def ortho_parse_hu(filename):
    src = open(filename)
    reslist = []
    for line in open(filename):
        if line[24:30] == "100.00":
            huval = line[35:43].strip()
            reslist.append((huval))
    return reslist

def ortho_parse_dict(filename):
    src = open(filename)
    reslist = []
    Can_Hu = {}
    for line in open(filename):
        if line[24:30] == "100.00":
            ca = line[0:6].strip()
            hu = line[35:43].strip()
            Can_Hu[ca] = {
                "ca":ca,
                "hu": hu,
                "geneid":[]}
            entry = Can_Hu[ca]
    return Can_Hu


dd = ortho_parse(filename)
ddd = ortho_parse_dict(filename)
ids = dd.values()
res = pull_uni(ids)

def add_geneid(ddp, res):
    """ddp is dictionary of Ca:Hu
    res is uniprot annotations for Hu gene in ddp"""
    for l, d in res:
        for k in ddp.iterkeys():
            entry = ddp[k]
            if len((entry['hu']).intersection(d['accs'])) != 0:
                entry['geneid'].update(d['geneid'])
    return ddp
