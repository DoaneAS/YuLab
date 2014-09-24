from pandas import Series, DataFrame
import pandas as pd

#filename= "~/YuLab/interlogs/HumanBinary_All.txt"
#filename= "HumanBinary_All.txt"

def get_seq_ppi(filename):
    """Return all the items in the file named filename; if testfn
    then include only those items for which testfn is true"""
    with open(filename) as file:
        #return get_pairs(file)
        p= get_Ps(file)
        pp = p[0]
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
        'CA_ortholog': None,
        }
        #ppi[tuple([i[0], i[1]])] = val
    return ppi_dict








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




#filename= "CA_HU.txt"
#filename = "C.test.txt"

def orth_parse(filename):
    src = open(filename)
    ddp = {}
    for line in src:
        if "%" in line[24:31]:
            ca = line[0:6].strip()
            ddp[ca] = {
                    "ca":ca,
                    "hu": set(),
                    "geneid":set()
        }
            entry = ddp[ca]
        if "%" in line[65:67]:
            entry['hu'].add(line[35:43].strip())
    return ddp


def orthl_parse(filename):
    src = open(filename)
    ddp = {}
    for line in src:
        if "%" in line[24:31]:
            ca = line[0:6].strip()
            ddp[ca] = {
                    "ca":ca,
                    "hu": [],
                    "geneid":set()
        }
            entry = ddp[ca]
        if "%" in line[65:67]:
            entry['hu'].append(line[35:43].strip())
    return ddp

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

def add_geneid(ddp, res):
    """cahu is dictionary of Ca:Hu
    res is uniprot annotations for Hu gene in ddp"""
    for l, d in res:
        for k in ddp.iterkeys():
            entry = ddp[k]
            if len((entry['hu']).intersection(d['accs'])) != 0:
                if len(entry['geneid']) == 0:
                    entry['geneid'] = set()
                entry['geneid'].update(d['geneid'])
    return ddp





#filename2= "CA_HU.txt"

def uni_ann(ddp):
    ids = []
    '''ddp is ca:hu dict'''
    for k in ddp.keys():
        entry = ddp[k]
        for i in list(entry['hu']):
            ids.append(i)
    res = pull_uni(ids)
    d = add_geneid(dd, res)
    return d

def get_hu2ca_dict(ddp):
    hu_ca_dict = {}
    for k in ddp.keys():
        entry = ddp[k]
        for eid in entry['geneid']:
            hu_ca_dict[eid] = entry['ca']
    return hu_ca_dict

def add_orthologs(p, hu2ca):
    '''p is human ppi data,
    hurca is human entrex id to CA ortholog mapping'''
    for i, j in p.keys():
        entry = p[i,j]
        A = hu2ca.get(i)
        B = hu2ca.get(j)
        if (A != None and B != None):
            entry['CA_ortholog'] = (A, B)

def make_interalogs(p, ddp):
    '''p is hum interaction data,
    ddp is ca to hu dict with ann'''
    interalogs = {}
    for k in p.keys():
        entry = p[k]
        if entry['CA_ortholog'] != None:
            interalogs[(entry['CA_ortholog'])] = p[k]
    for x, y in interalogs.keys():
        entry = interalogs[x,y]
        entry['A,B'] = (ddp[x]['hu'], ddp[y]['hu'])
    return interalogs

def make_interalogs2(p, ddp):
    '''p is hum interaction data,
    ddp is ca to hu dict with ann'''
    inter2 = {}
    for k in p.keys():
        entry = p[k]
        if entry['CA_ortholog'] != None:
            inter2[(entry['CA_ortholog'])] = p[k]
    for x, y in inter2.keys():
        entry = inter2[x,y]
        entry['A'] = ddp[x]['hu']
        entry['B'] = ddp[y]['hu']
    return inter2


if __name__ == '__main__':
    filename= "/Users/ashleysdoane/YuLab/interlogs/HumanBinary_All.txt"
    p = get_seq_ppi(filename)
    print p
    filename2= "/Users/ashleysdoane/YuLab/interlogs/CA_HU.txt"
    dd = orth_parse(filename2)
    ddp = uni_ann(dd)
    print ddp
    F = open('CA_Hu_dict2.pkl', 'wb')
    import pickle
    pickle.dump(ddp, F)
    F.close()
    hu2ca = get_hu2ca_dict(ddp)
    add_orthologs(p, hu2ca)
    interalogs = make_interalogs2(p, ddp)
    import pprint
    pprint.pprint(interalogs)

    F = open('interalogsCA_HU2.pkl', 'wb')
    pickle.dump(interalogs, F)
    F.close()