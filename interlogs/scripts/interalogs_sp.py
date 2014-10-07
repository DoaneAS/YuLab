from pandas import Series, DataFrame
import pandas as pd

from bioservices.uniprot import UniProt
u = UniProt(verbose=False)
from collections import defaultdict

filename = "/Users/ashleysdoane/YuLab/interlogs/SpBinary_All.txt"


def get_seq_ppi_sp(filename):
    """Return all the items in the file named filename; if testfn
    then include only those items for which testfn is true"""
    with open(filename) as file:
        # return get_pairs(file)
        p = get_Ps(file)
        pp = p[0]
        ppid = get_ppi_sp(pp)
        return ppid


def get_Ps(src):
    Ps = [[line.split()[0:5] for line in src]]
    return Ps


def get_pairs(src):
    pairs = ([line.split()[0], line.split()[1]] for line in src if line[0] != '')
    return pairs


def get_ppi_sp(src):
    ppi_dict = {}
    ppi = None
    for i in src:
        symbol = (i[2], i[3])
        ppi = tuple([i[0], i[1]])
        ppi_dict[ppi] = {
            "uniprot_sp": ppi,
            'symbol': symbol,
            'CA_ortholog': None,
        }
    return ppi_dict


def get_intlist_sp(intdata):
    intlist = []
    for k in intdata.keys():
        intlist.append(k)
    return intlist


def orth_parse_sp2ca(filename):
    ca_list = {}
    sp_list = {}
    src = open(filename)
    ddp = {}
    for line in src:
        if "%" in line[24:31]:
            ca = line[0:6].strip()
        if "%" in line[54:58]:
            sp = line[28:36].strip()
            ddp[sp] = ca
    return ddp


def get_gids(ppi):
    geneids = ppi.keys()
    ids = []
    for i, j in geneids:
        ids.append(i)
        ids.append(j)
    ids = list(set(ids))
    return ids


def get_sp_ids(intlist):
    sp_ids = []
    for i in intlist:
        a = i[0]
        b = i[1]
        sp_ids.append(a)
        sp_ids.append(b)
    return sp_ids


from bioservices.uniprot import UniProt
u = UniProt(verbose=False)


def get_uni_mapping(ids, frdb, todb):
    query = ' '.join(ids)
    queries = list(set(query.split()))
    res = u.mapping(fr=frdb, to=todb, query=queries)
    return res


def get_orth_matches_sp(intlist, sp2uni, sp2ca):
    ddres = {}
    for i in intlist:
        indA = 0
        indB = 0
        a = i[0]
        b = i[1]
        if a in sp2uni:
            uniprots = sp2uni.get(a)
            if type(uniprots) != list:
                if uniprots in sp2ca:
                    una = sp2ca[uniprots]
                    indA += 1
            if type(uniprots) == list:
                for u in uniprots:
                    if u in sp2ca:
                        una = sp2ca[u]
                        indA += 1

        if b in sp2uni:
            uniprots = sp2uni.get(b)
            if type(uniprots) != list:
                if uniprots in sp2ca:
                    unb = sp2ca[uniprots]
                    indB += 1
            if type(uniprots) == list:
                for w in uniprots:
                    if w in sp2ca:
                        unb = sp2ca[w]
                        indB += 1
        if indA + indB >= 2:
            ddres[a, b] = {'A_sp_prot': a, 'B_sp_prot': b,
                           'Auni': sp2uni[a],
                           'Buni': sp2uni[b],
                           'A*': [],
                           'B*': [],
                           }
            entry = ddres[a, b]
            entry['A*'].append(una)
            entry['B*'].append(unb)
    return ddres


def sp_int_by_ca(res):
    ddnew = defaultdict(list)
    for k, v in res.items():
        pa = v['A_sp_prot']
        pb = v['B_sp_prot']
        pp = [pa, pb]
        A = flatten(v['A*'])
        B = flatten(v['B*'])
        for a in A:
            for b in B:
                if not (a, b) in ddnew:
                    ddnew[(a, b)] = [pp]
                else:
                    ddnew[(a, b)].append(pp)
    return ddnew


def get_interalogs_sp(ppi_filename, inpar_filename):
    ppi_filename = "/Users/ashleysdoane/YuLab/interlogs/spBinary_All.txt"
    p = get_seq_ppi_sp(ppi_filename)
    intlist = get_intlist_sp(p)
    inpar_sp = "//Users/ashleysdoane/YuLab/interlogs/C.albicans-S.pombe.txt"
    sp2ca = orth_parse_sp2ca(inpar_sp)
    frdb = "POMBASE_ID"
    todb = "ACC"
    sp_ids = get_sp_ids(intlist)
    sp2uni = get_uni_mapping(sp_ids, frdb, todb)
    res = get_orth_matches_sp(intlist, sp2uni, sp2ca)
    return res


ppi_filename = "/Users/ashleysdoane/YuLab/interlogs/spBinary_All.txt"
p = get_seq_ppi_sp(ppi_filename)
inpar_sp = "//Users/ashleysdoane/YuLab/interlogs/C.albicans-S.pombe.txt"
sp2ca = orth_parse_sp2ca(inpar_sp)
sp_interalogs = get_interalogs_sp(ppi_filename, inpar_sp)
ca_sp_int = sp_int_by_ca(sp_interalogs)
ca_sp = dict(ca_sp_int)

output = open('results/ca_sp.txt', 'w')
for k in ca_sp.keys():
    key = '%s %s' % k
    output.write(key)
    for v in ca_sp[k]:
        output.write('\t' + '-'.join(v))
    output.write('\n')
