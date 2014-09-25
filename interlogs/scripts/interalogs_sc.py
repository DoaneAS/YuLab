from pandas import Series, DataFrame
import pandas as pd

from bioservices.uniprot import UniProt
u = UniProt(verbose=False)
from collections import defaultdict

def get_seq_ppi_sc(filename):
    """Return all the items in the file named filename; if testfn
    then include only those items for which testfn is true"""
    with open(filename) as file:
        #return get_pairs(file)
        p= get_Ps(file)
        pp = p[0]
        ppid = get_ppi_sc(pp)
        return ppid

def get_Ps(src):
    Ps = [[line.split()[0:5] for line in src]]
    return Ps

def get_pairs(src):
    pairs = ([line.split()[0], line.split()[1]] for line in src if line[0] != '')
    return pairs

def get_ppi_sc(src):
    ppi_dict ={}
    ppi = None
    for i in src:
        symbol = (i[2], i[3])
        ppi = tuple([i[0], i[1]])
        ppi_dict[ppi] = {
        "uniprot_sc": ppi,
        'symbol': symbol,
        'CA_ortholog': None,
        }
    return ppi_dict

def get_intlist_sc(intdata):
    intlist = []
    for k in intdata.keys():
        intlist.append(k)
    return intlist

def orth_parse_sc2ca(filename):
    ca_list = {}
    sc_list = {}
    src = open(filename)
    ddp = {}
    for line in src:
        if "%" in line[24:31]:
            ca = line[0:6].strip()
        if "%" in line[54:58]:
            sc = line[28:36].strip()
            ddp[sc] = ca
    return ddp

def get_gids(ppi):
    geneids = ppi.keys()
    ids = []
    for i, j in geneids:
        ids.append(i)
        ids.append(j)
    ids = list(set(ids))
    return ids

def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

def get_sc_ids(intlist):
    sc_ids = []
    for i in intlist:
        a = i[0]
        b= i[1]
        sc_ids.append(a)
        sc_ids.append(b)
    return sc_ids

from bioservices.uniprot import UniProt
u = UniProt(verbose=False)

def get_uni_mapping_old(ids, frdb, todb):
    query = ' '.join(ids)
    queries= list(set(query.split()))
    res = u.mapping(fr=frdb, to= todb, query = queries)
    return res

def get_uni_mapping(ids, frdb, todb):
    query = ' '.join(ids)
    queries= list(set(query.split()))
    res = u.mapping(fr=frdb, to= todb, query = queries)
    newd = defaultdict(list)
    for k, v in res.items():
        kk = k.encode('ascii', 'ignore')
        vv = []
        for i in v:
            vv.append(i.encode('ascii', 'ignore'))
        newd[kk] = vv
    return res

def get_orth_matches_sc(intlist, sc2uni, sc2ca):
    ddres = {}
    for i in intlist:
        indA = 0
        indB = 0
        a = i[0]
        b= i[1]
        if a in sc2uni:
            uniprots = sc2uni.get(a)
            if type(uniprots) != list:
                if uniprots in sc2ca:
                    una = sc2ca[uniprots]
                    indA += 1
            if type(uniprots) == list:
                for u in uniprots:
                    if u in sc2ca:
                        una = sc2ca[u]
                        indA += 1

        if b in sc2uni:
            uniprots = sc2uni.get(b)
            if type(uniprots) != list:
                if uniprots in sc2ca:
                    unb = sc2ca[uniprots]
                    indB += 1
            if type(uniprots) == list:
                for w in uniprots:
                    if w in sc2ca:
                        unb = sc2ca[w]
                        indB += 1
        if indA + indB >= 2:
                ddres[a,b] = {'A_sc_prot':a, 'B_sc_prot' :b,
                              'Auni' : sc2uni[a],
                              'Buni': sc2uni[b],
                              'A*':[],
                              'B*':[],
                            }
                entry = ddres[a,b]
                entry['A*'].append(una)
                entry['B*'].append(unb)
    return ddres

def sc_int_by_ca(res):
    ddnew = defaultdict(list)
    for k, v in res.items():
        pa = v['A_sc_prot']
        pb = v['B_sc_prot']
        pp = [pa, pb]
        A= flatten(v['A*'])
        B = flatten(v['B*'])
        for a in A:
            for b in B:
                if not (a, b) in ddnew:
                    ddnew[(a, b)] = [pp]
                else: ddnew[(a,b)].append(pp)
    return ddnew

ppi_filename  = "/Users/ashleysdoane/YuLab/interlogs/ScBinary_All.txt"
p = get_seq_ppi_sc(ppi_filename)
intlist = get_intlist_sc(p)
inpar_sc = "//Users/ashleysdoane/YuLab/interlogs/C.albicans-S.cerevisiae.txt"
sc2ca = orth_parse_sc2ca(inpar_sc)
frdb = "ENSEMBLGENOME_PRO_ID"
todb="ACC"
sc_ids = get_sc_ids(intlist)
sc2uni = get_uni_mapping(sc_ids, frdb, todb)
sc2ca_interalogs = get_orth_matches_sc(intlist, sc2uni, sc2ca)
ca2sc_int = sc_int_by_ca(sc2ca_interalogs)
ca_sc = dict(ca2sc_int)

output = open('results/ca_sc.txt', 'w')
for k in ca_sc.keys():
    key = '%s %s' % k
    output.write(key)
    for v in ca_sc[k]:
        output.write('\t'+'-'.join(v))
    output.write('\n')
