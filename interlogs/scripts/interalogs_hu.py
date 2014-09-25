from pandas import Series, DataFrame
import pandas as pd

from bioservices.uniprot import UniProt
u = UniProt(verbose=False)
from collections import defaultdict
filename= "HumanBinary_All.txt"

def get_seq_ppi(filename):
    """Return all the items in the file named filename; if testfn
    then include only those items for which testfn is true"""
    with open(filename) as file:
        #return get_pairs(file)
        p= get_Ps(file)
        pp = p[0]
        ppid = get_ppi(pp)
        return ppid

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
    return ppi_dict

def get_gids(ppi):
    geneids = ppi.keys()
    ids = []
    for i, j in geneids:
        ids.append(i)
        ids.append(j)
    ids = list(set(ids))
    return ids

#g2uni = get_ids2unp(out)
def get_ddd(p, g2uni):
    ddd={}
    for i, j in p.keys():
        ddd[(i,j)] = {
            'A': set(),
            'B' : set()
        }
        entry = ddd[(i,j)]
        if i in g2uni:
            if j in g2uni:
                entry['A'] = g2uni[i].values()[0]
                entry['B'] = g2uni[j].values()[0]
    dddset = {}
    for k, v in ddd.items():
        dddset[k] = v['A'], v['B']
    return dddset
    #ddd[(i,j)] = (A, B)


def get_intlist(intdata):
    intlist = []
    for k in intdata.keys():
        intlist.append(k)
    return intlist

def get_entrez2uni(g2uni):
    entrez2uni= {}
    for k,v in g2uni.items():
        entrez2uni[k] = g2uni[k].values()[0]
    return entrez2uni

def orth_parse_hu2ca(filename):
    ca_list = {}
    hu_list = {}
    src = open(filename)
    ddp = {}
    for line in src:
        if "%" in line[24:31]:
            ca = line[0:6].strip()
        if "%" in line[53:58]:
            hu = (line[28:36].strip())
            ddp[hu] = ca
    return ddp

def get_orth_matches(intlist, entrez2uni, hu2ca):
    ddres = {}
    for i in intlist:
        indA = 0
        indB = 0
        a = i[0]
        b= i[1]
        if a in entrez2uni:
            uniprots = entrez2uni.get(a)
            if type(uniprots) != list:
                if uniprots in hu2ca:
                    una = hu2ca[uniprots]
                    indA += 1
            if type(uniprots) == list:
                for u in uniprots:
                    if u in hu2ca:
                        una = hu2ca[u]
                        indA += 1

        if b in entrez2uni:
            uniprots = entrez2uni.get(b)
            if type(uniprots) != list:
                if uniprots in hu2ca:
                    unb = hu2ca[uniprots]
                    indB += 1
            if type(uniprots) == list:
                for w in uniprots:
                    if w in hu2ca:
                        unb = hu2ca[w]
                        indB += 1
        if indA + indB >= 2:
                ddres[a,b] = {'Agene':a, 'Bgene' :b,
                              'Auni' : entrez2uni[a],
                              'Buni': entrez2uni[b],
                              'A*':[],
                              'B*':[],
                            }
                entry = ddres[a,b]
                entry['A*'].append(una)
                entry['B*'].append(unb)
    return ddres



def get_hu_ids(intlist):
    hu_ids = []
    for i in intlist:
        a = i[0]
        b= i[1]
        hu_ids.append(a)
        hu_ids.append(b)
    return hu_ids

def get_uni_mapping(ids, frdb, todb):
    query = ' '.join(ids)
    queries= list(set(query.split()))
    res = u.mapping(fr=frdb, to= todb, query = queries)
    return res


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

def int_by_ca(res):
    ddnew = defaultdict(list)
    for k, v in res.items():
        hpa = v['Agene']
        hpb = v['Bgene']
        A= flatten(v['A*'])
        B = flatten(v['B*'])
        hpp = [hpa, hpb]
        for a in A:
            for b in B:
                if not (a, b) in ddnew:
                    ddnew[(a, b)] = [hpp]
                else: ddnew[(a,b)].append(hpp)
    return ddnew


ppi_filename= "/Users/ashleysdoane/YuLab/interlogs/HumanBinary_All.txt"
p = get_seq_ppi(ppi_filename)
intlist = get_intlist(p)
ids = get_hu_ids(intlist)
frdb = "P_ENTREZGENEID"
todb = "ACC"
entrez2uni = get_uni_mapping(ids, frdb,todb)
inpar_filename = "C.albicans-H.sapiens.txt"
hu2ca = orth_parse_hu2ca("C.albicans-H.sapiens.txt")
res = get_orth_matches(intlist, entrez2uni, hu2ca)
ca_res = int_by_ca(res)
ca_hu = dict(ca_res)

output = open('results/ca_hu.txt', 'w')
for k in ca_hu.keys():
    key = '%s %s' % k
    output.write(key)
    for v in ca_hu[k]:
        output.write('\t'+'-'.join(v))
    output.write('\n')

