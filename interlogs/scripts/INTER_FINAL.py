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
        #return get_pairs(file)
        p= get_Ps(file)
        pp = p[0]
        ppid = get_ppi_sp(pp)
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

def get_ppi_sp(src):
    ppi_dict ={}
    ppi = None
    for i in src:
        symbol = (i[2], i[3])
        ppi = tuple([i[0], i[1]])
        ppi_dict[ppi] = {
        "uniprot_sp": ppi,
        'symbol': symbol,
        'CA_ortholog': None,
        }
        #ppi[tuple([i[0], i[1]])] = val
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


    #import mygene
    #mg = mygene.MyGeneInfo()
    #out = mg.querymany(ids, spopes='entrezgene', fields='uniprot', species='human')
#outdf = mg.querymany(ids, spopes='entrezgene', fields='uniprot', species='human')





def get_sp_ids(intlist):
    sp_ids = []
    for i in intlist:
        a = i[0]
        b= i[1]
        sp_ids.append(a)
        sp_ids.append(b)
    return sp_ids


from bioservices.uniprot import UniProt
u = UniProt(verbose=False)

def get_uni_mapping(ids, frdb, todb):
    query = ' '.join(ids)
    queries= list(set(query.split()))
    res = u.mapping(fr=frdb, to= todb, query = queries)
    return res

def get_orth_matches_sp(intlist, sp2uni, sp2ca):
    ddres = {}
    for i in intlist:
        indA = 0
        indB = 0
        a = i[0]
        b= i[1]
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
                ddres[a,b] = {'A_sp_prot':a, 'B_sp_prot' :b,
                              'Auni' : sp2uni[a],
                              'Buni': sp2uni[b],
                              'A*':[],
                              'B*':[],
                            }
                entry = ddres[a,b]
                entry['A*'].append(una)
                entry['B*'].append(unb)
    return ddres

def sp_int_by_ca(res):
    ddnew = defaultdict(list)
    for k, v in res.items():
        pa = v['A_sp_prot']
        pb = v['B_sp_prot']
        pp = [pa, pb]
        A= flatten(v['A*'])
        B = flatten(v['B*'])
        for a in A:
            for b in B:
                if not (a, b) in ddnew:
                    ddnew[(a, b)] = [pp]
                else: ddnew[(a,b)].append(pp)
    return ddnew

def get_interalogs_sp(ppi_filename, inpar_filename):
    ppi_filename  = "/Users/ashleysdoane/YuLab/interlogs/spBinary_All.txt"
    p = get_seq_ppi_sp(ppi_filename)
    #gids = get_gids(p)
    intlist = get_intlist_sp(p)
    #g2uni = get_ids2unp(gids)
    #entrez2uni = get_entrez2uni(g2uni)
    inpar_sp = "//Users/ashleysdoane/YuLab/interlogs/C.albicans-S.pombe.txt"
    sp2ca = orth_parse_sp2ca(inpar_sp)
    frdb = "POMBASE_ID"
    todb="ACC"
    sp_ids = get_sp_ids(intlist)
    sp2uni = get_uni_mapping(sp_ids, frdb, todb)
    res = get_orth_matches_sp(intlist, sp2uni, sp2ca)
    #res = get_orth_matches_sp(intlist, sp2ca)
    return res


ppi_filename  = "/Users/ashleysdoane/YuLab/interlogs/spBinary_All.txt"
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
        output.write('\t'+'-'.join(v))
    output.write('\n')

ddd_all = {}
for k, v in ca_hu.items():
    if not k in ddd_all:
        ddd_all[k] = {"hu":v}
    else:
        entry = ddd_all[k]
        entry['hu'].append(v)
        #dd_all[k].append((v))
for i, j in ca_sc.items():
    if not i in ddd_all:
        ddd_all[i] = {'sc':j}
    else:
        entry = ddd_all[i]
        if not 'sc' in entry:
            entry['sc'] = []
        entry['sc'].append(j)
for m, n in ca_sp.items():
    if not m in ddd_all:
        ddd_all[m] = {'sp':n}
    else:
        entry = ddd_all[m]
        if not 'sp' in entry:
            entry['sp'] = []
        entry['sp'].append(n)

#write to file, dict above

output = open('results/ca_hu.sc.sp.txt', 'w')
output.write('ca'+'\t'+'hu'+'\t'+'sc'+'\t'+'sp'+'\n')
for k in ddd_all.keys():
    key = '%s %s' % k
    output.write(key)
    if 'hu' in ddd_all[k]:
        output.write("\t")
        count = 0
        for v in ddd_all[k]['hu']:
            count += 1
            if count >1:
                output.write("; "+"-".join(flatten(v)))
            else:
                output.write('-'.join(flatten(v)))
    else:
        output.write('\t')
    if 'sc' in ddd_all[k]:
        output.write("\t")
        count = 0
        for v in ddd_all[k]['sc']:
            count += 1
            if count >1:
                output.write("; "+"-".join(flatten(v)))
            else:
                output.write('-'.join(flatten(v)))
    else:
        output.write('\t')

    if 'sp' in ddd_all[k]:
        output.write("\t")
        count = 0
        for v in ddd_all[k]['sp']:
            count += 1
            if count >1:
                output.write("; "+"-".join(flatten(v)))
            else:
                output.write('-'.join(flatten(v)))
    output.write('\n')
output.close()