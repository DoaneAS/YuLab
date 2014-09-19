from pandas import Series, DataFrame
import pandas as pd

#filename= "~/YuLab/interlogs/HumanBinary_All.txt"
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

#p = get_seq_ppi(filename)


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
    #out = mg.querymany(ids, scopes='entrezgene', fields='uniprot', species='human')
#outdf = mg.querymany(ids, scopes='entrezgene', fields='uniprot', species='human')

def get_ids2unp(ids):
    import mygene
    mg = mygene.MyGeneInfo()
    out = mg.querymany(ids, scopes='entrezgene', fields='uniprot', species='human')
    g2unip = {}
    for g in out:
        if 'uniprot' in g:
            v = g.get('uniprot')
            k = g.get('query')
        g2unip[k] = v
    return g2unip

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


def get_unikey(out):
    uni_key= {}
    for g in out:
        gid = g.get('query')
        if 'uniprot' in g:
            if 'Swiss-Prot' in g['uniprot']:
                v = g['uniprot']['Swiss-Prot']
                if type(v) == list:
                    for i in v:
                        uni_key[i] = {
                            'gid': gid}
                if type(v) != list:
                    uni_key[v] =  {'gid': gid}
    return uni_key


#filename = "C.test.txt"

def orth_parse_u(inpar_filename):
    src = open(filename)
    ddp = {}
    for line in src:
        if "%" in line[24:31]:
            ca = line[0:6].strip()
            ddp[ca] = {
                    "hu": set(),
        }
            entry = ddp[ca]
        if "%" in line[65:67]:
            entry['hu'].add(line[35:43].strip())
    return ddp




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
        if "%" in line[65:67]:
            hu = (line[35:43].strip())
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
                    una = uniprots
                    indA += 1
            if type(uniprots) == list:
                for u in uniprots:
                    if u in hu2ca:
                        una = uniprots
                        indA += 1

        if b in entrez2uni:
            uniprots = entrez2uni.get(b)
            if type(uniprots) != list:
                if uniprots in hu2ca:
                    unb = uniprots
                    indB += 1
            if type(uniprots) == list:
                for u in uniprots:
                    if u in hu2ca:
                        unb = uniprots
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

def get_interalogs(ppi_filename, inpar_filename):
    p = get_seq_ppi(ppi_filename)
    gids = get_gids(p)
    intlist = get_intlist(p)
    g2uni = get_ids2unp(gids)
    entrez2uni = get_entrez2uni(g2uni)
    hu2ca = orth_parse_hu2ca(inpar_filename)
    res = get_orth_matches(intlist, entrez2uni, hu2ca)
    return res





if __name__ == '__main__':
    ppi_filename= "/Users/ashleysdoane/YuLab/interlogs/HumanBinary_All.txt"
    inpar_filename= "CA_HU.txt"
    res = get_interalogs(ppi_filename, inpar_filename)
    import pprint
    pprint.pprint(res)

