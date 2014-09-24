filename = "/Users/ashleysdoane/YuLab/interlogs/ScBinary_All.txt"
def get_seq_ppi_sc(filename):
    """Return all the items in the file named filename; if testfn
    then include only those items for which testfn is true"""
    with open(filename) as file:
        #return get_pairs(file)
        p= get_Ps(file)
        pp = p[0]
        ppid = get_ppi_sc(pp)
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
        #ppi[tuple([i[0], i[1]])] = val
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


    #import mygene
    #mg = mygene.MyGeneInfo()
    #out = mg.querymany(ids, scopes='entrezgene', fields='uniprot', species='human')
#outdf = mg.querymany(ids, scopes='entrezgene', fields='uniprot', species='human')





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
    return newd

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
        hpp = [(v['Auni'], v['Buni'])]
        A= flatten(v['A*'])
        B = flatten(v['B*'])
        for a in A:
            for b in B:
                if not (a, b) in ddnew:
                    ddnew[(a, b)] = hpp
                else: ddnew[(a,b)].append(hpp)
    return ddnew


# def get_interalogs_sc(ppi_filename, inpar_filename):
#     ppi_filename  = "/Users/ashleysdoane/YuLab/interlogs/ScBinary_All.txt"
#     p = get_seq_ppi_sc(ppi_filename)
#     #gids = get_gids(p)
#     intlist = get_intlist_sc(p)
#     #g2uni = get_ids2unp(gids)
#     #entrez2uni = get_entrez2uni(g2uni)
#     inpar_sc = "//Users/ashleysdoane/YuLab/interlogs/C.albicans-S.cerevisiae.txt"
#     sc2ca = orth_parse_sc2ca(inpar_sc)
#     frdb = "ENSEMBLGENOME_PRO_ID"
#     todb="ACC"
#     sc_ids = get_sc_ids(intlist)
#     sc2uni = get_uni_mapping(sc_ids, frdb, todb)
#     res = get_orth_matches_sc(intlist, sc2uni, sc2ca)
#     #res = get_orth_matches_sc(intlist, sc2ca)
#     return res

ppi_filename  = "/Users/ashleysdoane/YuLab/interlogs/ScBinary_All.txt"
p = get_seq_ppi_sc(ppi_filename)
    #gids = get_gids(p)
intlist = get_intlist_sc(p)
    #g2uni = get_ids2unp(gids)
    #entrez2uni = get_entrez2uni(g2uni)
inpar_sc = "//Users/ashleysdoane/YuLab/interlogs/C.albicans-S.cerevisiae.txt"
sc2ca = orth_parse_sc2ca(inpar_sc)
frdb = "ENSEMBLGENOME_PRO_ID"
todb="ACC"
sc_ids = get_sc_ids(intlist)
sc2uni = get_uni_mapping(sc_ids, frdb, todb)
sc2ca_interalogs = get_orth_matches_sc(intlist, sc2uni, sc2ca)
ca2sc_int = sc_int_by_ca(sc2ca_interalogs)
ca_sc = dict(ca2sc_int)
ind = ca2sc_int.keys()
ca_sc_series = Series(ca2sc_int, index=ind)