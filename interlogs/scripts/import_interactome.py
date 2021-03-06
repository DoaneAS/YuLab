filename = "/Users/ashleysdoane/YuLab/interlogs/imports/InteractomePilot.txt"
dd=[]
with open(filename) as src:
    for line in src:
        #dd = [line.strip().split("\t")[1]],
        sym = line.split("\t")[2]
        dd.append(sym) #.upper())

iids = dd[1:]
i_ids = []
for i in iids:
    if i.startswith("Sc"):
        i_ids.append(i[2:])
    else:
        i_ids.append(i)
collab_ids = []
for i in i_ids:
    collab_ids.append(i.upper())

Cids = []
for k in ddd_all.keys():
    Cids.append(k[0])
    Cids.append(k[1])
cids = set(Cids)
cids = list(cids)

count = 0
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
    count = 0
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
                'accs': [],
                'gene': [],
                'refseq': []
            }
            entry = metadata_by_acc[uniprot_id]


        if tag == "AC":
            accs = [w.replace(";", "") for w in words]
            entry['accs'].extend(accs)

        if tag == "GN":
            if 'gene' not in entry:
                entry['gene'] = []
            parts = words[0].split("=")
            count += 1
            if "Name" in parts[0]:
                #print parts[1], count
                entry['gene'].append(parts[1].replace(";",""))
            else:
                entry['gene'].append(parts)
                print parts

        if tag == "DE":
            # if "RecName" in words[0]:
            if 'full_name' not in entry:
                fid = [w.replace("Full=", " ") for w in words[1:]]
            entry['full_name'] = [" ".join(fid[0:])]

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
    #return uni_txt
    res = parse_uniprot(uni_txt)
    return res


import sys
import pprint


ids = cids
res = pull_uni(ids)

uni2sym = {}
for k in res.keys():
    for x in res[k]['accs']:
        for v in res[k]['gene']: #res[k]['gene']):
            if not x in uni2sym:
                uni2sym[x] = v
            else:
                if uni2sym[x] != None:
                    entry = uni2sym[x]
                    uni2sym[x] = ''.join(entry)
uni= {}
for k in uni2sym.keys():
    kk = k.encode('ascii', 'ignore')
    key = kk
    val = uni2sym[k]
    uni[key] = val.encode('ascii', 'ignore')

inv_map = {v: k for k, v in uni.items()}

matches = {}
for c in collab_ids:
    if c in inv_map:
        matches[c] = inv_map.get(c)

cd = {z: w for w, z in matches.items()}

#all dict

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

#hiconf dict
dd_hc = {}
dd_hc = {k: v for k, v in ddd_all.iteritems() if len(v) >1}

#allI Int with candida genes

ddd_allAlb = {k: v for k, v in ddd_all.iteritems()}

for k, v in ddd_allAlb.items():
    a = k[0]
    b = k[1]
    indA = 0
    indB = 0
    if a in cd:
        c_a = cd.get(a)
        indA += 1
    if b in cd:
        c_b = cd.get(b)
        indB += 1
    if indA != 0:
        entry = ddd_allAlb[k]
        if not 'c_A' in entry:
            entry['c_A'] = []
        entry['c_A'].append(c_a)
    if indB == 1:
        entry = ddd_allAlb[k]
        if not 'c_B' in entry:
            entry['c_B'] = []
        entry['c_B'].append(c_b)


##add albicans symbols to dd_hc
# for k, v in dd_hc.items():
#     a = k[0]
#     b = k[1]
#     indA = 0
#     indB = 0
#     if a in cd:
#         c_a = cd.get(a)
#         indA += 1
#     if b in cd:
#         c_b = cd.get(b)
#         indB += 1
#     if indA != 0:
#         entry = dd_hc[k]
#         if  'c_A' not in entry:
#             entry['c_A'] = []
#         entry['c_A'].append(c_a)
#     if indB != 0:
#         entry = dd_hc[k]
#         if  'c_B' not in entry:
#             entry['c_B'] = []
#         entry['c_B'].append(c_b)




output = open('results/intAll.txt', 'w')
output.write('alb_interalog'+'\t'+'alb_A'+'\t' +'alb_B'+'\t'+'hu_int'+'\t'+'cerv_int'+'\t'+'pombe_int'+'\n')
for k in ddd_allAlb.keys():
    key = '%s-%s' % k
    output.write(key)

    if 'c_A' in ddd_allAlb[k]:
        output.write("\t")
        count = 0
        for v in ddd_allAlb[k]['c_A']:
            count += 1
            if count >1:
                output.write("; "+"-".join(v))
            else:
                output.write(v)
    else:
        output.write('\t')

    if 'c_B' in ddd_allAlb[k]:
        output.write("\t")
        count = 0
        for v in ddd_allAlb[k]['c_B']:
            count += 1
            if count >1:
                output.write("; "+"-".join(v))
            else:
                output.write(v)
    else:
        output.write('\t')
    if 'hu' in ddd_allAlb[k]:
        output.write('\t')
        count = 0
        for v in ddd_allAlb[k]['hu']:
            count += 1
            if count > 1:
                output.write("; " + "-".join(flatten(v)))
            else:
                output.write('-'.join(flatten(v)))
    else:
        output.write('\t' + 'none')

    if 'sc' in ddd_allAlb[k]:
        output.write('\t')
        count = 0
        for v in ddd_allAlb[k]['sc']:
            count += 1
            if count > 1:
                output.write("; " + "-".join(flatten(v)))
            else:
                output.write('-'.join(flatten(v)))

    else:
        output.write('\t' + 'none')

    if 'sp' in ddd_allAlb[k]:
        output.write('\t')
        count = 0
        for v in ddd_allAlb[k]['sp']:
            count += 1
            if count > 1:
                output.write("; " + "-".join(flatten(v)))
            else:
                output.write('-'.join(flatten(v)))
    else:
        output.write('\t' + 'none')

    output.write('\n')
output.close()



######
#hi conf with candA and candB##
# dd_hcAlb = {}

dd_hcAlb = {k: v for k, v in dd_hc.iteritems() if 'c_A' in v and 'c_B' in v}

output = open('results/candidaAandB_HiConf.txt', 'w')
output.write('alb_interolog' + '\t' + 'alb_A' + '\t' + 'alb_B' +
             '\t' + 'hu_int' + '\t' + 'cerv_int' + '\t' + 'pombe_int' + '\n')
for k in dd_hcAlb.keys():
    key = '%s-%s' % k
    output.write(key)

    if 'c_A' in dd_hcAlb[k]:
        output.write("\t")
        for v in dd_hcAlb[k]['c_A']:
            output.write(v)
    else:
        output.write('\t' + 'none')

    if 'c_B' in dd_hcAlb[k]:
        output.write("\t")
        count = 0
        for v in dd_hcAlb[k]['c_B']:
            count += 1
            if count > 1:
                output.write("; " + "-".join(v))
            else:
                output.write(v)
    else:
        output.write('\t' + 'none')

    if 'hu' in dd_hcAlb[k]:
        output.write('\t')
        count = 0
        for v in dd_hcAlb[k]['hu']:
            count += 1
            if count > 1:
                output.write("; " + "-".join(flatten(v)))
            else:
                output.write('-'.join(flatten(v)))
    else:
        output.write('\t' + 'none')

    if 'sc' in dd_hcAlb[k]:
        output.write('\t')
        count = 0
        for v in dd_hcAlb[k]['sc']:
            count += 1
            if count > 1:
                output.write("; " + "-".join(flatten(v)))
            else:
                output.write('-'.join(flatten(v)))

    else:
        output.write('\t' + 'none')

    if 'sp' in dd_hcAlb[k]:
        output.write('\t')
        count = 0
        for v in dd_hcAlb[k]['sp']:
            count += 1
            if count > 1:
                output.write("; " + "-".join(flatten(v)))
            else:
                output.write('-'.join(flatten(v)))
    else:
        output.write('\t' + 'none')

    output.write('\n')
output.close()
