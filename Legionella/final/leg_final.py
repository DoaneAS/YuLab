
# import legionella genome with biopython
# def  genes as leg_genes dict

import os
from Bio import SeqIO
from Bio import Entrez

filename = "Legionella.gbk"
record = SeqIO.read(filename, "genbank")

# parse seq format to dict
# extract selected legionella genes from genome
predicted_genes = []
leg_genes = {}
for i, feature in enumerate(record.features):
    if feature.type == 'CDS':
        product = feature.qualifiers['product'][0].lower()
        seq = feature.extract(record.seq)
        locus = feature.qualifiers['locus_tag'][0]
        predicted_genes.append(feature)
        trans = feature.qualifiers['translation']
        leg_genes[locus] = {
            "locus": locus, "seq": seq, "product": product, "translation": trans}


# Process blast result
# original selected leg genes
###


filename = "substrates.txt"
with open(filename) as src:
    p = [line.split('\r') for line in src][0]
legdict = {}
for line in p:
    pp = line.split()
    lid = pp[0].lower(),
    alias = pp[1]
    legdict[lid[0]] = {
        "lid": lid[0],
        "alias": alias}

from Bio import SearchIO
result_handle = open("results/blastp_resultAD.txt")
# comments=True)
idx = SearchIO.index('results/blastp_resultAD.txt', 'blast-tab')
idss = (idx.keys())
homs = {}
rec_list_all = []

for ids in idss:
    rec_list = []
    for rec in idx[ids].hsps:
        rec_list.append((rec.hit_id, rec.evalue, rec.ident_pct))
        rec_list_all.append(
            (rec.query_id, rec.hit_id, rec.evalue, rec.ident_pct))
    homs[ids] = {ids: rec_list}
leg_uniq = []
ids = []
for p in sorted(rec_list_all, key=lambda x: x[0][1]):
    if p[2] < 10E-3:
        if p[1] not in p[0]:
            if p[1] not in ids:
                ids.append(p[1])
                leg_uniq.append(p)
ids_nomatch = []
for k in legdict.keys():
    if k not in ids:
        ids_nomatch.append(k)

# some ids have more than one homologe and therefore are written more than once
output = open('results/leg_blastp_all.txt', 'w')
output.write(
    'id' + '\t' + 'blastp_hit' + '\t' + 'e-value' + '\t' + 'ident_pct' + '\n')
for p in leg_uniq:
    if p[0] in legdict.keys():  # check to make sure
        qid = '%s' % p[0]
        hit = '%s' % p[1]
        ev = '%s' % p[2]
        per = '%s' % p[3]
        output.write(qid + "\t" + hit + "\t" + ev + '\t' + per + "\n")
for q in ids_nomatch:
    output.write(q + "\n")
output.close()


ids_homs = []
for p in leg_uniq:
    ids_homs.append(p[1])
ids_orig = legdict.keys()

ids_all = list(set(ids_orig + ids_homs))
output = open('results/ids_all_unique.txt', 'w')
for h in ids_all:
    output.write(h + '\n')
output.close()

####
####
####
####

###
# import all leg genes in biopython
# extract genes in ids_all
# define new dict with DNA sequence for each gene from ids_all
###

import os
from Bio import SeqIO
from Bio import Entrez


input_handle = open("Legionella_pneumophila_Phil/CP003885.gbk")
record = SeqIO.read(input_handle, "genbank")

####
filename = "Legionella.gbk"
record = SeqIO.read(filename, "genbank")
predicted_genes = []
leg_genes = {}
for i, feature in enumerate(record.features):
    if feature.type == 'CDS':
        product = feature.qualifiers['product'][0].lower()
        seq = feature.extract(record.seq)
        locus = feature.qualifiers['locus_tag'][0]
        predicted_genes.append(feature)
        trans = feature.qualifiers['translation']
        leg_genes[locus] = {
            "locus": locus, "seq": seq, "product": product, "translation": trans}
leg_genes
# now, can use sequence info from this
# below, create seq object with selected leg genes for PCR
filename = "Legionella.gbk"
record = SeqIO.read(filename, "genbank")
predicted_genes = []
leg_genes = {}
for i, feature in enumerate(record.features):
    if feature.type == 'CDS':
        product = feature.qualifiers['product'][0].lower()
        seq = feature.extract(record.seq)
        locus = feature.qualifiers['locus_tag'][0]
        predicted_genes.append(feature)
        trans = feature.qualifiers['translation']
        leg_genes[locus] = {
            "locus": locus, "seq": seq, "product": product, "translation": trans}

seq_list = []
for g in ids_all:
    if g in leg_genes:
        entry = leg_genes[g]
        seq_list.append(entry['seq'])
# print str(seq_list[1])

seq_dict = {}
for g in ids_all:
    if g in leg_genes:
        entry = leg_genes[g]
        seq_dict[g] =  entry['seq']



from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import primer3

test_seq= str(seq_dict.get(ids_all[1]))
#can iterate through with a list of ids

ftail53 = Seq("GGGGACAACTTTGTACAAAAAAGTTGGCACC", IUPAC.unambiguous_dna)
rtail53 = Seq("GGGGACAACTTTGTACAAGAAAGTTGGCAA", IUPAC.unambiguous_dna)

seq_arg = {
    'SEQUENCE_ID': 'test_leg',
    'SEQUENCE_TEMPLATE': test_seq,
    # 'SEQUENCE_INCLUDED_REGION': [3, len(test_seq)],
    # 'SEQUENCE_PRIMER': str(fp),
    # 'SEQUENCE_PRIMER_REVCOMP': str(rp),
    'SEQUENCE_FORCE_LEFT_START': 3,
    'SEQUENCE_FORCE_RIGHT_START': len(test_seq)-1,
}

global_arg = {
    # 'PRIMER_OPT_SIZE': 51,
    'PRIMER_PICK_ANYWAY': 1,
    'PRIMER_PICK_LEFT_PRIMER': 1,
    'PRIMER_PICK_RIGHT_PRIMER': 1,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    # 'PRIMER_INTERNAL_MAX_SELF_END': 8,
    'PRIMER_MIN_SIZE': 16,
    'PRIMER_MAX_SIZE': 28,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 55.0,
    'PRIMER_MAX_TM': 68.0,
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    # 'PRIMER_MAX_POLY_X': 100,
    # 'PRIMER_INTERNAL_MAX_POLY_X': 100,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    # 'PRIMER_MAX_NS_ACCEPTED': 0,
    # 'PRIMER_MAX_SELF_ANY': 12,
    # 'PRIMER_MAX_SELF_END': 8,
    # 'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    # 'PRIMER_PAIR_MAX_COMPL_END': 8,
    # 'PRIMER_PRODUCT_SIZE_RANGE': [[75, 100], [100, 125], [125, 150], [150, 175], [175, 200], [200, 225]],
}

res = primer3.designPrimers(seq_arg, global_arg)
res = primer3.designPrimers(seq_arg, global_arg)

import pprint
# pprint.pprint(res)
# print res['PRIMER_LEFT_0_SEQUENCE']


fp = ftail53 + Seq(res.get('PRIMER_LEFT_0_SEQUENCE'), IUPAC.unambiguous_dna)
rp = rtail53 + Seq(res.get('PRIMER_RIGHT_0_SEQUENCE'), IUPAC.unambiguous_dna)
print fp, rp

primer3.calcTm(str(fp)), primer3.calcTm(str(rp))

# output = open('results/primers.txt', 'w')
# for k, v in ...
#    output.write(h + '\n')
# output.close()

seq_arg_batch = {
    'SEQUENCE_ID': 'test_leg',
    'SEQUENCE_TEMPLATE': test_seq,
    # 'SEQUENCE_INCLUDED_REGION': [3, len(test_seq)],
    # 'SEQUENCE_PRIMER': str(fp),
    # 'SEQUENCE_PRIMER_REVCOMP': str(rp),
    'SEQUENCE_FORCE_LEFT_START': 3,
    'SEQUENCE_FORCE_RIGHT_START': len(test_seq)-1,
}

global_arg_batch = {
    # 'PRIMER_OPT_SIZE': 51,
    'PRIMER_PICK_ANYWAY': 1,
    'PRIMER_PICK_LEFT_PRIMER': 1,
    'PRIMER_PICK_RIGHT_PRIMER': 1,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    # 'PRIMER_INTERNAL_MAX_SELF_END': 8,
    'PRIMER_MIN_SIZE': 16,
    'PRIMER_MAX_SIZE': 28,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 55.0,
    'PRIMER_MAX_TM': 68.0,
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    # 'PRIMER_MAX_POLY_X': 100,
    # 'PRIMER_INTERNAL_MAX_POLY_X': 100,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    # 'PRIMER_MAX_NS_ACCEPTED': 0,
    # 'PRIMER_MAX_SELF_ANY': 12,
    # 'PRIMER_MAX_SELF_END': 8,
    # 'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    # 'PRIMER_PAIR_MAX_COMPL_END': 8,
    # 'PRIMER_PRODUCT_SIZE_RANGE': [[75, 100], [100, 125], [125, 150], [150, 175], [175, 200], [200, 225]],
}



ids_test = ids_all[0:9]

primers_list = []
for i in ids_test:
    target = str(seq_dict.get(i))
    seq_arg_batch = {
        'SEQUENCE_ID': 'test_leg',
        'SEQUENCE_TEMPLATE': target,
        # 'SEQUENCE_INCLUDED_REGION': [3, len(test_seq)],
        # 'SEQUENCE_PRIMER': str(fp),
        # 'SEQUENCE_PRIMER_REVCOMP': str(rp),
        'SEQUENCE_FORCE_LEFT_START': 3,
        'SEQUENCE_FORCE_RIGHT_START': len(target)-1,
    }

    global_arg_batch = {
        # 'PRIMER_OPT_SIZE': 51,
        'PRIMER_PICK_ANYWAY': 1,
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        # 'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 16,
        'PRIMER_MAX_SIZE': 28,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 68.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        # 'PRIMER_MAX_POLY_X': 100,
        # 'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        # 'PRIMER_MAX_NS_ACCEPTED': 0,
        # 'PRIMER_MAX_SELF_ANY': 12,
        # 'PRIMER_MAX_SELF_END': 8,
        # 'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        # 'PRIMER_PAIR_MAX_COMPL_END': 8,
        # 'PRIMER_PRODUCT_SIZE_RANGE': [[75, 100], [100, 125], [125, 150], [150, 175], [175, 200], [200, 225]],
    }
    res = primer3.designPrimers(seq_arg_batch, global_arg_batch)
    fp = ftail53 + Seq(res.get('PRIMER_LEFT_0_SEQUENCE'), IUPAC.unambiguous_dna)
    rp = rtail53 + Seq(res.get('PRIMER_RIGHT_0_SEQUENCE'), IUPAC.unambiguous_dna)
    prms = [i, str(fp), str(rp)]
    primers_list.append(prms)




primers_list