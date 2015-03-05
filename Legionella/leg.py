import os
from Bio import SeqIO
from Bio import Entrez

###########
####Retrieve Genome###
####

Entrez.email = "Ashley.Doane@gmail.com"     # Always tell NCBI who you are
handle = Entrez.efetch(db="Protein", id="CP003885",
                       rettype="gp",#or gp for prot sequence,
                       retmode="text")


filename = "Legionella.gbk"
# if not os.path.isfile(filename):
#     # Downloading...
#     net_handle = Entrez.efetch(db="protein", id="CP003885",
#                        rettype="gp",#or gp for prot sequence,
#                        retmode="text")
#     out_handle = open(filename, "w")
#     out_handle.write(net_handle.read())
#     out_handle.close()
#     net_handle.close()
#     print("Saved")

# print("Parsing...")
# record = SeqIO.read(filename, "genbank")
##local input
record = SeqIO.read(filename, "genbank")
##record = SeqIO.read("Legioneslla_pneumophila_Phil/CP003885.gbk")
#print record

#local source input_handle = open("Legionella_pneumophila_Phil/CP003885.gbk")



predicted_genes=[]
leg_genes = {}
for i,feature in enumerate(record.features):
    if feature.type=='CDS':
        product=feature.qualifiers['product'][0].lower()
        seq=feature.extract(record.seq)
        locus = feature.qualifiers['locus_tag'][0]
        predicted_genes.append(feature)
        trans = feature.qualifiers['translation']
        leg_genes[locus] = {"locus": locus, "seq": seq, "product": product, "translation": trans}


filename = "substrates.txt"
with open(filename) as src:
    p= [line.split('\r') for line in src][0]
legdict = {}
for line in p:
    pp = line.split()
    lid = pp[0].lower(),
    alias = pp[1]
    legdict[lid[0]] = {
    "lid": lid[0],
    "alias": alias}

leg_comb = {}
for k, v in legdict.items():
    #print v['lid']
    if k in leg_genes:
        #print leg_genes
        leg_comb[k] = {"lid":v['lid'], "product": leg_genes[k]['product'], "sequence": leg_genes[k]['seq'], "translation": leg_genes[k]['translation']}
print leg_comb

fa = [rec.seq for rec in SeqIO.parse("Legionella_pneumophila_Phil/CP003885.faa", "fasta")]
hand = open("leg.fasta", "w")
sequence = SeqIO.parse("Legionella_pneumophila_Phil/CP003885.faa", "fasta")
SeqIO.write(sequence, hand, 'fasta')
hand.close()
leg_p = SeqIO.parse("Legionella_pneumophila_Phil/CP003885.faa", "fasta")




######
##Create protein fasta file from gbk input####
#for whole genome#

from Bio import GenBank

gbk_filename = "Legionella.gbk"
faa_filename = "Legionella_all_converted.faa"

output_handle = open(faa_filename, "w")
input_handle  = open(gbk_filename, "r")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print "Dealing with GenBank record %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s from %s\n%s\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()

input_handle = open(faa_filename)
recfa = SeqIO.to_dict(SeqIO.parse(faa_filename, "fasta"))


##Create protein fasta file from gbk input####
#for selected genes only#

gbk_filename = "Legionella.gbk"
faa_filename = "results/Legionella_selected.faa"

output_handle = open(faa_filename, "w")
input_handle  = open(gbk_filename, "r")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print "Dealing with GenBank record %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            if not legdict.get(seq_feature.qualifiers['locus_tag'][0]) == None:
                output_handle.write(">%s from %s\n%s\n" % (
                       seq_feature.qualifiers['locus_tag'][0],
                       seq_record.name,
                       seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()

input_handle = open(faa_filename, "r")
rec_dict = SeqIO.to_dict(SeqIO.parse(input_handle, "fasta"))

print len(rec_dict)
