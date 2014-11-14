from Bio import SearchIO
result_handle = open("results/blastp_result.txt")
idx = SearchIO.index('results/blastp_result.txt', 'blast-tab')# comments=True)
idss = (idx.keys())
homs = {}
rec_list_all = []

for ids in idss:
    rec_list = []
    for rec in idx[ids].hsps:
        rec_list.append((rec.hit_id, rec.evalue))
        rec_list_all.append((rec.hit_id, rec.evalue))
    homs[ids] = {ids: rec_list}

output = open('results/leg_blastp_hits.txt', 'w')
output.write('id'+'\t'+'blastp_hit'+ '\t' + 'e-value' + '\n')
for k in homs.keys():
    key = '%s' % k
    for v in homs[k].values():
        for h in v:
            hit = '%s' % h[0]
            ev = '%s' % h[1]
            output.write(key + "\t" + hit + "\t" + ev + "\n")
output.close()