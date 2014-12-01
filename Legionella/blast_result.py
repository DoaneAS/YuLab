from Bio import SearchIO
result_handle = open("results/blastp_result.txt")
idx = SearchIO.index('results/blastp_result.txt', 'blast-tab')# comments=True)
idss = (idx.keys())
homs = {}
rec_list_all = []

for ids in idss:
    rec_list = []
    for rec in idx[ids].hsps:
        rec_list.append((rec.hit_id, rec.evalue, rec.ident_pct))
        rec_list_all.append((rec.query_id, rec.hit_id, rec.evalue, rec.ident_pct))
    homs[ids] = {ids: rec_list}
leg_uniq = []
ids = []
for p in sorted(rec_list_all, key=lambda x: x[0][1]):
    if p[2] < 10E-5 and p[1] not in p[0]:
        if p[1] not in ids:
            ids.append(p[1])
            leg_uniq.append(p)

output = open('results/leg_blastp_hits.txt', 'w')
output.write('id'+'\t'+'blastp_hit'+ '\t' + 'e-value' + '\t' + 'ident_pct' + '\n')
for p in leg_uniq:
    qid = '%s' % p[0]
    hit = '%s' % p[1]
    ev = '%s' % p[2]
    per = '%s' % p[3]
    output.write(qid + "\t" + hit + "\t" + ev + '\t' + per + "\n")
output.close()