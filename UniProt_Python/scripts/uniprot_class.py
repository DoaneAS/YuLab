import requests
import sys
import argparse
import json


class UniProtEntry:

    """Represents protein sequence and annotation data for specified UniProt ID
    """
    #instances = {}

    # def __init__(self, data):
    #     self.accs, self.uniprot_id = data[0][1].get('accs')[0], data, data[0][0]
    #     self.sequence = data[0][1].get('sequence')
    #     self.accs_alt = data[0][1].get('accs')[1:]
    #     self.gene_name = data[0][1].get('gene')
    #     self.refseq = data[0][1].get('refseq')
    #     self.length = data[0][1].get('length')
    #     self.pdb, self.pdb_all = data[0][1].get('pdb'), data[0][1].get('pdbs')
    #     self.full_name = data[0][1].get('full_name')
    #     self.organism = data[0][1].get('organism')

    def __init__(self, data):
        self.accs, self.uniprot_id = accs, uniprot_id
        self.accs_alt = accs_alt
        self.gene_name = gene_name
        self.refseq = refseq
        self.length = length
        self.pdb, self.pdb_all = dpdb, pdb_all
        self.full_name = full_name
        self.organism = organism
        self.sequence = sequence
        self.Instances[self.accs] = self

    def __str__(self):
        return "UniProtEntry-" + self.get_accs()

    def get_accs(self):
        return self.accs

    def get_accs_alt(self):
        return self.accs_alt

    def sequence(self):
        return self.sequence

    def get_gene_name(self):
        return self.gene_name

    def get_refseq(self):
        return self.refseq

    def get_pdb(self):
        return self.pdb, self.pdb_all

    def get_full_name(self):
        return self.full_name

    def get_organism(self):
        return self.organism


class UniProtParser:

    """Retrieve and parse a UniProtEntry file, returning a UniProtEntry instance"""


    # def __init__(self):
        #self.accs = accs
        #self.uniprot_id =  uniprot_id
        #self.metadata_by_acc = metadata_by_acc
        #self.tag = tag

    def __init__(self):
        self.ids = ids
        self.url = 'http://www.uniprot.org/'
        #self.uniprot_id = None
        #self.metadata_by_acc = {}
        #self.tag = None

# Action


    def parse(self, ids):
        ids = self.ids
        uni_txt = self.pull_uni(ids)
        data = self.parse_uniprot()
        accs, uniprot_id = data[0][1].get('accs')[0], data, data[0][0]
        sequence = data[0][1].get('sequence')
        accs_alt = data[0][1].get('accs')[1:]
        gene_name = data[0][1].get('gene')
        refseq = data[0][1].get('refseq')
        length = data[0][1].get('length')
        pdb, pdb_all = data[0][1].get('pdb'), data[0][1].get('pdbs')
        full_name = data[0][1].get('full_name')
        organism = data[0][1].get('organism')
        return UniProtEntry(accs, uniprot_id, accs_alt, gene_name, refseq, length, pdb, pdb_all, full_name, organism, sequence)


# Access Methods
    #def uni_txt(self):
      #  return self.pull_uni(ids)


#Action support
    def pull_uni(self, ids, format='txt'):
        """ request entries by uniprot acc using batch retrieval

        Args:
            ids: list of ids to retrieve
            format: txt by default
            possible formats:
            txt, xml, rdf, fasta, gff"""
        ids = self.ids
        url = self.url
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
        return self.uni_txt






    def parse_uniprot(self):
        """Parse metadata text from uniprot.org

        Returns a dictionary with the UNIPROT ACC as keys.
        """
        uni_txt = self.uni_txt
        uniprot_id = None
        metadata_by_acc = {}
        tag = None
        # with get_parts
        for l in uni_txt.splitlines():
            ttag = l[:5].strip()
            if ttag and ttag != tag:
                tag = ttag
            words = (l[5:].strip()).split()
            get_ids(tag)



        def get_id(self):   #filter
            if tag == "ID":
                uniprot_id = words[0]
                is_reviewed = words[1].startswith('Reviewed')
                length = int(words[2])
                metadata_by_acc[uniprot_id] = {
                    'id': uniprot_id,
                    'is_reviewed': is_reviewed,
                    'length': length,
                    'sequence': '',
                    'accs': [],
                    'gene': [],
                    'refseq': []
                }
                entry = metadata_by_acc[uniprot_id]

            if tag == "SQ":
                if words[0] != "SEQUENCE":
                    entry['sequence'] += ''.join(words)

            if tag == "AC":
                accs = [w.replace(";", "") for w in words]
                entry['accs'].extend(accs)

            if tag == "GN":
                if 'gene' not in entry:
                    entry['gene'] = []
                parts = words[0].split("=")
                entry['gene'].append(parts[1])

            if tag == "DR":
                if 'PDB' in words[0]:
                    if 'pdb' not in entry:
                        entry['pdb'] = [words[1][:-1]]
                    if 'pdbs' not in entry:
                        entry['pdbs'] = []
                    entry['pdbs'].append(words[1][:-1])
                if 'RefSeq' in words[0]:
                    if 'refseq' not in entry:
                        entry['refseq'] = []
                    self = [w[:-1] for w in words[1:]]
                    entry['refseq'].extend(self)

            if tag == "DE":
                # if "RecName" in words[0]:
                if 'full_name' not in entry:
                    fid = [w.replace("Full=", " ") for w in words[1:]]
                entry['full_name'] = [" ".join(fid[0:])]

            if tag == "OS":
                if 'organism' not in entry:
                    oid = [w.replace("Full=", " ") for w in words[0:]]
                entry['organism'] = [" ".join(oid[0:])]

        return metadata_by_acc

if __name__ == '__main__':
    ids = ["P01116"]
    up = UniProtParser().parse(ids)
    print(id)
# class UsingDictionary:
#     Instances = {}
#     def __init__(self[, arg, ...]):
#         self.Instances[self.accs] = self
