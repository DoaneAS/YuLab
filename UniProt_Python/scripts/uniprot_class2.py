import requests
import sys
import argparse
import json




class UniProtEntry:

    Instances = {}

    def __init__(self, udata):
        self.accs = udata[0][1].get('accs')[0]
        self.udata = udata
        self.sequence = udata[0][1].get('sequence')
        self.accs_alt = udata[0][1].get('accs')[1:]
        self.gene_name = udata[0][1].get('gene')
        self.refseq = udata[0][1].get('refseq')
        self.length = udata[0][1].get('length')
        self.pdb = udata[0][1].get('pdb')
        self.pdb_all = udata[0][1].get('pdbs')
        self.full_name = udata[0][1].get('full_name')
        self.organism = udata[0][1].get('organism')
        self.uniprot_id = udata[0][0]
        self.Instances[self.accs] = self

    @classmethod
    def GetInstances(self): # returning generator as discussed earlier
        return (value for value in self.Instance.keys())

    def __str__(self):
        return "UniProtEntry-" + self.get_uniprot_id()

    def accs(self):
        return self.udata[0][1].get('accs')[0]

    def uniprot_id(self):
        return self.udata[0][0]

    def accs_alt(self):
        return self.udata[0][1].get('accs')[0], udata, udata[0][0]

    def gene_name(self):
        return self.udata[0][1].get('gene')

    def length(self):
        self.udata[0][1].get('length')

    def full_name(self):
        return udata[0][1].get('full_name')

    def sequence(self):
        return self.udata[0][1].get('sequence')

    def organism(self):
        return self.udata[0][1].get('organism')

    def pdb(self):
        return self.udata[0][1].get('pdb')

    def pdb_all(self):
        return self.udata[0][1].get('pdbs')


class UniProtParser:

    """Retrieve and parse a UniProtEntry file, returning a UniProtEntry instance"""


    def __init__(self):
        #self.ids = ids
        self.tag = None
        #self.udata = {}


#Action

    def parse(self, ids):
        ids = self.ids

        udata = self.get_udata()

        return UniProtEntry(udata)



#Supporting


    def pull_uni(self, format='txt'):
        """ request entries by uniprot acc using batch retrieval

        Args:
            ids: list of ids to retrieve
            format: txt by default
            possible formats:
            txt, xml, rdf, fasta, gff"""
        ids = self.ids
        url = 'http://www.uniprot.org/'
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
        return uni_txt



    def get_udata(self):
        """Parse metadata text from uniprot.org

        Returns a dictionary with the UNIPROT ACC as key.
        If >1 id, return list of dictionaries
        """
        # with get_parts
        uni_txt = self.pull_uni()
        uniprot_id = None
        metadata_by_acc = {}
        tag = None
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

        return metadata_by_acc.items()





if __name__ == '__main__':
    #ids = ["P01116"]
    ids = ['P01116', 'P15056']
    import sys
    import pprint
    up = UniProtParser().parse(ids)
    inst = up.Instances
    print inst
    #pprint.pprint(up)

