class UniProtEntry:

    Instances = {}

    def __init__(self, udata):
        self.accs = accs
        self.source = source


    def __str__(self):
        return "UniProtEntry-" + self.get_accs()

    info = {}

    def get_accs(self):
        return self.info['accs']

    def get_gene_name(self):
        return self.info['gene_name']

    def get_uniprot_id(self):
        return self.info['uniprot_id']

    def get_accs_alt(self):
        return self.source['accs_alt']

    def sequence(self):
        return self.source['sequence']

    def get_refseq(self):
        return self.source['refseq']

    def get_pdb(self):
        return self.source['pdb']

    def get_pdb_all(self):
        return self.source['pdb_all']

    def get_full_name(self):
        return self.source['full_name']

    def get_organism(self):
        return self.source['organism']













class UniProtParser:

    """Retrieve and parse a UniProtEntry file, returning a UniProtEntry instance"""

    def __init__(self):
        self.ids = ids
        self.tag = None
        #self.udata = {}

# Action
    def parse(self, ids):
        ids = self.ids

        udata = self.get_udata()

        return UniProtEntry(udata)

# Supporting
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


    def make_feature_generator(self):
        """Return a generator that produces instances of GenBankFeature
        using the data found in the features section of the file"""
        ttag = l[:5].strip()
        while ttag and ttag != tag:
                tag = ttag
                words = (l[5:].strip()).split()


         uni_txt = self.pull_uni()
        uniprot_id = None
        metadata_by_acc = {}
        tag = None
        for l in uni_txt.splitlines():
            ttag = l[:5].strip()
            if ttag and ttag != tag:
                tag = ttag
            words = (l[5:].strip()).split()






    def is_at_info(self):

    def is_at_source(self):




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
            words = (l[5:].strip()).split())


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




        tag = (line[:5].strip() for line in uni_txt.splitlines()) # generator

        words = (line[5:].strip() for line in uni_txt.splitlines()) ##words gen

        [(line[:5].strip(), line[5:].strip()) for line in (uni_txt.splitlines())] ##tuples
        ((line[:5].strip(), line[5:].strip()) for line in (uni_txt.splitlines())) #generator for tuples
