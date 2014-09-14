class InParndOrthologs:
    """docstring for InParndOrthologs"""

    Instances = {}

    #def __init__(self, ref_sp, orth_sp, ref_prot, orth_prot, ortholog, source):
    def __init__(self, ortholog_d):
        #self.ref_sp = ref_sp
        #self.orth_sp =orth_sp
        #self.source = source
        #self.ref_prot = ref_prot
        #self.orth_prot = orth_prot
        self.ortholog_d = ortholog_d
        #self.ortholog = (ref_prot, orth_prot)
        #self.Instances[(self.ref_prot, self.orth_prot)] = self  #???

class Ortholog_d(object):
    '''dict of orthologs'''

    def __init__(self, ortholog_d):
        self.ortholog_d = ortholog_d

    def get_ortholog_d(self):
        return self.ortholog_d



class Ortholog(object):
    """docstring for Ortholog"""

    def __init__(self, ortholog):
        self.ortholog = ortholog


    def info(self):
        return (self.ortholog)

    def ref_Prot(self):
        return self.ortholog[0]

    def orth_Prot(self):
        return self.ortholog[1]

    def get_ortholog(self):
        return self.ortholog

    def __repr__(self):
        return '[Ortholog: %s, %s]' % (self.ortholog)



class InParanoidParser(object):
    """Parse a InParanoid txt file, returning an InParanoid instance"""

    def __init__(self):
        self.line = None
        self.linenum = 0
        self.src = None

#Predicates

    def is_at_ref(self):
        return self.line and "%" in self.line[24:31]

    def is_at_orth(self):
        return self.line and "%" in self.line[65:67]

    def is_at_orth_end(self):
        return self.line and "%" not in self.line[65:67]

    def is_at_end(self):
        return self.line and self.line.startswith('//') ##???


#Action

    def parse(self, filename):
        """Use the data in file named filename to create and
        return an instance of InParndOrthologs"""
        with open(filename) as self.src:
            while True:
                self.line = self.src.readline()
                line = self.line
                ortholog = self.get_ortholog()
                if not self.line: break
            return Ortholog(ortholog)

    def parser(self, filename):
        with open(filename) as self.src:
            while True:
                self.line = self.src.readline()
                ortholog = self.get_ortholog()
                ortholog_d = self.get_ortholog_d()
                return Ortholog_d(ortholog_d)

                              # until end-of-file
                                  # call a function object



    def read_line(self):
        self.line = self.src.readline()

#Action support

    def get_ortholog(self):
        ref_prot = self.get_ref_prot()
        orth_prot = self.get_orth_prot()
        ortholog = (ref_prot, orth_prot)
        return Ortholog(ortholog)


    def get_ref_prot(self):
        self.read_line()
        while not self.is_at_ref():
            self.read_line()
        ref = self.line[0:6].strip()
        return ref

    def get_orth_prot(self):
        orth = []
        #self.read_next_line()
        while not self.is_at_orth():
            self.read_line()
        while not self.is_at_orth_end():
            orth.append(self.line[35:43].strip())
            self.read_line()
        return orth

    # def get_ortholog_d(self):
    #     ddp = {}
    #     for k,v in self.orth_gen():
    #         ddp[k] = {
    #         'ref': k,
    #         'orth': v
    #         }
    #     return ddp





    def orth_parse(self, filename):
        with open(filename) as self.src:
            src = self.src
            ca = None
            ddp = {}
            for line in src:
                if "%" in line[24:31]:
                    ca = line[0:6].strip()
                ddp[ca] = {
                        "ca":ca,
                        "hu": set(),
                        "geneid":set()
            }
                entry = ddp[ca]
            if "%" in line[65:67]:
                entry['hu'].add(line[35:43].strip())
        ortholog_d= ddp
        return Ortholog_d(ortholog_d)

    #def get_ref_prot(self):



if __name__ == '__main__':
    ca_hu = InParanoidParser().orth_parse("CA_HU.txt")
    import pprint
    pprint.pprint(ca_hu.ortholog_d)


