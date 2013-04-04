'''
Created on Dec 13, 2012

@author: jcmorgan
'''
from Bio import SeqIO
import csv, os, difflib, pickle, sys
import simplejson as json

def findall(sub, string):
    index = 0 - len(sub)
    try:
        while True:
            index = string.index(sub, index + len(sub))
            yield index
    except ValueError:
        pass
    
def list_to_str(l):
    s = ''
    for i in l:
        s += i+', '
    return s

def refresh_generator(file_name, format_type):
    f = open(file_name)
    return SeqIO.parse(f, format_type)

class Full_miRNA_Archive():
    def __init__(self, name, reads, lib, thres=0, sim=0):
        self.thres = thres
        self.sim = sim
        self.id = id(self)
        self.reads = reads
        self.lib = lib
        self.repository = {}
        self.name = name
        self.working_dir = '.'
        self.matched = os.path.join(self.working_dir, 'matched.csv')
        self.unique = os.path.join(self.working_dir, 'unique.csv')
        self._preprocess()
        
    def _preprocess(self):
        self.reads = self.clip_with_thres(self.reads)
        experimental_gen = refresh_generator(self.reads, 'fasta')
        for e in experimental_gen:
            r = Experimental_miRNA(e)
            mi = miRNA_Record()
            mi.set_SEQUENCE(str(r.seq))
            self.repository[str(r.id)] = mi
            del r
        
    def start(self):
        print 'Starting\t', self.lib
        mature = refresh_generator(self.lib, 'fasta')        
        for m in mature:
            experimental = refresh_generator(self.reads, 'fasta')
            com_list = list()
            for e in experimental:
                com = self._shiftFrameAndCheck(m, e)
                if com:
                    com_list.append(com)

            for c in com_list:
                obj = self.repository[c[1]]
                obj.add_SPECIES((m.description, str(m.seq)))
                obj.set_RECORDS([x[1] for x in com_list])
                obj.add_COUNT(sum([x[0] for x in com_list]))            
        #self.display()
    
    def display(self):
        experimental = refresh_generator(self.reads, 'fasta')
        for e in experimental:
            r = Experimental_miRNA(e)
            obj = self.repository[r.id]
            if len(obj.get_RECORDS()) is not 0:
                print obj.str_repr()
   
    def clip_with_thres(self,unclipped):
        u = 'clipped_%s.fasta' % (str(self.id))
        name = os.path.join(self.working_dir,u)
        exp_gen = SeqIO.parse(unclipped, 'fasta')
        with open(name, 'a+') as f:
            for e in exp_gen:
                if int(e.description.split('-')[1]) >= self.thres:
                    SeqIO.write(e, f, 'fasta')
                else:
                    break
        return name

    def get_initial_conditions(self,reads,lib):
        mature_gen = None
        m_total = 0
        if lib:
            mature_gen = SeqIO.parse(open(lib), 'fasta')
            m_total = sum([1 for m in mature_gen])
        
        exp_gen = SeqIO.parse(open(reads), 'fasta')
        e_total = sum([1 for e in exp_gen])
             
        return (m_total, e_total)
             
    def _shiftFrameAndCheck(self, mature, experimental):
        m = Mature_miRNA(mature)
        e = str(experimental.seq)
        s = difflib.SequenceMatcher()
        s.set_seq1(m.body)
        exp_count = 0
        exp_records = []
        frames = len(e) - len(m.seq) + 1
        l = experimental.description.split('-')
        if frames > 0:
            for f in range(frames):
                b = len(m.body)
                eb = e[9+f:9+f+b]
                s.set_seq2(eb)
                if m.header==e[1+f:9+f] and s.ratio() >= float(self.sim):
                    exp_count = int(l[1])
                    exp_records = l[0]
                    return (exp_count, exp_records)
        else:
            return None
        
    def preserve(self):
        val = pickle.dumps(self)
        with open(os.path.join('Results', self.name, '%s.pyk'%(self.name) ), 'w+') as f:
            f.write(val)

class miRNA_Record():
    '''Design of data structure:
    {   
        <ID>:
        <<<<{
            SEQUENCE: <EXPERIMENTAL miRNA>,
            INFO:
            [{   
                SPECIES:    [(
                                <SPECIES1>,    
                                <SPECIESmiRNA>
                            )],
                COUNT:      [...],
                RECORDS:    [...]
            }]>>>>
    }
    '''
    def __init__(self):
        self._repo = {'SEQUENCE':None,'INFO':{'SPECIES':[],'COUNT':0,'RECORDS':[]}}
        
    def set_SEQUENCE(self, seq):
        self._repo['SEQUENCE'] = str(seq)
        
    def get_SEQUENCE(self):
        return self._repo['SEQUENCE']
        
    def set_SPECIES(self, species):
        if type(species)==type(list()):
            self._repo['INFO']['SPECIES'] = species
        else:
            raise ValueError,'Species defined must be a list.'
        
    def add_SPECIES(self, species):
        if species not in self._repo['INFO']['SPECIES']:
            self._repo['INFO']['SPECIES'].append(species)
        
    def get_SPECIES(self ):
        return self._repo['INFO']['SPECIES']
    
    def set_RECORDS(self, records):
        if type(records)==type(list()):
            self._repo['INFO']['RECORDS'] = records
        else:
            raise ValueError,'Records defined must be a list.'
        
    def add_RECORD(self, record):
        if record not in self._repo['INFO']['RECORDS']:
            self._repo['INFO']['RECORDS'].append(record)
            return True
        else:
            return False
        
    def get_RECORDS(self ):
        return self._repo['INFO']['RECORDS']
        
    def add_COUNT(self, count):
        self._repo['INFO']['COUNT'] += count
        
    def reset_COUNT(self):
        self._repo['INFO']['COUNT'] = 0
        
    def get_COUNT(self):
        return self._repo['INFO']['COUNT']
    
    def str_repr(self):
        return json.dumps(self._repo)

class Mature_miRNA():
    
    def __init__(self, seq_record):
        self.species = seq_record.description
        self.seq = str(seq_record.seq)
        self.header = str(seq_record.seq)[1:9]
        self.body = str(seq_record.seq)[9:]

class Experimental_miRNA():
    
    def __init__(self, record):
        l = record.description.split('-')
        self.id = l[0]
        self.count = int(l[1])
        self.seq = str(record.seq)
        
class TableFile():
    
    def __init__(self, name, archive, erase=False):
        if os.path.isfile(name) and not erase:
            self.file = open(name, 'ab')
        else:
            self.file = open(name, 'wb')
        self.archive = archive

    def get_archive(self):
        return self.__archive

    def set_archive(self, value):
        self.__archive = value

    def del_archive(self):
        del self.__archive

    archive = property(get_archive, set_archive, del_archive, "archive's docstring")
    
class ReadsCountFile(TableFile):
    def __init__(self, name, archive):
        TableFile.__init__(self, name, archive)
        
    def render(self):
        writer = csv.writer(self.file, delimiter='\t')
        for e in self.archive.repository.keys():
            r = self.archive.repository[e]
            writer.writerow([r.get_SEQUENCE(), e, r.get_COUNT()])
    
class MatureCountFile(TableFile):
    def __init__(self, name, lib, archive):
        TableFile.__init__(self, name, archive)
        self.lib = lib
        
    def render(self):
        f = open(self.lib,'rb')
        s = SeqIO.parse(f,'fasta')
        
        file_list = dict()
        writer = csv.writer(self.file, delimiter='\t')
        
        for m in s:
            ident = str(m.description).split()[0]
            file_list[ident]=[str(m.seq),0,[]]
        
        for m in self.archive.repository.keys():
            spec = self.archive.repository[m].get_SPECIES()
            if len(spec) is not 0:
                for s in spec:
                    ident = s[0].split()[0]
                    file_list[ident][1] = self.archive.repository[m].get_COUNT()
                    file_list[ident][2] = self.archive.repository[m].get_RECORDS()
        
        for fl in file_list.keys():
            if file_list[fl][1] is not 0:
                info = file_list[fl]
                writer.writerow([fl, info[0], info[1], info[2]])
    
class ReadsMatureFile(TableFile):
    def __init__(self, name, archive):
        TableFile.__init__(self, name, archive)
        
    def render(self):
        writer = csv.writer(self.file, delimiter='\t')
        for e in self.archive.keys():
            r = self.archive[e]
            writer.writerow([r.get_SEQUENCE(), e, r.get_COUNT()])

def main(lib, reads, thres, sim):
    
    info = 'mature_libraries'
    species = os.path.basename(lib)[:-3]
    name = species

    try:
        os.mkdir('Results/%s'%species)
    except:
        pass
    
    
    archive = Full_miRNA_Archive(name, reads, lib, thres, sim)
    archive.start()
    archive.preserve()

    path = os.path.join('Results',species,'ReadsCountFile_%s.tab'%species)
    rct = ReadsCountFile(path , archive)
    rct.render()

    path = os.path.join('Results',species,'MatureCountFile_%s.tab'%species)
    mct = MatureCountFile(path, lib, archive)
    mct.render()
    
    clean_up = True
    if clean_up:
        rm = 'rm Results/clipped_*' 
        os.system(rm)
            
if __name__=='__main__':
    
##    library = os.path.join(info,'Acanthopterygii')
##    reads = os.path.join(info,'rna_dna.fasta')
##    thres = 500
##    sim = .95
    library = sys.argv[1]
    reads = sys.argv[2]
    thres = sys.argv[3]
    sim = sys.argv[4]
    main(library, reads, thres, sim)

    
