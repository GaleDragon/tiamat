from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqRecord
import difflib, csv, os, sys
from time import gmtime, strftime

# define globals
matched_records = []
matched_seqs = []
species = []

class Mature_miRNA():
    
    def __init__(self, seq):
        self.seq = str(seq)
        self.header = str(seq)[1:9]
        self.body = str(seq)[9:]


def findall(sub, string):
    index = 0 - len(sub)
    try:
        while True:
            index = string.index(sub, index + len(sub))
            yield index
    except ValueError:
        pass

   
def closer():
##def improvedSFAC(mature, experimental, thres):
##    ## Uses the findall method to return a list of indices of a specific substring ie the mature header
##    m = Mature_miRNA(mature.seq)
##    e = str(experimental.seq)
##
##    s = difflib.SequenceMatcher()
##    s.set_seq1(m.body)
##
##    e_len = len(e)
##
##    exp_count = 0
##    exp_records = []
##
##    l = experimental.description.split('-')
##
##    for ind in findall(m.header, e):
##        offset = ind - 1
##        s.set_seq2(e[offset:len(m.seq)+offset])
##        a = '-'*offset+m.seq
##        tail = e_len - len(a)
##        offshoot = False
##        if len(a+'-'*tail) > e_len:
##            offshoot = True
##        else:
##            print a+'-'*tail
##            print e
##
##        if s.ratio() > .95 and not offshoot:
##            exp_count = int(l[1])
##            exp_records = l[0]
##            break
##    return (exp_count, exp_records)
    pass


def shiftFrameAndCheck(mature, experimental, thres, similarity):
    m = Mature_miRNA(mature.seq)
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
            if m.header==e[1+f:9+f] and s.ratio() >= float(similarity):
                exp_count = int(l[1])
                exp_records = l[0]
                break
    return (exp_count, exp_records)


def getLaunchPrep(lib,reads,thres):
    mature_gen = None
    if lib:
        mature_gen = SeqIO.parse(open(lib), 'fasta')
    exp_gen = SeqIO.parse(open(reads), 'fasta')
    m_total = 0
    e_total = 0
    if lib:
        for m in mature_gen:
            m_total += 1
    for e in exp_gen:
        e_total += 1
    return (m_total, e_total)


def list_to_str(list):
    s = ''
    for i in list:
        s += i+', '
    return s


def begin_taxon(lib,reads,thres, sim):
    # bring globals into namespace
    global matched_records
    global matched_seqs
    global species
    # setup the csv.writer
    writer = csv.writer(open('../Results/matched.csv','w+'), delimiter='\t')
    writer.writerow(['DESCRIPTION','SEQUENCE','TOTAL','RECORDS'])
    # create the generator for the
    mature_gen = SeqIO.parse(open(lib), 'fasta')
    #
    for m in mature_gen:
        exp_gen = SeqIO.parse(open(reads), 'fasta')
        total = 0
        record_list = []
        records = list()
        for e in exp_gen:
            r = shiftFrameAndCheck(m,e, thres, sim)
            if bool(r[0]):
               record_list.append(r)
        # Adds the total miRNA counts up and creates the list of records
        # those counts were created from.
        for t in record_list:
            total += t[0]
            records.append(t[1])
            # Doesn't add the record to the final list of all matched
            if t[1] not in matched_records:
                matched_records.append(t[1])
        if total > 0:
            writer.writerow([m.description, str(m.seq), str(total), list_to_str(records)[:-2]])
    writer2 = csv.writer(open('../Results/unique.csv','w+'), delimiter='\t')
    writer2.writerow(['DESCRIPTION','SEQUENCE','TOTAL','RECORDS'])
    exp_gen = SeqIO.parse(open(reads), 'fasta')
    for e in exp_gen:
        id_num = e.description.split('-')
        if id_num[0] not in matched_records:
            writer2.writerow(["Unique_miRNA_{0:0>10}".format(id_num[0]), str(e.seq), str(id_num[1]), id_num[0]])


def clip_with_thres(tol, unclipped):
    name = '../Results/clipped.fasta'
    exp_gen = SeqIO.parse(unclipped, 'fasta')
    with open(name, 'w+') as f: pass
    with open(name, 'a+') as f:
        for e in exp_gen:
            if int(e.description.split('-')[1]) > int(tol):
                SeqIO.write(e, f, 'fasta')
            else:
                break
    return name


if __name__=='__main__':
#    library = sys.argv[1]
#    reads = sys.argv[2]
#    tol = sys.argv[3]
#    similarity = sys.argv[4]
    print os.getcwd()
    library = '../mature_libraries/Acanthopterygii'
    reads = '../mature_libraries/rna_dna.fasta'
    tol = 750
    similarity = .95
    clipped = clip_with_thres(tol, reads)
    #print clipped
    
    (mt, et) = getLaunchPrep(library,clipped,tol)
    print "Mature sequences:",mt
    print "Experimental sequences:",et
    
    begin_taxon(library, clipped, tol, similarity)
    ROOT = "Headers"
    #begin_traverse(ROOT, reads, tol)