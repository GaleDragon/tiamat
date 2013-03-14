from Bio import SeqIO
import os

def create_headers():
    mat = SeqIO.parse('mature_libraries/mature.fa','fasta')
    s_header = ''
    prefix = ''
    f = ''
    for s in mat:
        s_header = str(s.seq)[1:8]
        prefix = s_header[:3]
        suffix = s_header[3:]
        p_folder = os.path.isdir('Headers/%s'%(prefix))
        s_folder = os.path.isdir('Headers/%s/%s'%(prefix,suffix))
        
        if not p_folder:
            os.system('mkdir Headers/%s'%(prefix))

        if not s_folder:
            os.system('mkdir Headers/%s/%s'%(prefix,suffix))

        
        f = open('Headers/%s/%s/%s.fa'%(prefix,suffix,s_header),'a+')
        SeqIO.write(s,f,'fasta')

if __name__=='__main__':
    create_headers()