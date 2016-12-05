# coding:utf-8
import os
from Bio import SeqIO
path = "../Human/sequences"

lens = []

for filename in os.listdir(path):
    for sec in SeqIO.parse(os.path.join(path, filename), 'fasta'):
        length = len(sec.seq)
        idd = sec.id.split("|")[1]
    lens.append('{:5d} {:s}\n'.format(length, idd))

lens.sort()
fw = open("seq_len.txt", 'w')
fw.writelines(lens)
fw.close()

