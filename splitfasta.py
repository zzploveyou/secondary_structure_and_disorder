#coding:utf-8
"""
split fastafile into single sequence fastafile.
length>30 and no X in seq.
"""

import os
import log
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def splitfasta(fastafile, sourcedir):
    mylog = log.Terminal_log()
    if not os.path.exists(sourcedir):
        os.makedirs(sourcedir)
        mylog.info("make sourcedir.")
    if not os.path.exists(fastafile):
        mylog.fatal("fastafile not existed.")
        sys.exit(1)
    for seq_record in SeqIO.parse(fastafile, 'fasta'):
        seq = seq_record.seq
        idd = seq_record.id.split("|")[1]
        des = seq_record.description
        # filter
        if len(seq) >= 30 and 'X' not in seq:
            rec = SeqRecord(seq,
            id=seq_record.id,
            description=des)
            target = os.path.join(sourcedir, "%s.fasta" %idd)
            if os.path.exists(target):
                mylog.error("%s existed." % target)
            SeqIO.write([rec], target, 'fasta')
            mylog.debug("%s splited." % idd)
    mylog.done("done.")

if __name__ == '__main__':
    fastafile = "/home/biolab/zzp/Human/UP000005640_9606.fasta"
    sourcedir = "/home/biolab/zzp/Human/sequences/"
    splitfasta(fastafile, sourcedir)
