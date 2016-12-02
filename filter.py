# coding:utf-8
"""
filter some id from fasta into a new fasta.

filter uniref90 from uniref100
according ids.
"""
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def filter_fasta(fastafile, filters):
    recs = []
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        ID = seq_record.id
        seq = seq_record.seq
        des = seq_record.description
        rec = SeqRecord(seq, id=ID, description=des)
        idd = ID.split("|")[1]
        if idd in filters:
            recs.append(rec)
    path = os.path.dirname(fastafile)
    SeqIO.write(recs, os.path.join(path, "result.fasta"), "fasta")

def main():
    filters = set()
    for line in open(""):
        filters.add(line.strip())

    filter_fasta("/home/zzp/DATABASE/Human/.fasta", filters)

if __name__ == "__main__":
    main()
