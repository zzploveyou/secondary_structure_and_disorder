# coding:utf-8
from Bio import SeqIO
import os
HD = "/home/biolab/zzp/Human/"
d2p2_protein_to_uniprot_filename = os.path.join(HD, "D2P2/d2p2_protein_to_uniprot.tsv")
d2p2_fastafile = os.path.join(HD, "D2P2/human.fasta")
uni_fastafile = os.path.join(HD, "UP000005640_9606.fasta")
ptm_filename = os.path.join(HD, "D2P2/ptm_assignment.tsv")
ptm_result = os.path.join(HD, "D2P2/PTM_human.txt")

def getidlens(fastafile, tag="normal"):
    m = {}
    for rec in SeqIO.parse(fastafile, "fasta"):
        if tag == "d2p2":
            m[rec.id] = len(rec.seq)
        elif tag == "normal":
            m[rec.id.split("|")[1]] = len(rec.seq)
    return m

def getmap():
    """return a dic, map uniprot id to d2p2 id"""
    human_IDs = ""
    m = {}
    for line in open(d2p2_protein_to_uniprot_filename):
        tmp = line.strip().split()
        m[tmp[1]] = tmp[0]
        # print tmp[0], tmp[1]
    return m

def main():
    m_d2p2 = getidlens(d2p2_fastafile, tag="d2p2")
    ma = getmap()
    m_uni = getidlens(uni_fastafile, tag="normal")
    id_map = {}
    
    for idd in m_uni.keys():
        if idd in ma:
            if ma[idd] in m_d2p2:
                if m_d2p2[ma[idd]] == m_uni[idd]:
                    id_map[ma[idd]]=idd
    lines = []
    for line in open(ptm_filename):
        tmp = line.strip().split("\t")
        """在ptm file中搜索符合human的条目"""
        if tmp[1] in id_map:
            lines.append("%s\t%s\t%s\t%s\n" % (id_map[tmp[1]], tmp[3], tmp[4], tmp[5]))
    fw = open(ptm_result, 'w')
    fw.writelines(lines)
    fw.close()

if __name__ == '__main__':
    main()
    
