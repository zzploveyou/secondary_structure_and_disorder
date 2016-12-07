# coding:utf-8
"""
compare evolutionary rate in 
different structures(disorder(D), ambigous(A), coil(C), sheet(E), helix(H))
for each sequence.
"""
import os
import sys
from glob import glob
import json
import log
import math
from Bio import SeqIO

mylog = log.Terminal_log()
PATH = "/home/biolab/zzp/Human"
SS_DIR = os.path.join(PATH, "secondary_structure")
DIS_DIR = os.path.join(PATH, "disorder")
D2P2_DIR = os.path.join(PATH, "D2P2")

RET_DIR = os.path.join(PATH, "RET")




class SS:

    def __init__(self, idd):
        self.idd = idd
        self.matfile = os.path.join(SS_DIR, idd + ".mat")
        self.ssfile = os.path.join(SS_DIR, idd + ".ss2")
        self.disfile = os.path.join(DIS_DIR, idd + ".json")
        # self.ss存放结构信息, D,A,C,E,H
        self.ss = []
        # self.er存放保守率
        self.er = []
        # self.m存放各种结构的长度以及平均er值
        self.m = {'D':[0,0],'A':[0,0],'C':[0,0], 'E':[0,0], 'H':[0,0]}
        
        self.tag = 1
        self.result = {}
        if not os.path.exists(self.ssfile):
            mylog.debug("ssfile not existed.")
            self.tag = 0
        if not os.path.exists(self.disfile):
            mylog.debug("disfile not existed.")
            self.tag = 0

    def readss(self):
        """读取二级结构ss2文件"""
        for line in open(self.ssfile):
            tmp = line.strip().split()
            if len(tmp) == 6:
                self.ss.append(tmp[2])

    def changess(self, start, end, tag):
        """指定位置加入disorder信息"""
        mylog.debug('length:%s start:%s end:%s'%(len(self.ss), start, end))
        for i in range(start, end+1):
            try:
                self.ss[i] = tag
            except IndexError as e:
                """uniprot上序列长度不全.例如Q8N5G0只有67长度,但是disorder却很长"""
                pass

    def readdis(self):
        """读取MobiDB中的disorder信息"""
        #self.disfile = "123.json"
        js = json.load(open(self.disfile))
        disorder = js['consensus']['full']
        for son in disorder:
            ann = son['ann']
            start = int(son['start'])-1
            end = int(son['end'])-1
            if ann == 'A' or ann == 'D' or ann == 'd':
                ann = ann.upper()
                self.changess(start, end, ann)
        
    def printss(self):
        """返回加入disorder信息的序列结构"""
        mylog.debug('length: %d, ss: %s' %(len(self.ss), ''.join(self.ss)))
    
    def conser(self, args):
        """计算保守性"""
        s = 0.
        for i in args:
            if i != 0:
                s += (-i*math.log(i, 2))
        return s 

    def caler(self):
        """读取PSSM矩阵,并计算序列各个残基的保守性"""
        for line in open(self.matfile):
            tmp = line.strip().split()
            if tmp != [] and tmp[0].isdigit():
                pro = (float(i)/100 for i in tmp[22:42])
                self.er.append(self.conser(pro))                
        
    def getmap(self):
        """计算不同结构的长度,以及平均的保守性变化"""
        assert len(self.ss) == len(self.er)
        for ss,er in zip(self.ss, self.er):
            self.m[ss][0] += 1
            self.m[ss][1] += er
        for i in self.m:
            if self.m[i][0] != 0:
                self.m[i][1] /= self.m[i][0]
    def d2p2_exam(self):
        d2p2id = ""
        for line in open(os.path.join(D2P2_DIR, "d2p2_protein_to_uniprot.tsv")):
            tmp = line.strip().split()
            if tmp[1] == self.idd:
                d2p2id = tmp[0]
                break
        if d2p2id == "":mylog.error("not exists d2p2id_uniprot.")
        mylog.debug("d2p2id: %s" %(d2p2id))
        tag = 0
        for rec in SeqIO.parse(os.path.join(D2P2_DIR, "human.fasta"), "fasta"):
            if rec.id == d2p2id:
                tag = 1
                if len(rec.seq) != len(self.ss):
                    mylog.debug("d2p2 length: %s, ss length: %s" % (len(rec.seq), len(self.ss)))
                    
                    mylog.error("d2p2 length not equals us")
                    return False
        if tag == 1:
            mylog.debug("d2p2 length equals us.")
        else:
            mylog.warn("can found d2p2")
        return True

    def run(self):
        self.readss()
        self.readdis()
        self.caler()
        self.getmap()
        self.result['ss'] = self.ss
        self.result['er'] = self.er
        self.result['map'] = self.m
        self.d2p2_exam()
        # print json.dumps(self.result, sort_keys=False, indent=4)
        fp = open(os.path.join(RET_DIR, self.idd+".json"), 'w')
        json.dump(self.result, fp, indent=4)
        

if __name__ == '__main__':        
    for filename in glob(os.path.join(SS_DIR, "Q9NQG1.mat")):
        idd = os.path.basename(filename)[:-4]
        mylog.info(idd)
        s = SS(idd)
        if s.tag == 1:
            s.run()
        raw_input()
        


