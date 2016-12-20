# coding:utf-8
"""
compare evolutionary rate in 
different structures(disorder(D), ambigous(A), coil(C), sheet(E), helix(H))
for each sequence.
in mobidb, C means Anbiguous.
"""
import os
import sys
from glob import glob
import json
import log
import math
from Bio import SeqIO

mylog = log.Terminal_log(brief=True)
PATH = "/home/biolab/zzp/Human"
SS_DIR = os.path.join(PATH, "secondary_structure")
DIS_DIR = os.path.join(PATH, "disorder")
D2P2_DIR = os.path.join(PATH, "D2P2")
ANCHOR_DIR = os.path.join(PATH, "ANCHOR")

RET_DIR = os.path.join(PATH, "RET")




class SS:

    def __init__(self, idd):
        self.idd = idd
        self.matfile = os.path.join(SS_DIR, idd + ".mat")
        self.ssfile = os.path.join(SS_DIR, idd + ".ss2")
        self.disfile = os.path.join(DIS_DIR, idd + ".json")
        # self.seq存放序列
        self.seq = []
        # self.ss存放结构信息, D,A,C,E,H
        self.ss = []
        # self.er存放保守率
        self.er = []
        # self.anchor存放anchor概率
        self.anchor = []
        # self.m存放各种结构的长度以及平均er值
        self.m = {'D':[0,0],'A':[0,0],'C':[0,0], 'E':[0,0], 'H':[0,0]}
        # self.ptm 存放PTM信息
        self.ptm = []
        
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
                self.seq.append(tmp[1])
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
        m = {'C':'A', 'd':'D', 's':'S', 'D':'D', 'S':'S'}
        for son in disorder:
            ann = son['ann']
            start = int(son['start'])-1
            end = int(son['end'])-1
            ann = m[ann]
            if ann == 'A' or ann == 'D':
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
        return s/math.log(20, 2) 

    def caler(self):
        """读取PSSM矩阵,并计算序列各个残基的保守性"""
        for line in open(self.matfile):
            tmp = line.strip().split()
            if tmp != [] and tmp[0].isdigit():
                pro = (float(i)/100 for i in tmp[22:42])
                self.er.append(self.conser(pro))                
        
    def getmap(self):
        """计算不同结构的长度,以及平均的保守性变化"""
        try:
            assert len(self.ss) == len(self.er)
        except Exception as e:
            mylog.fatal(e)
        for ss,er in zip(self.ss, self.er):
            self.m[ss][0] += 1
            self.m[ss][1] += er
        for i in self.m:
            if self.m[i][0] != 0:
                self.m[i][1] /= self.m[i][0]
    
    def readanchor(self):
        for line in open(os.path.join(ANCHOR_DIR, self.idd+".out")):
            if line[0] != "#":
                self.anchor.append(float(line.strip().split()[2]))
    
    def getptm(self):
        self.ptm = ['None'] * len(self.seq)
        for line in open(os.path.join(D2P2_DIR, "PTM_human.txt" )):
            tmp = line.strip().split("\t")
            idd, ty, pos, ami = tmp
            pos = int(pos) - 1
            if idd == self.idd:
                if self.seq[pos] == ami and len(self.seq) > pos:
                    self.ptm[pos] = ty
                    mylog.info("pos right")
                else:
                    mylog.error("pos wrong.")

    def run(self):
        # 二级结构
        self.readss()
        # disorder信息
        self.readdis()
        # 保守率
        self.caler()
        # anchor
        self.readanchor()
        # PTM
        self.getptm()

        # self.getmap()
        result = ["#Residue;Structure;Conserved-Rate;Anchor-Rate;PTM type;\n"]
        assert len(self.seq) == len(self.ss)
        assert len(self.seq) == len(self.er) 
        try:
            assert len(self.seq) == len(self.anchor)
        except Exception as e:
            mylog.error("len(self.seq):%s, len(self.anchor):%s" %(len(self.seq), len(self.anchor)))
        for X, ss, er, an, ptm in zip(self.seq, self.ss, self.er, self.anchor, self.ptm):
            result.append("{:s} {:s} {:.3f} {:.3f} {:s}\n".format(X, ss, er, an, ptm))
        fw = open(os.path.join(RET_DIR, self.idd+".ret"), 'w')
        fw.writelines(result)
        fw.close()
        
        #fp = open(os.path.join(RET_DIR, self.idd+".json"), 'w')
        #json.dump(self.result, fp, indent=4)
        
def run(filename):
    idd = os.path.basename(filename)[:-4]
    mylog.info(idd)
    s = SS(idd)
    if s.tag == 1:
        s.run()
    else:
        mylog.warn("%s not run." %(idd))

if __name__ == '__main__':        
    import threadpool
    pool = threadpool.ThreadPool(1)
    requests = threadpool.makeRequests(run, glob(os.path.join(SS_DIR, "*.mat")))
    [pool.putRequest(req) for req in requests]
    pool.wait()
