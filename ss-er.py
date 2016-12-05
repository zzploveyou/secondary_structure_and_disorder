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
reload(sys)
sys.setdefaultencoding('utf-8')
mylog = log.Terminal_log()
PATH = "/home/biolab/zzp/Human"
SS_DIR = os.path.join(PATH, "secondary_structure")
DIS_DIR = os.path.join(PATH, "disorder")
ER_DIR = os.path.join(PATH, "ER")
RET_DIR = os.path.join(PATH, "SS-ER")

class SS:

    def __init__(self, idd):
        self.idd = idd
        self.ssfile = os.path.join(SS_DIR, idd + ".ss2")
        self.disfile = os.path.join(DIS_DIR, idd + ".json")
        self.ss = []
        self.tag = 1
        if not os.path.exists(self.ssfile):
            mylog.debug("ssfile not existed.")
            self.tag = 0
        if not os.path.exists(self.disfile):
            mylog.debug("disfile not existed.")
            self.tag = 0

    def readss(self):
        for line in open(self.ssfile):
            tmp = line.strip().split()
            if len(tmp) == 6:
                self.ss.append(tmp[2])

    def readdis(self):
        #self.disfile = "123.json"
        data = open(self.disfile).read().encode('utf-8')
        # print data
        js = json.loads(data)
        print js

        

    def getss(self):
        return self.ss       
        
    def result(self):
        if self.tag == 0:
            return None
        else:
            return 123

s = SS("P59796")
s.readss()
s.readdis()
#print s.getss()

