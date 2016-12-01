# coding:utf-8
'''
get secondary structure and disorder info
using d2p2
'''
import json
import os
import subprocess
import sys
import time
import urllib
import urllib2
from optparse import OptionParser

from Bio import SeqIO
from Bio.Seq import Seq

import log


class Autod2p2:

    def __init__(self, log, fastaName):
        """初始化构造函数"""
        self.log = log
        if not os.path.exists(fastaName):
            self.log.fatal("%s is not existed." % fastaName)
            sys.exit(1)
        self.fastaName = os.path.realpath(fastaName)
        self.log.info("found fastafile: %s" % os.path.basename(self.fastaName))
        self.path = os.path.dirname(self.fastaName)
        self.path_structure = os.path.join(self.path, "disorder_d2p2")
        if not os.path.exists(self.path_structure):
            self.log.info("make disorder dir.")
            os.makedirs(self.path_structure)
        self.num_download = 0
        self.IDs = []

    def getIDs(self, parseType='fasta'):
        """从fasta格式文件中获取需要的UniProt ID"""
        for seq_record in SeqIO.parse(self.fastaName, parseType):
            if len(seq_record.seq) >= 30 and 'X' not in seq_record.seq:
                ID = seq_record.id.split("|")[1]
                self.log.debug("ID: %s" % ID)
                self.IDs.append(ID)
        self.log.info("total numbers: %d" % (len(self.IDs)))

    def getdis(self, ID):
        """获取disorder信息"""
        data = 'seqids=["%s"]' %ID
        request = urllib2.Request('http://d2p2.pro/api/seqid', data)
        response = json.loads(urllib2.urlopen(request).read())
        if response[ID] == []:
            return ""
        else:
            result = response[ID][0][2]['disorder']['consensus']
            result = [str(i) for i in result]
            return " ".join(result)


    def download(self, ID):
        """对单个ID做处理, 如果不存在或者文件大小为0才会下载, 否则不做处理"""
        disfile = os.path.join(self.path_structure, ID + '.txt')
        if not os.path.exists(disfile) or os.path.getsize(disfile) == 0:
            self.num_download += 1
            result = str(self.getdis(ID))
            if result != "":
                with open(disfile, 'w') as f:
                    f.write(result)
                return ID
            else:
                self.log.warn("cannot fetching ID: %s." %ID)
        return None

    def callback(self, requests, result):
        """callback of download"""
        if result != None:
            self.log.done("%s done, download num: %d" %(result, self.num_download))

    def download_all(self):
        """download all"""
        import threadpool
        pool = threadpool.ThreadPool(8)
        requests = threadpool.makeRequests(self.download, self.IDs, self.callback)
        for req in requests: pool.putRequest(req)
        pool.wait()
        return True

    def downloads(self):
        """auto retry"""
        tag = False
        while tag != True:
            try:
                tag = self.download_all()
            except Exception as e:
                pass

# json.load(open("filename.json"))
# json.get('consensus').get('full')


def main(fastafile):
    start = time.time()

    mylog = log.Terminal_log(brief=True)
    amd = Autod2p2(mylog, fastafile)
    amd.getIDs()

    # 从d2p2下载fasta中各条序列的disorder, uniprot的json文件
    amd.downloads()

    end = time.time()
    mylog.info("Total Time: %.2f s" % (end - start))
    mylog.info("Average Time for each download: %.2f s" % ((end - start) / amd.num_download))
    mylog.done("done.")

if __name__ == '__main__':
    #filename = "/home/zzp/DATABASE/Human/test.fasta"
    filename = "/home/zzp/DATABASE/Human/UP000005640_9606.fasta"
    main(filename)
