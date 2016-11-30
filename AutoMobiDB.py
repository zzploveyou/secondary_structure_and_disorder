# coding:utf-8
'''
get secondary structure and disorder info
using MobiDB
'''
import json
import os
import subprocess
import sys
import time
import urllib
from optparse import OptionParser

from Bio import SeqIO
from Bio.Seq import Seq

import log


class AutoMobiDB:

    def __init__(self, log, fastaName):
        """初始化构造函数"""
        self.log = log
        if not os.path.exists(fastaName):
            self.log.fatal("%s is not existed." % fastaName)
            sys.exit(1)
        self.fastaName = os.path.realpath(fastaName)
        self.log.info("found fastafile: %s" % os.path.basename(self.fastaName))
        self.path = os.path.dirname(self.fastaName)
        self.path_structure = os.path.join(self.path, "structure")
        if not os.path.exists(self.path_structure):
            self.log.info("make structure dir.")
            os.makedirs(self.path_structure)
        self.IDs = []

    def getIDs(self, parseType='fasta'):
        """从fasta格式文件中获取需要的UniProt ID"""
        for seq_record in SeqIO.parse(self.fastaName, parseType):
            if len(seq_record.seq) >= 30 and 'X' not in seq_record.seq:
                ID = seq_record.id.split("|")[1]
                self.log.debug("ID: %s" % ID)
                self.IDs.append(ID)
        self.log.info("total numbers: %d" % (len(self.IDs)))

    def getData(self, url):
        """获取网页数据"""
        response = urllib.urlopen(url)
        time.sleep(1)
        data = json.loads(response.read())
        return data

    def getdis(self, ID, cool=False):
        """获取disorder信息"""
        url = "http://mobidb.bio.unipd.it/ws/entries/%s/disorder" % (ID)
        data = self.getData(url)
        if cool:
            return json.dumps(data, sort_keys=False, indent=4)
        else:
            return data

    def getuni(self, ID, cool=False):
        """获取UniProt信息(包括二级结构)"""
        url = "http://mobidb.bio.unipd.it/ws/entries/%s/uniprot" % (ID)
        data = self.getData(url)
        if cool:
            return json.dumps(data, sort_keys=False, indent=4)
        else:
            return data

    def download(self, ID):
        """对单个ID做处理, 如果不存在或者文件大小为0才会下载, 否则不做处理"""
        disfile = os.path.join(self.path_structure, ID + '_dis.json')
        unifile = os.path.join(self.path_structure, ID + '_uni.json')
        if not os.path.exists(disfile) or os.path.getsize(disfile) == 0:
            result = self.getdis(ID)
            with open(disfile, 'w') as f:
                f.write(str(result))
            self.log.debug("disorder json  %s has been downloaded." %
                           os.path.basename(disfile))

        if not os.path.exists(unifile) or os.path.getsize(unifile) == 0:
            result = self.getuni(ID)
            with open(unifile, 'w') as f:
                f.write(str(result))
            self.log.debug("uniprot json  %s has been downloaded." %
                           os.path.basename(unifile))

    def downloads(self):
        import threading
        th = []
        for ID in self.IDs:
            self.log.debug("work ID %s start." % ID)
            # 多线程添加下载条目
            if len(th) != 8:
                th.append(threading.Thread(target=self.download, args=(ID,)))
            else:
                for t in th:
                    t.start()
                for t in th:
                    t.join()
                th = []

# json.load(open("filename.json"))
# json.get('consensus').get('full')


def main(fastafile):
    start = time.time()

    mylog = log.Terminal_log(brief=True)
    amd = AutoMobiDB(mylog, fastafile)
    amd.getIDs()

    # 从MobiDB下载fasta中各条序列的disorder, uniprot的json文件
    amd.downloads()

    end = time.time()
    mylog.info("Total Time: %.2f s" % (end - start))
    mylog.info("Average Time: %.2f s" % ((end - start) / len(amd.IDs)))
    mylog.done("done.")

if __name__ == '__main__':
    # filename = "/home/zzp/DATABASE/Human/test.fasta"
    filename = "/home/zzp/DATABASE/Human/uniprot-human.fasta"
    main(filename)