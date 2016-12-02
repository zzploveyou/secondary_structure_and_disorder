# coding:utf-8
"""
threadpool psipred cal ss.
"""
import os
import threadpool
from glob import glob
import log

mylog = log.Terminal_log()

def psipred(fastafile):
    idd = os.path.basename(fastafile).split('.')[0]
    order = " ".join(('psipred', fastafile))
    os.system(order)
    return idd

def callback(requests, result):
    if result == True:
        mylog.done("%s done." % result)
    else:
        mylog.error("error.")

def cal(seqdir):
    filenames = []
    for filename in glob("*.fasta"):
        f1 = os.path.splitext(filename)+".ss2"
        f2 = os.path.splitext(filename)+".horiz"
        if not os.path.exists(f1) or not os.path.exists(f2):
            filenames.append(filename)

    pool = threadpool.ThreadPool(8)
    requests = threadpool.makeRequests(psipred, filenames, callback)
    for req in requests: pool.putRequest(req)
    pool.wait()

if __name__ == "__main__":
    seqdir = "/home/zzp/DATABASE/Human/sequences/"
    cal(seqdir)
