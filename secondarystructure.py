# coding:utf-8
"""
threadpool psipred cal ss.
"""
import os
import multiprocessing
from glob import glob
import log

mylog = log.Terminal_log()
PSIPRED = '/home/zzp/program/release/PSIPRED/4.01/runpsipred'
def psipred(fastafile):
    idd = os.path.basename(fastafile).split('.')[0]
    order = " ".join((PSIPRED, fastafile))
    mylog.debug("order: %s" % order)
    os.system(order)
    return idd

def cal(seqdir):
    pool = multiprocessing.Pool(processes=4)
    results = []
    for filename in glob(os.path.join(seqdir, "*.fasta")):
        f1 = os.path.splitext(filename)[0] + ".ss2"
        f2 = os.path.splitext(filename)[0] + ".horiz"
        mylog.debug("f1: %s, f2: %s" % (f1, f2))
        if not os.path.exists(f1) or not os.path.exists(f2):
            results.append(pool.apply_async(psipred, (filename,)))
    for r in results:
        mylog.info("%s done." % r.get())
    pool.close()
    pool.join()
    mylog.done("done.")

    
if __name__ == "__main__":
    seqdir = "/home/zzp/DATABASE/Human/sequences/test/"
    cal(seqdir)
