# coding:utf-8
"""
threadpool psipred cal ss.
"""
import os
import multiprocessing
import time
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

    start = time.time()
    pool = multiprocessing.Pool(processes=2)
    results = []
    '''
    for filename in glob(os.path.join(seqdir, "*.fasta")):
        f1 = os.path.splitext(filename)[0] + ".ss2"
        f2 = os.path.splitext(filename)[0] + ".horiz"
        # mylog.debug("f1: %s, f2: %s" % (f1, f2))
        if not os.path.exists(f1) or not os.path.exists(f2):
            results.append(pool.apply_async(psipred, (filename,)))
    '''
    filename = os.path.join(seqdir, "W5XKT8.fasta")
    results.append(pool.apply_async(psipred,(filename,)))
    filename = os.path.join(seqdir, "W5XKT8.fasta")
    results.append(pool.apply_async(psipred,(filename,)))
    num = 0
    for r in results:
        num += 1
        t = time.time()-start
        mylog.info("num: %d, %s done, already %.2f min, ave time: %.2f min." % (num, r.get(), t/60., t*1.0/num/60.))

    pool.close()
    pool.join()
    mylog.done("done.")
    # 注意psipred结果位置.

    
if __name__ == "__main__":
    seqdir = "/home/zzp/DATABASE/Human/sequences/"
    cal(seqdir)
