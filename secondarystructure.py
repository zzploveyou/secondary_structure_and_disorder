# coding:utf-8
"""
threadpool psipred cal ss.
"""
import os
import multiprocessing
from Bio import SeqIO
import time
from glob import glob
import log
import sys

mylog = log.Terminal_log(brief=False)
PSIPRED = '/home/biolab/zzp/program/PSIPRED/4.01/runpsipred'
RESULT_DIR = '/home/biolab/zzp/Human/secondary_structure/'

def psipred(fastafile, seqlen):
    idd = os.path.basename(fastafile).split('.')[0]
    order = " ".join((PSIPRED, fastafile))
    mylog.info("add JOB: %s, length: %d" % (os.path.basename(fastafile), seqlen))
    os.system(order)
    return idd

def seq_len(filename):
    for rec in SeqIO.parse(filename, "fasta"):
        return len(rec.seq)

def cal(seqdir, process):

    start = time.time()
    pool = multiprocessing.Pool(processes=process)
    results = []

    filenames = []
    for filename in glob(os.path.join(seqdir, "*.fasta")):
        idd = os.path.basename(filename)[:-6]
        f1 = os.path.join(RESULT_DIR, idd + ".ss2")
        f2 = os.path.join(RESULT_DIR, idd + ".horiz")
        """not predicted."""
        if not os.path.exists(f1) or not os.path.exists(f2):
            filenames.append((seq_len(filename), filename))
        else:
            mylog.debug("%s already predicted." %filename)
    """short length sequence first"""
    filenames.sort()
    for seqlen, filename in filenames:
        results.append(pool.apply_async(psipred, (filename, seqlen,)))
    
    num = 0
    for r in results:
        num += 1
        t = time.time() - start
        mylog.done("num: %d, %s done, already %.2f min, ave time: %.2f min." % (num, r.get(), t/60., t*1.0/num/60.))
        sys.stdout.flush()

    pool.close()
    pool.join()
    mylog.done("done.")
    # 注意psipred结果位置.
    order = " ".join(("mv", "*.ss", "*.ss2", "*.horiz", "%s" % RESULT_DIR))
    os.system(order)

    
if __name__ == "__main__":
    seqdir = "/home/biolab/zzp/Human/sequences/"
    """make sure assert below in case of errors."""
    assert os.path.realpath(os.path.dirname(__file__)) == os.getcwd()
    cal(seqdir, process=24)
