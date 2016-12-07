import os
from glob import glob
import log

SEQ_DIR = "/home/biolab/zzp/Human/sequences/"
ANCHOR_DIR = "/home/biolab/zzp/Human/ANCHOR/"
mylog = log.Terminal_log()

assert os.path.exists(SEQ_DIR)

num = 0
for filename in glob(os.path.join(SEQ_DIR, "*.fasta")):
    num += 1
    mylog.info(filename)
    idd = os.path.basename(filename)[:-6]
    mylog.info(idd)
    outfile = os.path.join(ANCHOR_DIR, idd+".out")
    mylog.debug(outfile)
    if not os.path.exists(outfile):
        order = "perl ~/zzp/program/ANCHOR/anchor.pl %s -m ~/zzp/program/ANCHOR/motifs_elm.txt > %s" %(filename, outfile)
        os.system(order)
    mylog.done("job %d: %s done." % (num, idd))    
