# coding:utf-8
"""
For what:
    ls job calculated now, sort by modificated time.
    find process that cannot get result.
Usage:
    start this script, if some process mtime don't change,
    then you can kill the process.
WT:
    wait time.
MT:
    modificated time.
"""
from Bio import SeqIO
import os
from glob import glob
import datetime
import time
import sys

FORMAT = '%H:%M:%S'
PFORMAT = '  {:0>2s} {:05.1f} {:^11s} {:^12s} {:^11s} {:^20s}'

LENs = {}

def htime(timestamp):
    tt = datetime.datetime.fromtimestamp(timestamp)
    return tt.strftime(FORMAT)

def jobnow(dirs="."):
    global LENs
    print
    print '  {:0>2s} {:^6s} {:^10s} {:^12s} {:^10s} {:^20s}'.format('id', 'WT', 'MT', 'uni-id', 'seqlen', 'blastfile')
    idds = []
    for filename in glob(os.path.join(dirs, "psitmp*.fasta")):
        if filename not in LENs:
            for sec in SeqIO.parse(filename, 'fasta'):
                idd = sec.id.split("|")[1]
                length = len(sec.seq)
                LENs[filename] = (idd, length)
        idd = LENs[filename][0]
        length = LENs[filename][1]
        blastfile = os.path.splitext(filename)[0] + ".blast"
        # working time: wtime
        idds.append((time.time()-os.path.getctime(filename),
            os.path.getmtime(blastfile), idd, length, filename))
    idds.sort(reverse=False)
    idx = 1
    for wtime, mtime, idd, length, filename in idds:
        print PFORMAT.format(str(idx), wtime*1.0/60, htime(mtime), idd, str(length), filename)
        idx += 1
    print

def main(refresh):
    while True:
        os.system('clear')
        print 'now: %s, refresh time: %s s' % (datetime.datetime.now().strftime(FORMAT), refresh)
        jobnow()
        time.sleep(refresh)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print
        print 'now: %s' % datetime.datetime.now().strftime(FORMAT)
        jobnow()
    elif len(sys.argv) == 2:
        try:
            refresh = int(sys.argv[1])
            main(refresh)
        except Exception as e:
            print e
