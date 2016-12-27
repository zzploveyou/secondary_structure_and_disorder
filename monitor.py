# coding:utf-8

"""
For what:
    ls job calculated now, sort by modificated time.
    find process that cannot get result.
Usage:
    start this script, if some process mtime don't change,
    then you can kill the process.
PER:
    percent progress, search, round1, round2
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
PFORMAT = ' {:0>2s} {:^7s} {:^11s} {:^11s} {:>4s} {:^30s} {:^11s}'

LENs = {}

def htime(timestamp):
    tt = datetime.datetime.fromtimestamp(timestamp)
    return tt.strftime(FORMAT)

def jobnow(dirs="."):
    global LENs
    print
    print PFORMAT.format('Id', 'Wait', 'Modified', 'UniID', 'len', 'blastfile', 'Progress(50)')
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
        f = open(blastfile).read()
        n = f.count("Searching")
        per = f.strip().split("\n")[-1].count(".")
        per = str('{:02d}'.format(per))
        
        if 'Searching' not in f.strip().split("\n")[-1]:
            """the last line confused"""
            per = "pre"
            n += 1
        
        if n == 3:
            per = "\033[31mstep3 : \033[0m" + per
        elif n == 2:
            per = "\033[33mstep2 : \033[0m" + per
        elif n == 1:
            per = "\033[32mstep1 : \033[0m" + per
        else:
            per = "0"
        # working time: wtime
        idds.append((time.time()-os.path.getctime(filename), per,
            os.path.getmtime(blastfile), idd, length, filename))
    idds.sort(reverse=False)
    idx = 1
    for wtime, per, mtime, idd, length, filename in idds:
        print PFORMAT.format('%02d' % idx, '{:05.1f}'.format(wtime*1.0/60), htime(mtime), idd, str(length), filename, per)
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
