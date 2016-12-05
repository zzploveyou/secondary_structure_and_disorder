# coding:utf-8
"""
For what:
    ls job calculated now, sort by modificated time.
    find process that cannot get result.
Usage:
    start this script, if some process mtime don't change,
    then you can kill the process.

"""
from Bio import SeqIO
import os
from glob import glob
import datetime
import time
import sys

FORMAT = '%H:%M:%S'
PFORMAT = '  {:0>2s} {:^8s} {:^12s} {:^10s} {:^20s}'

def htime(timestamp):
    tt = datetime.datetime.fromtimestamp(timestamp)
    return tt.strftime(FORMAT)

def jobnow(dirs="."):
    print
    print PFORMAT.format('id', 'mtime', 'uni-id', 'seqlen', 'blastfile')
    idds = []
    for filename in glob(os.path.join(dirs, "*.fasta")):
        for sec in SeqIO.parse(filename, 'fasta'):
            idd = sec.id.split("|")[1]
            length = len(sec.seq)
        blastfile = os.path.splitext(filename)[0] + ".blast"
        idds.append((os.path.getmtime(blastfile), idd, length, filename))
    idds.sort(reverse=True)
    idx = 1
    for mtime, idd, length, filename in idds:
        print PFORMAT.format(str(idx), htime(mtime), idd, str(length), filename)
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
