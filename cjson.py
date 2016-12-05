# coding:utf-8
"""
修复json文件错误
"""
from glob import glob
import os
PATH = "/home/biolab/zzp/Human/disorder/"
for filename in glob(os.path.join(PATH, "*.json")):
    print filename
    new = open(filename).read().replace("u\'","\'").replace("\'", "\"")
    fw = open(filename, 'w')
    fw.write(new)
    fw.close()

