# 目标

进化速率与结构, 与功能的关系. 依赖结构更多,而非功能.

验证进化速率:

    disorder > 中间态 > coil > sheet > helix


# 过程

## 从uniprot下载人类蛋白

人类蛋白组UP000005640_9606.fasta含有21032条序列,

```
/home/zzp/DATABASE/Human/UP000005640_9606.fasta
```

## 预处理

满足序列长度不小于30,且不含有位置残基X的序列有20761条.

```
# splitfasta.py 将 UP000005640_9606.fasta分解到如下文件夹:
/home/zzp/DATABASE/Human/sequences
```

## 抓取disorder信息

编程实现了根据uniprot id从MobiDB, d2p2抓取disorder信息(8个线程的线程池)

```
# MobiDB data
/home/zzp/DATABASE/Human/disorder/
# d2p2 data
/home/zzp/DATABASE/Human/disorder_d2p2/
```

|数据库|数量|单条时间|格式|
|:-:|:-:|:-:|:-:|
|MobiDB | 20169 |  0.26s |  (格式为json格式)|
|d2p2|    18325 |  0.11s |  (格式为每条序列各个残基位置为disorder的打分,从0~9)|

## 计算二级结构

PSIPRED计算二级结构信息.(更改runpsipred脚本中结果保存目录.)
biolab@222.30.48.144
passwd: qaztmko*330
```
任务: 320长度序列PSIPRED预测二级结构 (W5XKT8.fasta)

开1个进程跑一次任务,每个任务分配8个线程, 平均一个用时8分钟
开1个进程跑一次任务,每个任务分配1个线程, 平均一个用时11分钟

开2个进程同时跑两次任务,每个任务分配1个线程, 平均一个用时6分钟
开2个进程同时跑两次任务,每个任务分配2个线程, 平均一个用时3分钟
开2个进程同时跑两次任务,每个任务分配4个线程, 平均一个用时2.2分钟

mathbio上运行:
python 进程池中 开16进程, blastpgp运行参数为默认一个进程.
对sequence length进行排序, 先对长度短的预测结构, 平均0.35min一条序列.



```
