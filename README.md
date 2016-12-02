进化速率与结构, 与功能的关系. 依赖结构更多,而非功能.

验证进化速率:

    disorder > 中间态 > coil > sheet > helix


人类蛋白组UP000005640_9606.fasta含有21032条序列,

```
/home/zzp/DATABASE/Human/UP000005640_9606.fasta
```

满足序列长度不小于30,且不含有位置残基X的序列有20761条.

```
/home/zzp/DATABASE/Human/sequences
```

编程实现了根据uniprot id从MobiDB, d2p2抓取disorder信息(8个线程的线程池)

|数据库|数量|单条时间|格式|
|:-:|:-:|:-:|:-:|
|MobiDB | 20169 |  0.26s |  (格式为json格式)|
|d2p2|    18325 |  0.11s |  (格式为每条序列各个残基位置为disorder的打分,从0~9)|

PSIPRED计算二级结构信息.(更改runpsipred脚本中结果保存目录.)
