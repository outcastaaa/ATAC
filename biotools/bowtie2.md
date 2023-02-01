# 建立索引
## 用法
```bash
bowtie2-build [options]* <reference_in> <bt2_base>
```
<reference_in>：如果此处使用-f 参数，则指明index的参考fasta 文件；如果使用-c参数，则指明index的参考序列，例如，GGTCATCCT,ACGGGTCGT,CCGTTCTATGCGGCTTA.  
<bt2_base>：指的是生成的index文件的前缀，默认情况，bowtie2-build产生NAME.1.bt2, NAME.2.bt2, NAME.3.bt2, NAME.4.bt2, NAME.rev.1.bt2, and NAME.rev.2.bt2, where NAME is <bt2_base>.  
--threads 使用的线程数  

## 举例
```bash
bowtie2-build -f /public/Reference/GRCh38.primary_assembly.genome.fa --threads 24 GRCh38
```

# 比对
[参考文章1](https://www.jianshu.com/p/f84ffba2ec1e)    [参考文章2](https://cloud.tencent.com/developer/article/1772432)  

## 单端比对
```bash
bowtie2 [options]* -x <bt2-idx> -U <fq> -S <sam_output> -p <threads> 2>Align.summary

-x：参考基因组index文件的前缀（包括路径）
-U：单端测序的fastq文件
-S：输出的SAM文件，包含比对结果
-p：使用的线程数
"2>Align.summary"：将输出到屏幕的标准误(standard error)重导向到"Align.summary"文件

# 其格式通常如下
## Single-end
20000 reads; of these:
  20000 (100.00%) were unpaired; of these:
    1247 (6.24%) aligned 0 times
    18739 (93.69%) aligned exactly 1 time
    14 (0.07%) aligned >1 times
93.77% overall alignment rate

## Paired-end
10000 reads; of these:
  10000 (100.00%) were paired; of these:
    650 (6.50%) aligned concordantly 0 times
    8823 (88.23%) aligned concordantly exactly 1 time
    527 (5.27%) aligned concordantly >1 times
    ----
    650 pairs aligned concordantly 0 times; of these:
      34 (5.23%) aligned discordantly 1 time
    ----
    616 pairs aligned 0 times concordantly or discordantly; of these:
      1232 mates make up the pairs; of these:
        660 (53.57%) aligned 0 times
        571 (46.35%) aligned exactly 1 time
        1 (0.08%) aligned >1 times
96.70% overall alignment rate
The indentation indicates how subtotals relate to t
```

## 双端比对
```bash
bowtie2 [options]* -x <bt2-idx> -1 <fq1> -2 <fq2> -S <sam_output> -p <threads> 2>Align.summary

双端比对模式基本与单端一致，只需替换fastq文件传入的参数即可
-1：一链fastq文件  
-2：二链fastq文件
```
