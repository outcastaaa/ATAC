# samtools
## 参数  

```bash
$ samtools

Program: samtools (Tools for alignments in the SAM format)
Version: 1.16 (using htslib 1.16)

Usage:   samtools <command> [options]

Commands:
  -- Indexing
     dict           create a sequence dictionary file
     faidx          index/extract FASTA
     fqidx          index/extract FASTQ
     index          index alignment

  -- Editing
     calmd          recalculate MD/NM tags and '=' bases
     fixmate        fix mate information
     reheader       replace BAM header
     targetcut      cut fosmid regions (for fosmid pool only)
     addreplacerg   adds or replaces RG tags
     markdup        mark duplicates
     ampliconclip   clip oligos from the end of reads

  -- File operations
     collate        shuffle and group alignments by name
     cat            concatenate BAMs
     consensus      produce a consensus Pileup/FASTA/FASTQ
     merge          merge sorted alignments
     mpileup        multi-way pileup
     sort           sort alignment file
     split          splits a file by read group
     quickcheck     quickly check if SAM/BAM/CRAM file appears intact
     fastq          converts a BAM to a FASTQ
     fasta          converts a BAM to a FASTA
     import         Converts FASTA or FASTQ files to SAM/BAM/CRAM
     reference      Generates a reference from aligned data

  -- Statistics
     bedcov         read depth per BED region
     coverage       alignment depth and percent coverage
     depth          compute the depth
     flagstat       simple stats
     idxstats       BAM index stats
     phase          phase heterozygotes
     stats          generate stats (former bamcheck)
     ampliconstats  generate amplicon specific stats

  -- Viewing
     flags          explain BAM flags
     head           header viewer
     tview          text alignment viewer
     view           SAM<->BAM<->CRAM conversion
     depad          convert padded BAM to unpadded BAM
     samples        list the samples in a set of SAM/BAM/CRAM files

  -- Misc
     help [cmd]     display this help message or help for [cmd]
     version        detailed version information
```

## 常用参数  

```bash
1. samtools sort [options] input.bam  

options:
-n : 根据read的name进行排序，默认对最左侧坐标进行排序
-o : 设置排序后输出文件的文件名
-O : 后跟sam或bam，规定排序后输出文件的格式，默认是bam
-@ : 后跟正整数，指定分析所用线程数
-m : 后跟数字 + K/M/G，指定每个线程所使用内存量  


举例：对bam文件进行排序

samtools sort -o output.bam input.bam : 按默认方式排序并输出为bam
samtools sort input.bam > output.bam : 按默认方式排序并输出为bam

samtools sort -n -o output.bam input.bam : 按read name进行排序并输出为bam
samtools sort -n input.bam > output.bam : 按read name进行排序并输出为bam

链接：https://www.jianshu.com/p/e6bd632627a3


2. samtools index  

 为了能够快速访问bam文件，可以为已经基于坐标排序后bam或者cram的文件创建索引，生成以.bai或者.crai为后缀的索引文件。
 必须使用排序后的文件，否则可能会报错。
 另外，不能对sam文件使用此命令。如果想对sam文件建立索引，那么可以使用tabix命令创建。

index命令格式如下：

      samtools index [-bc] [-m INT] aln.bam |aln.cram [out.index]

      参数：

      -b 创建bai索引文件，未指定输出格式时，此参数为默认参数；

      -c 创建csi索引文件，默认情况下，索引的最小间隔值为2^14，与bai格式一致；

      -m INT 创建csi索引文件，最小间隔值2^INT；

原文链接：https://blog.csdn.net/u013553061/article/details/53192807
```
# bam文件
## 查看BAM文件内容：
使用samtools view查看BAM文件    
```bash
$ samtools view in.bam

* 如果不想从头开始看，希望快速地跳转到基因组的其它位置上，比如chr22染色体，
那么可以先用samtools index生成BAM文件的索引（如果已经有索引文件则不需该步骤），然后这样操作：  

$ samtools index in.bam  # 生成in.bam的索引文件in.bam.bai
$ samtools view in.bam chr22            # 跳转到chr22染色体
$ samtools view in.bam chr22:16050103   # 跳转到chr22:16050103位置
$ samtools view in.bam chr22:16050103-16050103  # 只查看该位置

* 这里发现，原始的.bam文件，和.sort.bam以及.name.sort.bam文件的大小不一致，并且.sort.bam小很多，检查三个文件的行数：

samtools view -c input.sort.bam
```

## 用代码查看BAM 文件内容   

```bash
xuruizhi@DESKTOP-HI65AUV:/mnt/d/project/rat6.0/output/align$ samtools view -h SRR2190795.fastq.gz.sort.bam | head -n 30

@HD     VN:1.0  SO:coordinate    # HD是必须的标准文件头  
# VN 格式版本  SO表示比对排序类型；有unknown（default），unsorted，queryname和coordinate几种


@SQ     SN:1    LN:260522016     # SQ参考序列染色体信息，顺序必须和参考序列一致，这些参考序列决定了比对结果sort的顺序


@PG     ID:hisat2       PN:hisat2       VN:2.2.1        CL:"/home/linuxbrew/.linuxbrew/bin/../Cellar/hisat2/2.2.1/bin/hisat2-align-s --wrapper basic-0 -t -x ../../genome/index/mRatBN7.2.chr1 -S ../align/SRR2190795.fastq.gz.sam --read-lengths 100,99,63,98,57,96,97,55,54,95,53,58,94,62,56,64,52,93,75,91,92,79,60,59,90,67,61,74,45,87,85,83,72,69,70,89,65,51,88,86,80,50,44,82,66,78,73,71,49,46,39,36,68,48,42,81,35,84,77,47,30,40,38,37,43,31,33,76,41,34,32 -U /tmp/131.unp"       
@PG     ID:samtools     PN:samtools     PP:hisat2       VN:1.16 CL:samtools sort -@ 4 SRR2190795.fastq.gz.sam            
 # PG是重要的read group信息，通常包含测序平台、测序文库、样本ID等信息  

@PG     ID:samtools.1   PN:samtools     PP:samtools     VN:1.16 CL:samtools view -h SRR2190795.fastq.gz.sort.bam


！ RECORD 每一行都是一条read比对信息 
注：CIGAR中的M，不能觉得它代表的是匹配就以为是百分百没有任何miss-match，这是不对的，多态性碱基或者单碱基错配也是用M标记！


SRR2190795.4538171  # read name    272  # flags比对信息位（272=0PE+16+256，根据表格对应比对情况）   1（参考序列名，或染色体编号）   \
20（从左边开始数，该基因在染色体上的具体位置）     1（比对质量值）       97M3S（CIGAR）   * (第二次比对)      0       0      \
TCCTTTTCAACAGAAGCAGAAGCTCATCTGAATATGCTCAAGGATGCTGACATCAACATTTAATCATCTCCTCACTCATCCAGGAAGAAGGGGAGATAAG  \
9?A@B@@A;BDDB?><ECEA?==.7FEC;DFF>IIGFF<FFDEFB9>B<DF@?0:*D?<@FDC<??C:*FBAA<A33A:?4<F>FEBFFFBC?4DDD8=:  \
AS:i:-3(匹配的得分) ZS:i:-3() XN:i:0(在参考序列上模糊碱基的个数)  XM:i:0(错配的个数)  XO:i:0(gap open的个数) \
XG:i:0(gap 延伸的个数)  NM:i:0(经过编辑的序列)  MD:Z:97(代表序列和参考序列错配的字符串)  YT:Z:UU(UU表示不是pair中一部分)  NH:i:5

SRR2190795.13114124     272     1       23      1       98M     *       0       0       TTTTCAACAGAAGCAGAAGCTCATCTGAATATGCTCAAGGATGCTGACATCAACATTTAATCATCTCCTCACTCATCCAGGAAGAAGGGGAGATCAGT      CCA@C@6.:DDDD>;?7777==7EC=7:EGEIFB4CFFB94?899*0?<DB:0B?*<C1?1:?EC?;A<229+A<;FEAF<G<BBC:DBFDD?B;;@@      AS:i:0  ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:98YT:Z:UU  NH:i:5
SRR2190795.5690499      256     1       151     1       100M    *       0       0       TCACCATGTTACAAAAATAGCAAGCTGCCATAATAAAAAATGAGGCTCCTCTATCCAGCACCAGATAGCATCATTTTACTTTCAAGCCTAGAAATTGCAC    @C@DDFEFHFHHHJIBHGEHHIJJJJJJJIIJEIIIIJJJJIJJJGGIJJJJGGIJIJJJJJGEACHHHHGHFFF;@@CCEEEECCDDDDDDDDDDDDCC    AS:i:-6 ZS:i:-6 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:41A58      YT:Z:UU NH:i:5
```

