
# ATAC-seq分析

- [0.Introduction](#0.Introduction)
- [1.Prepare](#1.Prepare)
- [2.Biotools](#2.Biotools)
    - [2.0 management](#20-management)
    - [2.1 sratoolkit](#21-sratoolkit)
    - [2.2 fastqc](#22-fastqc)
    - [2.3 multiqc](#23-multiqc)
    - [2.4 TrimGalore](#24-TrimGalore)



    - [2.6 hisat2](#26-hisat2)
    - [sortmerna](#sortmerna)
    - [2.7 samtools](#27-samtools)
    - [2.8 HTseq](#28-htseq)
    - [2.9 R](#29-r)
    - [2.10 Rstudio](#210-rstudio)
    - [2.11 parallel](#211-parallel)
    - [StringTie[可选]](#stringtie可选)

- [3.Data](#3.Data)
    - [3.1 sratoolkit](#31-sratoolkit)
    - [3.2 fastqc](#32-fastqc)



	- [Sequencing data](#sequencing-data)
	- [Reference genome data](#reference-genome-data)
- [Quality control and trimming](#quality-control-and-trimming)
- [Methylation analysis](#methylation-analysis)
	- [Genome indexing](#genome-indexing)
	- [Read alignment](#read-alignment)
	- [Aligned reads deduplication](#aligned-reads-deduplication)
	- [Methylation information extracting](#methylation-information-extracting)
- [Downstream analysis](#downstream-analysis)
	- [Input data preparation](#input-data-preparation)
	- [DML/DMR detection](#dmldmr-detection)
- [Practical methylation information analysis](#practical-methylation-information-analysis)
- [Reference](#reference)
- [Author](#author)





# 0.Introduction  

ATAC-seq（Assay for Transposase-Accessible Chromatin with high throughput sequencing） 是2013年由斯坦福大学William J. Greenleaf和Howard Y. Chang实验室开发的用于研究染色质可及性（通常也理解为染色质的开放性）的方法，原理是通过转座酶Tn5容易结合在开放染色质的特性，然后对Tn5酶捕获到的DNA序列进行测序。  


ATAC-seq利用DNA转座酶技术实现染色质可及性分析。DNA转座酶可以将自身结合的一段序列随机插入到基因组中。在ATAC-seq试验中，细胞或组织样本在核质分离后，将细胞核单独收集在一起，并通过转座酶Tn5对核内的染色质进行打断。紧密包裹的染色质DNA不会受到转座酶的打断，而开放区域的染色质DNA会被转座酶随机插入并打断。将这些打断后的DNA收集在一起，进行后续的建库、测序、分析，即可得到开放染色质的信息。ATAC-seq中的peak，往往是启动子、增强子序列，以及一些调控因子结合的位点。  

ATAC-seq可用于：  

- 生成表观基因组图谱  

- 得到在不同组织或不同条件下对应可及性区域    

- 得到核小体位置  

- 鉴定重要转录因子  

- 生成转录因子结合区域的特征(footprinting)  

[具体看该文章](https://github.com/outcastaaa/ATAC/blob/main/review%20of%20ATAC-seq.md)  


## 数据分析具体流程：    
预处理（Pre-analysis）包括比对前的质量控制 QC（Pre-alignment QC）、比对（Alignment）、比对后处理（Post alignment processing）、QC。  
核心分析（Core analysis）包括 Peak calling。  
高级分析（Advance analysis）包括 Peak、motif、footprint、nucleosome 分析。  
多组学整合包括与 ChIP-seq、RNA-seq 数据的整合以及调控网络的重建。    

[数据分析详细流程](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6996192/figure/Fig2/)   





# 1.Prepare  
```bash
#将目录建在d盘 
cd /mnt/d 
# 建立目录
mkdir biosoft  
mkdir ATAC    
cd ./ATAC
mkdir genome sequence qc align clean motif peaks
#查看结构
@xxx:/mnt/d/ATAC$ tree
```

# 2.Biotools 
软件详细用法记录在github的[biotools](https://github.com/outcastaaa/ATAC/tree/main/biotools)文件夹中。  


## 2.0 management  

Linux brew  
来源[wang-q Ubuntu -](https://github.com/wang-q/ubuntu#install-linuxbrew)

## 2.1 sratoolkit   
* 使用brew安装  
```bash
@xxx:~$ brew install sratoolkit
```
## 2.2 fastqc  
* 使用brew安装  
```bash
@xxx:~$ brew install fastqc
```
## 2.3 multiqc 
``` bash
# 使用python的安装器安装
pip install multiqc
```

## 2.4 TrimGalore  

! [作者GitHub](https://github.com/FelixKrueger/TrimGalore)已经更新至2021年7月的0.6.6版本
```bash
cd /mnt/d/biosoft

# Install Trim Galore
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o TrimGalore.tar.gz
tar xvzf TrimGalore.tar.gz

# Run Trim Galore
~/TrimGalore-0.6.6/trim_galore
```
结果：  
```
xuruizhi@DESKTOP-HI65AUV:/mnt/d/biosoft$ curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o TrimGalore.tar.gz
xuruizhi@DESKTOP-HI65AUV:/mnt/d/biosoft$ ls
TrimGalore.tar.gz  Trimmomatic-0.38  Trimmomatic-0.38.zip  hisat2-2.2.1  sortmerna-2.1  sortmerna-2.1.tar.gz  wget-log
xuruizhi@DESKTOP-HI65AUV:/mnt/d/biosoft$ tar xvzf TrimGalore.tar.gz
TrimGalore-0.6.6/
TrimGalore-0.6.6/.travis.yml
TrimGalore-0.6.6/Changelog.md
。。。
```
`！师兄的办法会得到一个单独的TrimGalore文件；作者的办法会得到包括trim_galore及其license在内的一个文件夹`

* fastp  
* trimmomatic  
```
cd /mnt/d/biosoft
# 先挂载到d盘相应文件 

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
unzip Trimmomatic-0.38.zip

cd Trimmomatic-0.38

# 导入临时环境变量
export PATH="$(pwd):$PATH"
```
## 2.6 hisat2  

1. [hisat2官网更改](https://daehwankimlab.github.io/hisat2/)
2. 右侧download下载,直接点击下载即可，不需要回到终端再下载。下载完成后剪切到d/biosoft文件夹内解压
```
Version: HISAT2 2.2.1
Release Date: 7/24/2020

Linux_x86_64	https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
```
3. 回到终端写入环境
```

# 导入临时环境变量
$ export PATH="~/biosoft/hisat2-2.1.0:$PATH"

# 测试是否可用
$ hisat2 -h

xuruizhi@DESKTOP-HI65AUV:/mnt/d/biosoft$  hisat2 -h
HISAT2 version 2.2.1 by Daehwan Kim (infphilo@gmail.com, www.ccb.jhu.edu/people/infphilo)
Usage:
```


## 2.7 samtools
最新版本为1.16  
本地下载时，在配制这步出错，使用`brew install samtools`安装

## 2.8 HTseq
```
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple HTseq
```

## 2.9 R
最新版本4.2.1_2  
先进入官网，用清华镜像源下载合适版本的R，再`brew install r`

```
xuruizhi@DESKTOP-HI65AUV:~$ brew install r
HOMEBREW_BREW_GIT_REMOTE set: using https://mirrors.tuna.tsinghua.edu.cn/git/homebrew/brew.git for Homebrew/brew Git remote.
HOMEBREW_CORE_GIT_REMOTE set: using https://mirrors.tuna.tsinghua.edu.cn/git/homebrew/homebrew-core.git for Homebrew/core Git remote.
Running `brew update --auto-update`...
Warning: r 4.2.1_2 is already installed and up-to-date.
To reinstall 4.2.1_2, run:
  brew reinstall r
  ```
* !R 安装时多次尝试，RStudio都识别不到，因此直接在官网选择`Download R for Windows; install R for the first time`下载安装包即可；注意可以将两个文件放在同一个文件夹内  
[参考](https://blog.csdn.net/m0_49354332/article/details/116059239)  


## 2.10 Rstudio
进入网站：`https://www.rstudio.com/products/rstudio/download/`  
R studio 可以在 Windows 下安装;
选择版本下载,下载完成之后双击安装。  
`
Windows 10/11	   
RStudio-2022.07.1-554.exe
`
## 2.11 parallel  
```
brew install parallel
```





# 3.Data
## 3.1 sequence data

1. 文中查找GEO数据库编号  

选取ATAC-seq测序文章 ATAC-seq Reveals an Isl1 Enhancer that Regulates Sinoatrial 
Node Development and Function (https://pubmed.ncbi.nlm.nih.gov/33044128/)，2020年12月发表于`CIRCULATION RESEARCH`期刊，文中搜索GSE，查找到数据集为`GSE148515`。GEO数据库详细介绍可在[这里](https://github.com/outcastaaa/bioinformatics-learning/blob/main/RNA-seq/RNA-SEQ%E6%B5%81%E7%A8%8B.md#32-%E6%B5%8B%E8%AF%95%E6%95%B0%E6%8D%AE%E5%AE%9E%E9%AA%8C%E6%95%B0%E6%8D%AE)查看。   


2. 登录NCBI网站查找所需数据    

登录[NCBI官网](https://www.ncbi.nlm.nih.gov/)搜索GSE编号，可看到本文所用测序技术为为RNA-seq及ATAC-seq。 点击	`SRP256236`,查看每个样本的编号；点击`Send results to Run selector`，下载所需数据。
![GSE](../ATAC/pictures/GSE.png)  

![SRA](../ATAC/pictures/SRA.png)    



3. 下载ATAC-seq数据

* 这里只下载ATAC-seq数据中的四个，`SRR11539111 SRR11539112 SRR11539115 SRR11539116`，先下载本页面的`Metadata`和`Accession List`，里面包含了数据详细信息，这两个文件下载到`/mnt/d/ATAC/sequence`。  
```bash
@xx:/mnt/d/ATAC/sequence$ ls
SRR_Acc_List.txt  SraRunTable.txt
```
* 将刚才在`Run selector`中查找到的数据的编号复制下来，之后下载测序数据，下载脚本如下，这里是采用`SRAtoolkit`工具包中的`prefetch`工具.
* 注：如果部分数据下载失败，那么再次执行下面的代码  
```bash
# 先把需要下载的文件名称写入一个单独的txt文件中
@xx:/mnt/d/ATAC/sequence$  cat >1.txt <<EOF
SRR11539111
SRR11539112
SRR11539115
SRR11539116
EOF
  
# 修改数据存储地址，我下载到了 ~/data/sra文件夹内
vdb-config --interactive
# 执行下列代码批量下载数据
@xx:/mnt/d/ATAC/sequence$  prefetch --option-file 1.txt
#或者
@xx:/mnt/d/ATAC/sequence$ cat 1.txt | while read id;do ( nohup prefetch $id & );done

# 单个文件 
@xx:/mnt/d/ATAC/sequence$  nohup prefetch SRR11539111 -o . &
```
* 数据下载到了`~/data/sra`文件夹下（可以选择使用mv命令将数据转移到ATAC项目下），`nohup.out`文件中存储下载进程。
```bash
#查看下载情况
@xx:~/data/sra$ ls -lh

```

⑤  格式转换  
下载得到`.sra`文件，使用SRAtoolkit工具包的`fastq-dump`工具，使用它来进行格式转化
```
xuruizhi@DESKTOP-HI65AUV:~/data/sra$ ls
SRR2190795.sra  SRR2240183.sra  SRR2240185.sra  SRR2240187.sra
SRR2240182.sra  SRR2240184.sra  SRR2240186.sra  SRR2240228.sra


$ parallel -j 4 "    # 用parallel多线程加快速度，并行任务数为4
    fastq-dump --split-3 --gzip {1}    # 将sra文件转化为fastq文件之后压缩为gz文件
" ::: $(ls *.sra)     # :::后接对象

# ls *.sra代表，列举出任何以.sra结尾的文件
--gzip 将转换出的fastq文件以gz格式输出，可以节省空间
--split-3 把pair-end测序分成两个文件输出，可用于双端测序转化为两个文件，本文举例为单端测序，删掉不影响
-O 输出文件夹名，不加直接放在该文件夹


# 删除sra文件
$ rm *.sra
```
结果：  
```
Academic tradition requires you to cite works you base your article on.
If you use programs that use GNU Parallel to process data for an article in a
scientific publication, please cite:

  Tange, O. (2022, July 22). GNU Parallel 20220722 ('Roe vs Wade').
  Zenodo. https://doi.org/10.5281/zenodo.6891516

This helps funding further development; AND IT WON'T COST YOU A CENT.
If you pay 10000 EUR you should feel free to use GNU Parallel without citing.

More about funding GNU Parallel and the citation notice:
https://www.gnu.org/software/parallel/parallel_design.html#citation-notice

To silence this citation notice: run 'parallel --citation' once.

Read 15107730 spots for SRR2190795.sra
Written 15107730 spots for SRR2190795.sra
Read 17622974 spots for SRR2240183.sra
Written 17622974 spots for SRR2240183.sra
Read 19779076 spots for SRR2240184.sra
Written 19779076 spots for SRR2240184.sra
Read 24510465 spots for SRR2240182.sra
Written 24510465 spots for SRR2240182.sra
Read 11837415 spots for SRR2240186.sra
Written 11837415 spots for SRR2240186.sra
Read 23017882 spots for SRR2240185.sra
Written 23017882 spots for SRR2240185.sra
Read 19519976 spots for SRR2240187.sra
Written 19519976 spots for SRR2240187.sra
Read 17296729 spots for SRR2240228.sra
Written 17296729 spots for SRR2240228.sra
```
网上找到的另一种循环语句的方法  [https://www.jianshu.com/p/bdfa8f7e5a61](https://www.jianshu.com/p/bdfa8f7e5a61)

```
#定义存放输出数据的文件夹，需要先创建这个文件夹‘fastq’
mkdir fastq
fqdir=/trainee2/Mar7/rna/project/fastq

#转换单个文件
fastq-dump --gzip --split-3 -X 25000 -O ${fqdir} SRR1039510


#批量转换，将样本名写成文件——sample.ID，echo是打印命令，while循环的意义是生成脚本
cat sample.ID | while read id
do
 echo "fastq-dump --gzip --split-3 -X 25000 -O ${fqdir} ${id}
done >sra2fq.sh
# 提交后台运行命令，脚本文件后缀为.sh，日志文件后缀为.log，运行脚本的命令为sh
nohup sh sra2fq.sh>sra2fq.log &

#查看输出的fastq的gz压缩文件，用zless命令
zless -S SRR1039510_1.fastq.gz
```




⑥ `parallel`用法补充 [parallel](https://www.jianshu.com/p/cc54a72616a1)  
```
Usage:

parallel [options] [command [arguments]] < list_of_arguments
parallel [options] [command [arguments]] (::: arguments|:::: argfile(s))...
cat ... | parallel --pipe [options] [command [arguments]]

常用选项：
::: 后面接参数
:::: 后面接文件
-j、--jobs   并行任务数
-N  每次输入的参数数量
--xargs会在一行中输入尽可能多的参数
-xapply 从每一个源获取一个参数（或文件一行）
--header  把每一行输入中的第一个值做为参数名
-m   表示每个job不重复输出“背景”（context）
-X   与-m相反，会重复输出“背景文本”
-q  保护后面的命令
--trim  lr 去除参数两头的空格，只能去除空格，换行符和tab都不能去除
--keep-order/-k   强制使输出与参数保持顺序 --keep-order/-k
--tmpdir/ --results   都是保存文件，但是后者可以有结构的保存
--delay  延迟每个任务启动时间
--halt  终止任务
--pipe    该参数使得我们可以将输入（stdin）分为多块（block）
--block  参数可以指定每块的大小
```


⑦ 格式介绍  
```
# 查看下载好的gz文件
   cd ~/data/sra
   gzip -d -c SRR2190795.fastq.gz | head -n 20

# gzip
-c或--stdout或--to-stdout 　把压缩后的文件输出到标准输出设备，不去更动原始文件。
-d或--decompress或----uncompress 　解开压缩文件。
```
结果：  

```
@SRR2190795.1 HWI-ST1147:240:C5NY7ACXX:1:1101:1320:2244 length=100
ATGCTGGGGGCATTAGCATTGGGTACTGAATTATTTTCAGTAAGAGGGAAAGAATCCATCTCCNNNNNNNNNNNNNNNNNNNNNNAAANAAAAATAAAAT
+SRR2190795.1 HWI-ST1147:240:C5NY7ACXX:1:1101:1320:2244 length=100
CCCFFFFFHHHHHJIJJJJJJJJDHHJJJIJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJHH#####################################
@SRR2190795.2 HWI-ST1147:240:C5NY7ACXX:1:1101:1598:2247 length=100
AACTTCGGTTCTCTACTAGGAGTATGCCTCATAGTACAAATCCTCACAGGCTTATTCCTAGCANNNNNNNNNNNNNNNNNNNNNNTAACAGCATTTTCAT
+SRR2190795.2 HWI-ST1147:240:C5NY7ACXX:1:1101:1598:2247 length=100
@@@7D8+@A:1CFG<C:23<:E<;FF<BHIIEHG:?:??CDF<9DCGGG?1?FEG@@<@CA#######################################
@SRR2190795.3 HWI-ST1147:240:C5NY7ACXX:1:1101:1641:2250 length=100
AGAAGGTCTTAGATCAGAAGGAGCACAGACTGGATGGTCGTGTCATTGACCCTAAAAAGGCTANNNNNNNNNNNNNNNNNNNNNTGAAGAAAATCTTTGT
+SRR2190795.3 HWI-ST1147:240:C5NY7ACXX:1:1101:1641:2250 length=100
BC@FFFDDHHHHHJJJJJJJJJJJJJJJJJJJJIJJJFHGHHEGHIIIHJIJJIJJIJIJJID#####################################
@SRR2190795.4 HWI-ST1147:240:C5NY7ACXX:1:1101:1851:2233 length=100
GGGATTTCATGGCCTCCACGTAATTATTGGCTCAACTTTCCTAATTGTCTGTCTACTACGACANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNNCNNN
+SRR2190795.4 HWI-ST1147:240:C5NY7ACXX:1:1101:1851:2233 length=100
@@?DDBDDFFDDDGHGGGGI?B;FFHGHA@FEHGHDDGHEGGFGHIGEHIIHIGGBGACD6AH#####################################
@SRR2190795.5 HWI-ST1147:240:C5NY7ACXX:1:1101:1957:2243 length=100
CAGCCATTGTGGCTCCCGATGGCTTTGACATCATTGACATGACAGCCGGAGGTCAGATAAACTNNNNNNNNNNNNNNNNNNNNNNATCNGTGGCAAAGGT
+SRR2190795.5 HWI-ST1147:240:C5NY7ACXX:1:1101:1957:2243 length=100
@CCFFFFFHHHHAHJJJIJJJJJJIJJIGGGIFIJIIHIIGGJJJJJJJFHIGIJHHHHHHFC#####################################
```
![p5](../RNA-seq/pictures/P5.webp)
```
1、sra数据
sra数据是SRA数据库用于储存二代测序数据的原始数据的一种压缩格式，这种数据格式不能直接进行处理，需要转换成fastq才能进行质控以及去adapt等处理——相当于解压缩

2、fastq文件（简称fq文件）
高通量测序得到的原始图像数据文件，经过碱基识别（base calling）分析转化为原始测序序列（sequenced reads），称之为raw data或raw reads，结果以fastq（简称fq）文件格式存储

链接：https://www.jianshu.com/p/bdfa8f7e5a61    

3. 为何转格式、将fq文件压缩？
因为sra是二进制文件，在Linux下如果用less去查看，它会显示这是个二进制文件，你是否确定打开它。一般我们分析测序数据，是用fastq文件打开分析，所以就需要转格式。没压缩的fq文件通常十几个G，文件一多硬盘就爆炸，所以希望能够以压缩好的gz文件存储，通常只有原始文件的1/8左右，只有原始SRA文件的2倍左右。如果利用gzip命令，处理是单线程，压缩起来很慢，因此需要parallel多线程提高速度

```

  

⑧ 一些尝试记录  
```
# 如果直接在随便一个文件夹下转换格式，不会成功
xuruizhi@DESKTOP-HI65AUV:~$ fastq-dump --split-3 SRR2190795.sra
2022-08-24T11:57:06 fastq-dump.3.0.0 err: item not found while constructing within virtual database module - the path 'SRR2190795.sra' cannot be opened as database or table
fastq-dump quit with error code 3

# 在存储sra文件的文件夹下去转换，ok
xuruizhi@DESKTOP-HI65AUV:~$ cd ~/data/sra

xuruizhi@DESKTOP-HI65AUV:~/data/sra$ fastq-dump --split-3 SRR2190795.sra
Read 15107730 spots for SRR2190795.sra
Written 15107730 spots for SRR2190795.sra

xuruizhi@DESKTOP-HI65AUV:~/data/sra$ ls
SRR2190795.fastq  SRR2240182.sra  SRR2240184.sra  SRR2240186.sra  SRR2240228.sra
SRR2190795.sra    SRR2240183.sra  SRR2240185.sra  SRR2240187.sra
```

## genome data

1. [Ensemble网址](https://asia.ensembl.org/)  
在左侧`All genomes`中，选择物种`Rat`; 在左侧`Download DNA sequence (FASTA)` 下载基因组序列数据; 在右侧的`Download GTF or GFF3 (files for genes, cDNAs, ncRNA, proteins)`下载基因注释文件     

![P1](../RNA-seq/pictures/P1.png) 


2. ensemble中[基因组数据集命名方式](http://ftp.ensembl.org/pub/release-107/fasta/rattus_norvegicus/dna/README)  

* 这些文件始终按照以下模式命名：
```
   <species>.<assembly>.<sequence type>.<id type>.<id>.fa.gz
例：
    Rattus_norvegicus.mRatBN7.2.dna.nonchromosomal.fa.gz  
    Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.1.fa.gz
    Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.fa.gz


<species>：物种的系统名称。

<assembly>：程序集构建名称。

<sequence type>：
  * 'dna' - unmasked未屏蔽的基因组 DNA 序列。
  * 'dna_rm' - masked掩蔽的基因组 DNA。使用 RepeatMasker 工具检测散布的 重复和低复杂性区域，并通过用“N”替换重复来掩盖。
  * 'dna_sm' -  soft-masked软掩蔽基因组 DNA。所有重复和低复杂性区域都已替换为 其核酸碱基的小写版本

<id 类型> 以下之一：
  * 'chromosome' - Ensembl 中大多数物种的顶级坐标系top-level coordinate system
  * 'nonchromosomal' - 包含尚未分配到染色体的 DNA
  * 'seqlevel' - 这通常是sequence scaffolds、块chunks或克隆clones。
     -- 'scaffold' - 短的测序reads（通常来自whole genome shotgun, WGS）组装成较大的序列contigs，但尚未组装成染色体。需要更多的基因组测序来缩小gaps并建立 a tiling path。
     -- 'chunk' - 虽然 contig 序列可以组装成更大区块，有时必须人为地将它们分解为更小的块，称为'chunks'。这是由于注释中的限制pipeline 和 MySQL 有限的记录大小，用于存储序列和注释信息。
     -- 'clone' - 通常这是最小的序列单位。它通常与一个 BAC 克隆的序列或一个 BAC 克隆的序列区域相同，后者形成了tiling path.
<id>：实际的序列标识符。根据 <id type> <id>
          可以代表A chromosome, a scaffold, a contig, a clone的名称..
          seqlevel 文件的字段为空

fa：这些目录中的所有文件都代表FASTA数据库文件

gz：所有文件都使用 GNU Zip 压缩以提高存储效率。
```
* TOPLEVEL    

这些文件包含了 在 Ensembl 模式中标记为toplevel的所有序列区域。 这包括染色体chromsomes、未组装成染色体not assembled into chromosomes的区域和 N 填充的单倍型haplotype/补丁patch区域。

* PRIMARY ASSEMBLY  

初级组装包含所有toplevel序列区域，不包括单倍型和补丁。  
该文件最适合用于`序列相似性搜索`，因为其中补丁和单倍型序列会混淆分析。 如果primary assembly文件不存在，则表明没有单倍型/补丁区域，此时与“toplevel”文件相同。

* special 注意  
一些染色体是单倍体，例如人类的X和Y染色体  
为了比对时能正确输出报告，这些单倍体的assembly和patch区域都会补上同等数量的N   
例如： A patch region with a start position of 1,000,001 will have 1e6 N's added，因这样对齐程序将报告相对于
整个染色体。

人类已对 Y 染色体进行了测序，并对 Y 上的伪常染色体区域pseudoautosomal region (PAR) 进行了注释。 根据定义，PAR 区域在 X 和 Y 染色体上是相同的。 Y染色体文件包含Y染色体减去这些重复的 PAR 区域，即 Y 的唯一部分。  


3. 基因组下载代码：  
目前大鼠的基因组测序版本到了7.2  
可以直接在网页下载，也可用代码  
```
# 下载
cd /mnt/d/project/rat/genome
wget http://ftp.ensembl.org/pub/release-107/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
gzip -d Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz

# 改名（方便后面使用，名字太长一来不方便输入，二来可能会输错）
mv Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa mRatBN7.2.fa
```
4. 对下载的基因组数据进行整理  

！ 和师兄的演示不一样的原因：因为师兄下载的测序版本是6，[旧版本](http://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/dna/)    
`下载的数据id类型是chromosome，新版本是PRIMARY ASSEMBLY`  

![tu](https://github.com/eternal-bug/RNA/blob/master/pic/Rat_genome.png)  


* 下载得到的基因组文件可以查看一下包含哪些染色体，确认文件是否下载正确。
```
cat mRatBN7.2.fa | grep "^>"
```
* 结果：除了1-20号+X+Y+MT之外还有很多别的ID名。这些都是scaffold
```
>19 dna:primary_assembly primary_assembly:mRatBN7.2:19:1:57337602:1 REF
>20 dna:primary_assembly primary_assembly:mRatBN7.2:20:1:54435887:1 REF
>X dna:primary_assembly primary_assembly:mRatBN7.2:X:1:152453651:1 REF
>Y dna:primary_assembly primary_assembly:mRatBN7.2:Y:1:18315841:1 REF
>MT dna:primary_assembly primary_assembly:mRatBN7.2:MT:1:16313:1 REF
>MU150191.1 dna:primary_assembly primary_assembly:mRatBN7.2:MU150191.1:1:1794995:1 REF

>JACYVU010000493.1 dna:primary_assembly primary_assembly:mRatBN7.2:JACYVU010000493.1:1:444596:1 REF
>MU150193.1 dna:primary_assembly primary_assembly:mRatBN7.2:MU150193.1:1:383091:1 REF
```
* 每一条primary_assembly的名称后面还跟了一些描述信息，这些描述信息就是当前组装版本，长度等等信息，但是这个信息会妨碍后面写脚本统计或者一些分析，所以这里最好去掉  
```
# 首先将之前的名称更改一下
mv mRatBN7.2.fa mRatBN7.2.raw.fa

# 然后去除染色体编号后的描述信息
$ cat mRatBN7.2.raw.fa | perl -n -e 'if(m/^>(.+?)(?:\s|$)/){ print ">$1\n";}else{print}' > mRatBN7.2.fa
#单行匹配，如果匹配到了 开头>多个字母空格或者$，将  >多个字母  打印出来

# 删除
$ rm mRatBN7.2.raw.fa
```  
结果：
```
>19
>20
>X
>Y
>MT
>MU150191.1
>MU150189.1
>MU150194.1
>MU150190.1
>MU150195.1
>JACYVU010000493.1
>MU150193.1
>MU150196.1
>MU150197.1
>JACYVU010000705.1
>JACYVU010000706.1
>MU150192.1
>JACYVU010000707.1
```

* 可以使用脚本统计每一条染色体的长度  
```
cat mRatBN7.2.fa | perl -n -e '
    s/\r?\n//;     #\r∶perl语言的转义，回车；删除回车
    if(m/^>(.+?)\s*$/){
        $title = $1;
        push @t, $title; 
    }elsif(defined $title){   #计算第二行数目
        $title_len{$title} += length($_);
    }
    END{
        for my $title (@t){
            print "$title","\t","$title_len{$title}","\n";
        }
    }
'
```
结果：  
```
。。。
19      57337602
20      54435887
X       152453651
Y       18315841
MT      16313
MU150191.1      1794995
MU150189.1      1402623
MU150194.1      648519
MU150190.1      573231
MU150195.1      529129
JACYVU010000493.1       444596
。。。
```  
* 以染色体1 举例  
```
cat mRatBN7.2.fa | perl -n -e '
  if(m/^>/){
    if(m/>1$/){
      $title = 1;
    }else{
      $title = 0;
    }
  }else{
    push @s, $_ if $title;  #title=1，即一号染色体
  }
  END{
    printf ">1\n%s", join("", @s);
  }
' > mRatBN7.2.chr1.fa
```
5. 下载基因组索引文件 - [可选]

方法1. 在[hisat2 官网](https://daehwankimlab.github.io/hisat2/download/#r-norvegicus)上可以找到现成的已经建立好索引的大鼠基因组文件, 点击`	https://genome-idx.s3.amazonaws.com/hisat/rn6_genome.tar.gz`, 下载到了`D:\database`文件夹内   
或通过代码下载  
```
cd /mnt/d/database
wget https://genome-idx.s3.amazonaws.com/hisat/rn6_genome.tar.gz
gzip -d rn6.tar
```

方法2. 自己用命令基于之前下载的基因组文件自行建立   


6.  下载注释信息  
```
# 下载 gff3 格式
cd /mnt/d/project/rat/annotation
wget http://ftp.ensembl.org/pub/release-107/gff3/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.107.gff3.gz
gzip -d Rattus_norvegicus.mRatBN7.2.107.gff3.gz
# 同样的也改名
mv Rattus_norvegicus.mRatBN7.2.107.gff3 mRatBN7.2.gff
# 使用head查看部分
head mRatBN7.2.gff

# gff查看结果：
xuruizhi@DESKTOP-HI65AUV:/mnt/d/project/rat/annotation$ head mRatBN7.2.gff
##gff-version 3
##sequence-region   1 1 260522016
##sequence-region   10 1 107211142
##sequence-region   11 1 86241447
##sequence-region   12 1 46669029
##sequence-region   13 1 106807694
##sequence-region   14 1 104886043
##sequence-region   15 1 101769107
##sequence-region   16 1 84729064
##sequence-region   17 1 86533673 


# 下载 gtf 格式  
cd /mnt/d/project/rat/annotation
wget http://ftp.ensembl.org/pub/release-107/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.107.gtf.gz
gzip -d 	Rattus_norvegicus.mRatBN7.2.107.gtf.gz

# gtf查看结果

xuruizhi@DESKTOP-HI65AUV:/mnt/d/project/rat/annotation$ head  mRatBN7.2.107.gtf

#gtf文件开头描述了这个注释数据的基本信息，比如版本号，更新时间，组装的NCBI的Assembly编号等等，后面每一行表示描述信息，说明了在哪条染色体的什么位置是什么东西。
#!genome-build mRatBN7.2
#!genome-version mRatBN7.2
#!genome-date 2020-11
#!genome-build-accession GCA_015227675.2
#!genebuild-last-updated 2021-02
1       ensembl gene    36112690        36122387        .       -       .       gene_id "ENSRNOG00000066169"; gene_version "1"; gene_source "ensembl"; gene_biotype "protein_coding";
# 比如该行表示在1号染色体负链上 36112690-36122387 这个范围内有一个基因编号为ENSRNOG00000066169的基因

1       ensembl transcript      36112690        36122387        .       -       .       gene_id "ENSRNOG00000066169"; gene_version "1"; transcript_id "ENSRNOT00000101581"; transcript_version "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding";
1       ensembl exon    36122324        36122387        .       -       .       gene_id "ENSRNOG00000066169"; gene_version "1"; transcript_id "ENSRNOT00000101581"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSRNOE00000618632"; exon_version "1";
1       ensembl CDS     36122324        36122387        .       -       0       gene_id "ENSRNOG00000066169"; gene_version "1"; transcript_id "ENSRNOT00000101581"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "ENSRNOP00000083062"; protein_version "1";
1       ensembl exon    36121478        36121512        .       -       .       gene_id "ENSRNOG00000066169"; gene_version "1"; transcript_id "ENSRNOT00000101581"; transcript_version "1"; exon_number "2"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSRNOE00000610554"; exon_version "1";
```

7. ensemble中  
[gff 注释信息命名方式](http://ftp.ensembl.org/pub/release-107/gff3/rattus_norvegicus/README)  
[gtf 注释信息命名方式](http://ftp.ensembl.org/pub/release-107/gtf/rattus_norvegicus/README)  

两者注释区别：gff是先判断该序列属于gene还是转录本等，呈现出不同的主要特征信息；但是gtf是所有的信息全部呈现出来

gff 注释信息：  

① gene features的类型：   

*  "gene" 代表 protein-coding genes 编码蛋白的基因  
*  "ncRNA_gene" 代表 RNA genes   RNA基因   
*  "pseudogene" 代表 pseudogenes假基因  

② transcript features的类型:
   * "mRNA" 代表 protein-coding transcripts 编码蛋白的转录本
   * a specific type or RNA transcript such as "snoRNA" or "lnc_RNA"  RNA转录本
   * "pseudogenic_transcript" for pseudogenes 假基因转录本
所有的转录本都和exon特征相关。编码蛋白的转录本和 "CDS", "five_prime_UTR", and "three_prime_UTR" 特征相关。


③ Attributes for feature types:
* region types: "five_prime_UTR", and "three_prime_UTR"
    * ID: 唯一识别号, 形式 "<region_type>:<region_name>"
    * [Alias]：别名，以逗号分隔的别名列表，通常包括 INSDC accession
    * [Is_circular]: 标志环形区域 circular regions
 * gene types:
    * ID:  唯一识别号, 形式 "gene:<gene_stable_id>"
    * biotype: Ensembl 生物型, e.g. "蛋白质编码", "假基因"
    * gene_id: Ensembl 基因稳定 ID  
    * version: Ensembl 基因版本   
    * [Name]： 基因名称
    * [description]： 基因描述
 * transcript types:  
    * ID: 唯一识别号, 形式 "transcript:<transcript_stable_id>"
    * Parent: 基因标识符, 形式 "gene:<gene_stable_id>"
    * biotype: Ensembl 生物型, e.g. "蛋白质编码", "假基因" 
    * transcript_id: Ensembl 转录本 stable ID
    * version: Ensembl 转录本版本
    * [Note]: 如果转录序列已被编辑 (i.e. 和基因组序列不同), 编辑在注释中描述。
 * exon
    * Parent: Transcript identifier, 形式 "transcript:<transcript_stable_id>"
    * exon_id: Ensembl 外显子 stable ID
    * version: Ensembl 外显子版本
    * constitutive: 组成型外显子，标志着该外显子在所有转录本中均存在
    * rank:  代表5'->3' ordering of exons的整数
 * CDS  
 CDS（Coding sequence）是指成熟mRNA中可以被翻译为蛋白质的编码序列区域，自起始密码子开始至终止密码子结束。
    * ID:  唯一识别号, 形式  "CDS:<protein_stable_id>"
    * Parent: Transcript identifier，形式"transcript:<transcript_stable_id>"
    * protein_id: Ensembl 蛋白 stable ID
    * version: Ensembl 蛋白版本

④ 元数据Metadata：
 * 基因组构建 - 构建assembly的标识符，例如GRCh37.p11
 * 基因组版本 - 此assembly的版本，例如GRCh37
 * 基因组日期 - 此assembly的发布日期，例如2009-02
 * 基因组构建加入 - 基因组加入，例如GCA_000001405.14
 * genebuild-last-updated - 最后一次genebuild更新的日期，例如2013-09


⑤ FILE NAMES： 
  <species>.<assembly>.<_version>.gff3.gz  

对于预测的基因集，在名称文件中添加了一个额外的 abinitio 标志。  
  <species>.<assembly>.<version>.abinitio.gff3.gz  

```
e.g. 
GL476399        Pmarinus_7.0    supercontig     1       4695893 .       .       .       ID=supercontig:GL476399;Alias=scaffold_71
GL476399        ensembl gene                  2596494   2601138 .       +       .       ID=gene:ENSPMAG00000009070;Name=TRYPA3;biotype=protein_coding;description=Trypsinogen A1%3B Trypsinogen a3%3B Uncharacterized protein  [Source:UniProtKB/TrEMBL%3BAcc:O42608];logic_name=ensembl;version=1
```


[gtf 注释信息](http://ftp.ensembl.org/pub/release-107/gtf/rattus_norvegicus/README)   


 GTF (General Transfer Format)   
① FILE NAMES：   
  <species>.<assembly>.<_version>.gtf.gz    

对于预测的基因集，在名称文件中添加了一个额外的 abinitio 标志。    
  <species>.<assembly>.<version>.abinitio.gtf.gz  

② Fields  

Fields are tab-separated. Also, all but the final field in each 
feature line must contain a value; "empty" columns are denoted 
with a '.'  
 
    seqname   - name of the chromosome or scaffold; chromosome names 
                without a 'chr' 
    source    - name of the program that generated this feature, or 
                the data source (database or project name)
    feature   - feature type name. Current allowed features are
                {gene, transcript, exon, CDS, Selenocysteine, start_codon,
                stop_codon and UTR}
    start     - start position of the feature, with sequence numbering 
                starting at 1.
    end       - end position of the feature, with sequence numbering 
                starting at 1.
    score     - a floating point value indiciating the score of a feature
    strand    - defined as + (forward) or - (reverse).
    frame     - one of '0', '1' or '2'. Frame indicates the number of base pairs
                before you encounter a full codon. '0' indicates the feature 
                begins with a whole codon. '1' indicates there is an extra
                base (the 3rd base of the prior codon) at the start of this feature.
                '2' indicates there are two extra bases (2nd and 3rd base of the 
                prior exon) before the first codon. All values are given with
                relation to the 5' end.
    attribute - a semicolon-separated list of tag-value pairs (separated by a space), 
                providing additional information about each feature. A key can be
                repeated multiple times.

③ Attributes  

The following attributes are available. All attributes are semi-colon
separated pairs of keys and values.  

- gene_id: The stable identifier for the gene
- gene_version: The stable identifier version for the gene
- gene_name: The official symbol of this gene
- gene_source: The annotation source for this gene
- gene_biotype: The biotype of this gene
- transcript_id: The stable identifier for this transcript
- transcript_version: The stable identifier version for this transcript
- transcript_name: The symbold for this transcript derived from the gene name
- transcript_source: The annotation source for this transcript
- transcript_biotype: The biotype for this transcript
- exon_id: The stable identifier for this exon
- exon_version: The stable identifier version for this exon
- exon_number: Position of this exon in the transcript
- ccds_id: CCDS identifier linked to this transcript
- protein_id: Stable identifier for this transcript's protein
- protein_version: Stable identifier version for this transcript's protein
- tag: A collection of additional key value tags
- transcript_support_level: Ranking to assess how well a transcript is supported (from 1 to 5)

④ Tags  

Tags are additional flags used to indicate attibutes of the transcript.  

- CCDS: Flags this transcript as one linked to a CCDS record
- seleno: Flags this transcript has a Selenocysteine edit. Look for the Selenocysteine
feature for the position of this on the genome
- cds_end_NF: the coding region end could not be confirmed
- cds_start_NF: the coding region start could not be confirmed
- mRNA_end_NF: the mRNA end could not be confirmed
- mRNA_start_NF: the mRNA start could not be confirmed.
- basic: the transcript is part of the gencode basic geneset





# Pre-alinment

## quality control checking

1. 目的：whether the sequencing quality is qualified or not
2. 使用软件：`FastQC`  
FastQC可用于可视化测序数据中的`碱基质量评分`、`GC含量`、序列长度分布、序列重复水平、k-mer的过度表达，及`引物、接头的污染`。


## pre-alinment QC
1. 目的：adapters and low quality reads trimming
2. 使用软件：`Trim Galore`  
Trim Galore可以自动检测接头序列，质控和去除接头两个步骤一起,适用于多种组学去接头  

```
mkdir -p ../trim/

trim_galore -o /mnt/d/ATAC/output/trim/ --fastqc /mnt/d/ATAC/sequence/*.fastq.gz

# 整合质控结果
cd /mnt/d/ATAC/output/trim/
multiqc .
```

 



# alignment 
1. 目的：将质控后的reads比对到目的基因组上
2. 使用软件： BWA-MEM or Bowtie2，本流程采用`BWA`

 

通常情况下，比对率大于80%视为比对成功。
对于哺乳动物物种，开放染色质检测和差异分析的建议最小mapped reads数为5000万，基于经验和计算估计的TF足迹为2亿。



# Post-alignment processing 

## mapping result sort as BAM 

## remove duplicate reads
1. 目的：将质控后的reads比对到目的基因组上
2. 使用软件： Picard and SAMtools，本流程采用`Picard`


unique mapping reads/rates唯一比对的reads或比例、duplicated read percentages 重复的reads百分比 fragment size distribution 和片段大小分布




去除没有匹配到的、匹配得分较低的、重复的reads；去除线粒体中染色质可及区域及ENCODE blacklisted regions。

1. ATAC-Seq与其他方法不同的一点是需要过滤去除线粒体（如果是植物，还需要过滤叶绿体），因为线粒体DNA是裸露的，也可以被Tn5酶识别切割。  
2. ENCODE blacklisted区域：基因中的重复序列，微卫星序列等，该片段GC含量不稳定，会特异性富集，会呈现假阳性   
Inconsistencies in the underlying annotation exist at regions where assembly has been difficult. For instance, repetitive regions may be collapsed or under-represented in the reference sequence relative to the actual underlying genomic sequence. Resulting analysis of these regions can lead to inaccurate interpretation, as there may be significant enrichment of signal because of amplification of noise.
在人基因组手动注释中发现，这种区域多为particularly rRNA, alpha satellites, and other simple repeats，长度covering on average 45 kb with the largest being 1.4 Mb。[参考文献The ENCODE Blacklist: Identification of Problematic Regions of the Genome](https://mp.weixin.qq.com/s/SS640LNI5QcvChmZNGEOmw)  

3. PCR过程中由于偏好性扩增出现的重复reads
## calculate insert size





# ATAC-seq质量评估

## ATACseqQC:给出国歌质量评估度量值，包括FRiP
还有其他需要评估的特定于 ATAC-seq 的质量度量。通常，一个成功的 ATAC-seq 实验应该生成一个片段大小分布图，其峰值与无核小体区域 (nucleosome-free regions: NFR) (<100 bp) 和单、二、三核小体 (~ 200、400、600 bp) (Fig. 1b) 相对应，呈递减和周期性。来自 NFR 的片段预计会在基因的转录起始位点 (transcription start site, TSS) 附近富集，而来自核小体结合区域的片段预计会在 TSS 附近被耗尽，在 TSS 附近的侧翼区域会有少量富集 (Fig. 1c)。  

![b](../ATAC/pictures/1b.png)    
b: 片段大小在 100bp 和 200bp 左右有明显的富集，表示没有核小体结合和单核小体结合的片段。
![c](../ATAC/pictures/1c.png)  
  c：TSS 富集可视化可以看出，没有核小体结合的片段在 TSS 处富集，而但核小体结合的片段在 TSS 上缺失，在 TSS 两侧富集。  



这些可以通过 ATACseqQC工具进行评估。最后，分别对正链和负链的 reads 进行 + 4bp 和 -5bp 的移位（目标DNA最后产生9bp的重复在ATAC-seq后续分析里要处理。这个长度近似于一个完整的DNA螺旋[参考文章](https://www.jianshu.com/p/13779b89e76b)），以解释 Tn5 转座酶修复损伤 DNA 所产生的 9bp 的重复，并实现 TF footprint 和 motif 相关分析的碱基对分辨率。 
## IDR：样本内的重复性检测，合并一致性peaks

## phantompeakqualtools：评估实验中信噪比、富集信号等

# 上面FastQC➔ trimmomatic➔BWA-MEM➔ATACseqQC

https://github.com/schmitzlab/The-prevalence-evolution-and-chromatin-signatures-of-plant-regulatory-elements/blob/master/Alighment_ATAC-seq_reads/Alighment_ATAC-seq_reads.sh
```
java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE \ 
 -threads $thread -phred33 \ 
 $input.R1.fastq.gz $input.R2.fastq.gz \ 
 ${name}.L.trim.fastq ${name}.L.trimU.fastq ${name}.R.trim.fastq ${name}.R.trimU.fastq \ 
 ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/NexteraPE-PE.fa:2:30:10 \ 
 SLIDINGWINDOW:3:20 LEADING:0 TRAILING:0 MINLEN:30

# bowtie mapping
bowtie $INDEX -t -p 4 -v 2 --best --strata -m 1 -X 1000 -S ${name}.sam \
 -1 ${name}.L.trim.fastq -2 ${name}.R.trim.fastq
 
# sort sam to bam 
samtools sort -O 'bam' -o ${name}.sorted.bam -T tmp ${name}.sam

# remove clonal
java -Xmx20g -classpath /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144 -jar \
  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
  INPUT=${name}.sorted.bam OUTPUT=${name}.clean.bam METRICS_FILE=XXX.txt \
  REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

# bam to bed
bedtools bamtobed -i ${name}.clean.bam > ${name}.bed

# Genome_coverage
bedtools genomecov -i ${name}.bed -split -bg -g $chrom_info > ${name}.bg
wigToBigWig ${name}.bg $chrom_info ${name}.bw

# Tn5_coverage
awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $2 + 5; else print $1, $3 - 6, $3 - 5}' ${name}.bed > ${name}.Tn5.bed
```

# shift
ATAC-seq关心的是在哪里切断，断点才是peak的中心，所以使用shift模型，--shift -75或-100.  



# peak calling 

ATAC-seq 数据分析的第二个主要步骤是识别开放区域（也称为 Peak），后续高级分析以此为基础。目前，MACS2 是 ENCODE ATAC-seq 流程的默认 Peak caller 程序。

  

# Peak differential analysis

csaw 是通过将 edgeR 框架扩展到将基因组分 bin 而开发的。滑动窗口方法被认为可以对基因组中的 reads 进行更多的无偏估计，但是需要严格的 FDR 控制才能正确合并相邻窗口。


# peak annotation

获得 Peak 后，Peak 的注释可将染色质的可及性与基因调控联系起来。通常，Peak 由最接近的基因或调控元件进行注释。ChIPseeker和 ChIPpeakAnno被广泛用于为最接近或重叠的基因、外显子、内含子、启动子、5'UTR、3'UTR 和其他基因组特征分配 Peak。ChIPseeker 和 ChIPpeakAnno 还具有丰富的可视化功能，可用于解释注释结果，例如带有注释的基因组特征的饼图。通常，来自 ATAC-seq 的 Peak 将代表不同的顺式调节元件的混合物，包括增强子和启动子 。在获得基因组特征列表之后，还可以使用 GO,KEGG和 Reactome等数据库进行功能富集分析。通常，Peak 注释会产生生物学和功能上有意义的结果，以供进一步研究。

# Motifs 
尽管 Peak 注释提供了功能解释，但它不能直接解释潜在的机制。开放的染色质可以通过影响 转录因子TF而影响转录，而 TF 通过识别并结合到 DNA 上的特定序列来促进转录。该序列称为 motif，结合位置称为 TF 结合位点（TFBS）。  

有两种类型的基于 motif 或基于 TF 的分析方法：基于序列的 motif 频率或活动预测以及针对 TF 占用的足迹。

`MEME suite`，其中包括`FIMO`用于搜索单个 motif，`MAST` 用于汇总来自多个 motif 的搜索结果，`MCAST` 用于推断由多个 motif 形成的调节模块。这些工具基于统计匹配生成推定的 TFBS 列表。由于 MEME suite 和 PWMScan 具有 Web 应用程序界面，因此更易于访问。

MEME-CentriMo 是一个广泛使用的 web 应用程序，它可以生成可视化报告，而 **chromVAR ** 可以作为 scATAC-seq 的替代方案。

到目前为止所提到的所有工具都间接地从 Peak 区域内发现的 motif 来预测假定的 TFBSs。这种 TFBSs 可能包含大量的误报，并且可能是不完整的和混淆的。这是因为并不是所有的 TFs 都有相同的 motif，来自同一家族的 TFs 可以共享非常相似的 motif。此外，预测的富集或活性变化可能具有微不足道的生物学意义，这妨碍了基于序列的 motif 分析结果的解释

# Footprints

解释 TF 调控的另一种方法是使用 footprint 。ATAC-seq 中的 footprint 指的是一个活跃的 TF 与 DNA 结合并阻止 Tn5 在结合位点内裂解的模式。这在开放的染色质区域留下了一个相对的消耗。因此，活性结合 TFs 的足迹可以用来重建特定样本的调控网络。(简单说就是，TF结合在开放染色质上影响了Tn5的结合)。  


足迹分析工具主要分为两类: De novo 和 motif-centric。

## De novo
De novo 方法根据典型足迹模式 (peak-dip-peak) 的特征，预测所有跨越 Peak 的足迹位置。然后这些假定的足迹位点被用来匹配已知的 motifs 或识别新的 motifs。  
对于  de novo 方法，重要的是数学上定义什么是 footprint 并从 Tn5 裂解偏差中去除 footprint 模式。 在多种工具中，目前只有 HINT-ATAC 处理 ATAC-seq 特定的偏差。 

## motif-centric
Motif-centric 以 motif 为中心的方法侧重于先验的 TFBSs（从motif里面找TFBS），与从头开始的方法相比，考虑到了TF特异性的footprint profiles。

大部分工具都是使用 DNase-seq 数据进行训练的，因此应该使用 ATAC-seq 数据进行再训练，以考虑不同数据的固有偏差。一般来说，由于 TF 和 cell 类型特定的足迹模式具有很大的可变性，因此对它们进行建模仍然很困难。如果对整体 TF 足迹模式在不同条件之间的变化感兴趣，可以使用 BaGFoot[132]。在序列深度归一化和偏差校正后，计算所有 TF  的足迹深度和侧翼可及性。该方法对分析类型 (DNase-seq 或 ATAC-seq)、Peak caller 和偏差校正方法都不错。


de novo 方法对于低质量和 novel motifs 仍然具有优势。尽管由于所选择的分析工具、参数设置和评价指标，不同研究对足迹方法的评价并不一致，但作者认为，由于 HINT-ATAC 具有特定于 ATAC-seq 的偏差校正，因此它可能是一个不错的选择。



# 核小体定位  

 HMMRATAC 和 NucleoATAC都行，是专门针对ATAC-seq 的核小体检测工具。










标准的生物信息学分析包括：

序列分析：将测序读段比对到基因组，并删除重复的reads。
峰发现： 使用MACS 2.1.0，将双端测序中的两个reads用于peak calling。
片段密度的确定： 为了鉴定基因组的转座事件的密度，将基因组分为32 bp的条带，并确定每个条带中的片段数。为此，将reads扩展到200 bp以使数据平滑。
通过随机采样对所有样本的比对reads数进行归一化，以使每一个样本包含所有样本中比对reads数最少的样本相同数目的reads。
活动区域分析： 为了比较两个或多个样本之间的峰，将重叠的峰分组为“活动区域”，这是由最上游峰的起始坐标和最下游峰的终止坐标定义的专有度量。




 作者建议研究人员可以建立一个有效的工作流程，结合 FastQC、trimmomatic 和 BWA-MEM 进行预处理，MACS2 进行 Peak calleing。对于高级分析，作者建议使用 csaw 进行差异 Peak 分析，使用 MEME suite 进行 motif 检测和富集，使用 ChIPseeker 进行注释和可视化，使用 HMMRATAC 进行核小体检测，使用 HINT-ATAC 进行足迹分析。如果 RNA-seq 数据可用，可以使用 PECA 方法重建调控网络。


# 整合多组学数据重建调控网络
## 与 ChIP-seq 进行整合
因为开放的染色质是大多数 TFs 结合的前提条件，所以 ATAC-seq Peak 通常与 TF ChIP-seq Peak 重叠，但通常更宽。因此，TF ChIP-seq 和 ATAC-seq 可以在同一个实验系统中相互验证彼此的质量和可靠性。


第1篇：ATAC-seq的背景介绍以及与ChIP-Seq的异同

第2篇：原始数据的质控、比对和过滤

第3篇：用MACS2软件call peaks

第4篇：对ATAC-Seq/ChIP-seq的质量评估（一）——phantompeakqualtools

第5篇：对ATAC-Seq/ChIP-seq的质量评估（二）——ChIPQC

第6篇：重复样本的处理——IDR

第7篇：用Y叔的ChIPseeker做功能注释

第8篇：用网页版工具进行motif分析

第9篇：差异peaks分析——DiffBind

第10篇：ATAC-Seq、ChIP-Seq、RNA-Seq整合分析

额外篇：文献推荐和解读

```bash
肿瘤资源文章的代码
1. Bowtie2 比对，移除比对到chrM和重复序列

-k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -X 2000 –rg-id # remove repeats的参数
--very-sensitive -X 2000 --rg-id # bowtie2参数
排序去除重复
使用Picard 的MarkDuplicates去除重复。

-f 2 -q 10 -b -@ 20 # 排序参数
VALIDATION_STRINGENCY =LENIENT REMOVE_DUPLICATES = true #去重参数

2. call peaks(MACS2)
这里他们选用固定宽度（fixed-width）的peaks,优点有：1）对大量的peaks进行counts和motif分析时可以减小误差；2）对于大量数据集的可以合并峰得到一致性的peaks; 
使用的是macs2 call peaks,参数如下：

--shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01
同时根据hg38 blacklist过滤，并除去染色体两端以外的峰。
一个样本的overlaps他们是通过迭代移除的方法，首先保留最显著的peak,然后任何与最显著peak有直接overlap的peaks都被移除；接着对另一个最显著性的peak进行相同的操作，最终保留所有更显著的peaks，移除与其有直接overlaps的peaks

3. ATAC-seq数据分析—— 构建counts矩阵并标准化
为了获得每个峰中独立的Tn5插入的数量，首先用RRsamtools “scanbam”对BAM文件矫正Tn5偏移量（“+” stranded +4 bp, “-” stranded -5 bp）并存入Genomic Ranges对象。然后用“countOverlaps”对矫正后的插入位点计数，最终得到 562,709 x 796 counts 矩阵。
counts矩阵用edgeR “cpm(matrix , log = TRUE,prior.count = 5)”标准化，然后用R中的preprocessCore’s “normalize.quantiles”做分位数标准化。

4. ATAC-seq data analysis – Transcription factor footprinting
TF足迹的分析：
一是参考了文章doi: 10.1016/j.celrep.2017.05.003：

首先确定peaks内的TF motif的位置，用pan-cancer peak set 结合CIS-BP motifs计算motif的位置，motifmatchr “matchMotifs(positions = “out”)

然后计算flanking accessibility 和 footprint depth

最后确定哪个TF的足迹与基因的表达是显著相关
通过将flanking accessibility or footprint depth与250个随机的TFs的关联分析生成零均值和标准偏差。

5. ATAC-seq data analysis – chromVAR for transcription factor activity
除了足迹分析，他们还用chromVAR包评估TF的活动，首先用chromVAR deviations函数计算GC矫正偏差，然后将矫正偏差与motif相关的TFs关联，最后5000个转录因子基序和非相关转录因子基因的RNA-seq基因表达之间的随机相关性，以计算每个相关性的FDR。具体参考：Week4— chromVAR:预测染色质可及性相关的转录因子

6. ATAC-seq data analysis – chromVAR for GWAS enrichment
首先从GWAS catalog（https://www.ebi.ac.uk/gwas/docs/file-downloads）下载SNPs位点，过滤和16种癌症类型相关的SNPs位点。

加上连锁不平衡（Linkage Disequilibrium ，LD) 信息（ r 2 > 0.8）
LD信息从haploreg 网站下载 http://archive.broadinstitute.org/mammals/haploreg/data/

移走位于exons或UTR区域的SNPs位点，得到最后的SNP列表

将最后的SNP列表与远端 binarization peak 集overlap，得到一个二元匹配矩阵。每列代表不同癌症癌症类型的GWAS SNP，每行代表一个peak，这个peak来自远端 binarization peak 集。

用chromVAR deviations函数计算GC矫正偏差

用PNAMER将“偏差分数”转换为p值，并使用Bejimi-HocHBG程序调整
```