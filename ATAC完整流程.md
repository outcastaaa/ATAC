
# ATAC-seq分析
- [ATAC-seq分析](#ATAC-seq分析)
- [introduction](#introduction)
- [Purpose](#purpose)
- [Data downloading](#data-downloading)
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





# introduction  

ATAC-seq（Assay for Transposase-Accessible Chromatin with high throughput sequencing） 是2013年由斯坦福大学William J. Greenleaf和Howard Y. Chang实验室开发的用于研究染色质可及性（通常也理解为染色质的开放性）的方法，原理是通过转座酶Tn5容易结合在开放染色质的特性，然后对Tn5酶捕获到的DNA序列进行测序。  


ATAC-seq利用DNA转座酶技术实现染色质可及性分析。DNA转座酶可以将自身结合的一段序列随机插入到基因组中。在ATAC-seq试验中，细胞或组织样本在核质分离后，将细胞核单独收集在一起，并通过转座酶Tn5对核内的染色质进行打断。紧密包裹的染色质DNA不会受到转座酶的打断，而开放区域的染色质DNA会被转座酶随机插入并打断。将这些打断后的DNA收集在一起，进行后续的建库、测序、分析，即可得到开放染色质的信息。ATAC-seq中的peak，往往是启动子、增强子序列，以及一些调控因子结合的位点。  

ATAC-seq可用于：  

- 生成表观基因组图谱  

- 得到在不同组织或不同条件下对应可及性区域    

- 得到核小体位置  

- 鉴定重要转录因子  

- 生成转录因子结合区域的特征(footprinting)  

[具体看该文章](https://github.com/outcastaaa/ATAC/blob/main/review%20of%20ATAC-seq.md)  


数据分析具体流程：  
![数据分析](../ATAC/pictures/%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90%E6%B5%81%E7%A8%8B.jpg)  


ATAC-Seq与其他方法不同的一点是需要过滤去除线粒体（如果是植物，还需要过滤叶绿体），因为线粒体DNA是裸露的，也可以被Tn5酶识别切割。

# raw FASTQ cutadapt

# FASTQ mapping with bowtie2

# mapping result sort as BAM

# remove PCR duplication  

# peak calling 

标准的生物信息学分析包括：

序列分析：将测序读段比对到基因组，并删除重复的reads。
峰发现： 使用MACS 2.1.0，将双端测序中的两个reads用于peak calling。
片段密度的确定： 为了鉴定基因组的转座事件的密度，将基因组分为32 bp的条带，并确定每个条带中的片段数。为此，将reads扩展到200 bp以使数据平滑。
通过随机采样对所有样本的比对reads数进行归一化，以使每一个样本包含所有样本中比对reads数最少的样本相同数目的reads。
活动区域分析： 为了比较两个或多个样本之间的峰，将重叠的峰分组为“活动区域”，这是由最上游峰的起始坐标和最下游峰的终止坐标定义的专有度量。


MACS2 进行 Peak calleing
csaw 进行差异 Peak 分析
MEME suite 进行 motif 检测和富集
ChIPseeker 进行注释和可视化
HMMRATAC 进行核小体检测
HINT-ATAC 进行足迹分析

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