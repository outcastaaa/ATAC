
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

标准的生物信息学分析包括：

序列分析：将测序读段比对到基因组，并删除重复的reads。
峰发现： 使用MACS 2.1.0，将双端测序中的两个reads用于peak calling。
片段密度的确定： 为了鉴定基因组的转座事件的密度，将基因组分为32 bp的条带，并确定每个条带中的片段数。为此，将reads扩展到200 bp以使数据平滑。
通过随机采样对所有样本的比对reads数进行归一化，以使每一个样本包含所有样本中比对reads数最少的样本相同数目的reads。
活动区域分析： 为了比较两个或多个样本之间的峰，将重叠的峰分组为“活动区域”，这是由最上游峰的起始坐标和最下游峰的终止坐标定义的专有度量。

# introduction  

ATAC-seq利用DNA转座酶技术实现染色质可及性分析。DNA转座酶可以将自身结合的一段序列随机插入到基因组中。在ATAC-seq试验中，细胞或组织样本在核质分离后，将细胞核单独收集在一起，并通过转座酶Tn5对核内的染色质进行打断。紧密包裹的染色质DNA不会受到转座酶的打断，而开放区域的染色质DNA会被转座酶随机插入并打断。将这些打断后的DNA收集在一起，进行后续的建库、测序、分析，即可得到开放染色质的信息。ATAC-seq中的peak，往往是启动子、增强子序列，以及一些调控因子结合的位点。  

[具体看该文章](https://github.com/outcastaaa/ATAC/blob/main/review%20of%20ATAC-seq.md)  



