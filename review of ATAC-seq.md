# ATAC-seq综述

## 研究背景

  绘制细胞状态的变化图谱是理解生物系统的一个关键方面。无论是在发育、分化还是疾病中，细胞状态都受变化的基因表达情况的调控，而基因表达的变化又由调控元件所协调。    

  在真核生物中，遗传物质DNA与组蛋白结合形成核小体，然后折叠和浓缩形成染色质（Zhou BR et al., 2019）。在DNA复制和转录过程中，染色质的某些区域被打开（Li B et al., 2007），即开放染色质（Open Chromatin）；各种转录因子等相关的调控蛋白便有机会结合到开放区域暴露的DNA位点上，调节DNA复制或转录（Maston GA et al., 2006），这种特性叫做染色质的可及性（chromatin accessibility）。该过程中染色质的潜在表观遗传状态也会发生改变，影响基因表达，包括DNA甲基化、组蛋白修饰等。  


  目前对染色质核蛋白结构中表观遗传信息的了解主要来源于高通量、全基因组的方法，用于分析染色质可及性、核小体定位和转录因子占用情况（Buenrostro JD et al.,2013）。最经典的方法是染色质免疫沉淀测序技术（chromatin immunoprecipitation-sequencing,  ChIP -seq），这种靶向分析方法可以为驱动特定细胞状态的表观遗传变化提供关键信息，但需要事先对调控机制有一定了解和预测。在某些情况下，获得更广泛的基因调控图景会更有帮助，特别是当观察到一种现象但表观遗传变化的具体性质仍然未知时。因此，人们开发了其他策略来探索基因调控景观。  


  许多这样的未知分析技术已经发展出来，包括脱氧核糖核酸酶I超敏测序技术（Deoxyribonuclease I hypersensitivity sequencing, DNase-seq）（Crawford GE et al.,2006）、微球菌核酸酶消化测序技术（micrococcal nuclease digestion with sequencing, MNase-seq）（Cui K et al.,2012）、调节元件的甲醛辅助分离技术（Formaldehyde-Assisted Isolation of Regulatory Elements, FAIRE-seq）（Giresi PG et al., 2009）等。但这些方法需要数百万个细胞作为起始材料，涉及复杂、耗时的样品制备，重复性较差，并且不能同时满足探测核小体定位、染色质可及性和TF结合之间的相互作用的需求。  


   为了解决这些限制，2013年美国Stanford大学的William Greenleaf教授研发了一种全新的方法来研究染色质的可进入性，即ATAC-seq（Buenrostro JD et al.,2013）。ATAC-seq全称Assay for Transposase Accessible Chromatin with high-throughput sequencing，即利用转座酶研究染色质可进入性的高通量测序技术。ATAC-seq生成的全基因组调控图谱与DNase-seq和MNase-seq衍生的调控图谱高度相似，同时减少了文库制备的复杂性和实际操作时间。ATAC-seq由于其低输入材料要求（< 50,000个单元）和快速处理等优点，有利于从大量样本中生成数据，而被广泛应用。  


## 实验原理及步骤

### ATAC-seq原理

   ATAC-seq利用DNA转座酶技术实现染色质可及性分析。DNA转座酶可以将自身结合的一段序列随机插入到基因组中。在ATAC-seq试验中，细胞或组织样本在核质分离后，将细胞核单独收集在一起，并通过转座酶对核内的染色质进行打断。紧密包裹的染色质DNA不会受到转座酶的打断，而开放区域的染色质DNA会被转座酶随机插入并打断。将这些打断后的DNA收集在一起，进行后续的建库、测序、分析，即可得到开放染色质的信息。ATAC-seq中的peak，往往是启动子、增强子序列，以及一些调控因子结合的位点。  

### 实验流程

1. 标记：ATAC-Seq技术核心便是孵育连接测序接头的Tn5转座酶。Tn5转座酶的活动可以用粘贴、复制机制来描述。首先提取细胞核，与连有接头的Tn5转座酶进行孵育，这时Tn5转座酶便会特异结合到染色质的开放区域，对开放区域进行切割、打断的同时插入测序接头，而螺旋缠绕的部分不会受到转座酶的影响。  

2. 开放区域的DNA两端与测序接头结合，形成可被测序的片段。  

3. 扩增：DNA片段被扩增。基因组中所有可供转座酶使用的片段都被大量扩增，并用于下游NGS。  

4. NGS：通过测序揭示哪些区域是染色质可及区域，有助于绘制核小体和开放染色质区域图谱。开放与否由激活子 promoters、增强子enhancers和其他调控元件regulatory elements以及是否能与转录机器accessible to transcription machinery 接触决定。染色质的压缩与动态的表观遗传密码有关：DNA甲基化、核小体位置、组蛋白结合以及修饰、转录因子、染色质重塑复合物和非编码RNA。  

5. 下游数据分析   

![process](./pictures/%E5%85%B7%E4%BD%93%E6%AD%A5%E9%AA%A4.png)  


## ATAC - seq技术优缺点  

###  ATAC - seq与其它技术的比较  

前文已经提到了几种全基因组开放染色质的识别技术。在DNase-seq中，DNase I具有内切酶活性，用于切割酶敏感部分的DNA；在MNase-seq中，酶直接消化开放区域，任何与蛋白质结合的DNA都被保留下来。前者可以根据核酸酶敏感位点的存在来识别开放染色质，而后者直接识别蛋白质-DNA结合区域。FAIRE-seq使用甲醛将暴露的DNA固定在染色质中，使用超声波破坏染色质，使用酚类氯仿来分离提取DNA片段。FAIRE-seq只使用物理方式，不使用酶，实验相对简单但带来更多的噪声比（Luo L et al.,2022）。  
将各技术的特点总结如下（Sun Y et al.,2019）：  

| **方法** | **MNase\-seq**             | **DNase\-seq**         | **FAIRE\-seq**         | **ATAC\-seq**                  |
|--------|----------------------------|------------------------|------------------------|---------------------------------|
| 细胞状态   | 细胞的任何状态                    | 细胞的任何状态                | 细胞的任何状态                | 新鲜细胞或缓慢冷却的冷冻细胞                  |
| 原理     | MNase消化不受蛋白质或染色质上核小体保护的DNA | DNase I优先切除没有核小体的DNA序列 | 基于甲醛固定和苯酚\-氯仿萃取的裸DNA分离 | Tn5转座酶插入不受蛋白质或核小体保护的DNA序列并将其切除  |
| 靶向区域   | 关注核小体定位                    | 可及性的染色质区域，集中于转录因子结合位点  | 可及性的染色质区域              | 全基因组可及性染色质区域，包括转录因子，组蛋白修饰       |
| 靶向区域 | 关注核小体定位                                                                              | 可及性的染色质区域，集中于转录因子结合位点                                 | 可及性的染色质区域                                           | 全基因组可及性染色质区域，包括转录因子，组蛋白修饰                                                     |
| 具体特点     | 1.大量的细胞作为起始材料；2.酶的数量需要准确；3.整个核小体和不活跃调控区域的定位；4.通过降解活性区域来检测非活性区域；5.标准分析需要150- 200 M的reads。 | 1.大量的细胞作为起始材料；2.样品制备过程复杂；3.酶的数量需要准确；4.标准分析需要20-50M的reads。 | 1.低信噪比使数据分析变得困难；2.结果在很大程度上依赖于甲醛的固定；3.标准分析需要20-50M的reads。 | 1.较少数量的起始材料；2.通过降低测序深度，标准分析需20–50M的reads；3.方便地获取全基因组中可及性的染色质区域； 4.线粒体数据对结果的准确性有影响  |
| 时间       | 2-3d                                                                                     | 2-3d                                                      | 3-4d                                                     | 3-4h                                                                               |  






![comparing](../ATAC/pictures/%E5%88%86%E6%9E%90%E6%96%B9%E6%B3%95%E6%AF%94%E8%BE%83.jpg)  



### TAC - seq的优点  

与传统的开放染色质研究方法相比， ATAC-seq具有样品需求量更少、制备时间更短、可靠性更高的优点。所需的细胞数量更少，只需500-50,000单位细胞；使用转座酶的特殊优点是，它可以切割DNA，然后将其直接连接到测序接头上，简化了实验过程，缩短了实验周期的同时降低了噪声比.  

### ATAC - seq的局限性  

1.  经常需要额外的步骤来消除线粒体DNA的污染，这可以通过实验和分析来完成。  
2.  一些高通量测序技术不可避免地会产生错误或偏差，包括ATAC-seq。研究发现，Tn5转座酶优先靶向核小体DNA（Sato S et al.,2019）的进出位点，这可能导致测序结果偏偏差。这种转座酶偏差可以通过开发计算工具或改进的统计模型来纠正。  

### ATAC - seq的应用

1. 鉴定重要转录因子：  

根据原理可以知道，ATAC所捕获染色质开放区一般是正在转录的那部分DNA序列的上下游，得到这些序列我们就可以对富集到的序列结合motif 分析，识别哪种转录因子参与了基因表达调控，最常见的就是去研究转录因子结合的启动子区域（对于抗体质量不好的转录因子，尤其有效）
2. 生成转录因子结合区域的特征(footprinting)：  


转录因子结合在DNA上后，它占有的空间阻碍了转座酶Tn5酶切在其他无核小体区域，这样就会留下一个一个小区域，称为足迹（footprint），在这些区域中，reads由高覆盖率峰值突然下降。所以ATAC-seq footprints可以帮助我们查看转录因子在全基因组上结合的状态，主要应用于研究细胞重编程机制，染色质重塑因子，表观修饰对疾病的作用域、T细胞耗竭等等。下面这张图就是已知motif的足迹分析，大概会看到有9个碱基作用的motif

## 参考文献  

1. Zhou BR, Bai Y. Chromatin structures condensed by linker histones. Essays Biochem. 2019;63(1):75-87. doi:10.1042/EBC20180056.   

2. Li B, Carey M, Workman JL. The role of chromatin during transcription. Cell. 2007;128(4):707-719. doi:10.1016/j.cell.2007.01.015.   

3. Maston GA, Evans SK, Green MR. Transcriptional regulatory elements in the human genome. Annu Rev Genomics Hum Genet. 2006;7:29-59. doi:10.1146/annurev.genom.7.080505.115623.   

4. Buenrostro JD, Giresi PG, Zaba LC, et al. Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position. Nat Methods. 2013;10(12):1213-1218. doi:10.1038/nmeth.2688.  

5. Crawford GE, Holt IE, Whittle J, et al. Genome-wide mapping of DNase hypersensitive sites using massively parallel signature sequencing (MPSS). Genome Res. 2006;16(1):123-131. doi:10.1101/gr.4074106.  

6. ui K, Zhao K. Genome-wide approaches to determining nucleosome occupancy in metazoans using MNase-Seq. Methods Mol Biol. 2012;833:413-419. doi:10.1007/978-1-61779-477-3_24.  

7. Giresi PG, Lieb JD. Isolation of active regulatory elements from eukaryotic chromatin using FAIRE (Formaldehyde Assisted Isolation of Regulatory Elements). Methods. 2009;48(3):233-239. doi:10.1016/j.ymeth.2009.03.003.  

8. Luo L, Gribskov M, Wang S. Bibliometric review of ATAC-Seq and its application in gene expression. Brief Bioinform. 2022;23(3):bbac061. doi:10.1093/bib/bbac061.  

9. Sun Y, Miao N, Sun T. Detect accessible chromatin using ATAC-sequencing, from principle to applications. Hereditas. 2019;156:29. doi:10.1186/s41065-019-0105-9.  

10. Sato S, Arimura Y, Kujirai T, et al. Biochemical analysis of nucleosome targeting by Tn5 transposase. Open Biol. 2019;9(8):190116. doi:10.1098/rsob.190116.  

