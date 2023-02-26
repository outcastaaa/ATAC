# Bioinformatics 


For RNA-seq, fastq files were trimmed for poly-A and poly-T sequences using Trimmomatic36 then aligned 
to mouse genome assembly mm10 using hisat2. FeatureCounts was used to obtain read count data and differential 
expression analysis was performed with edgeR using default filtering function (FilterbyExpr) to remove non-expressed 
genes and an FDR of 0.05 for differential expression. 


For ATAC-seq, fastq files were trimmed and aligned to mouse genome assembly mm9 using `Tophat`. Peak calling 
for each sample was performed with `MACS2`. To generate a list of consensus peaks, sequences were identified and filtered 
for loci with at least 5 reads in at least 2 of the 7 samples. These read counts were then normalized using `TMM` in the edgeR 
package, which was also used to generate a list of differentially accessible chromatin regions. Heatmaps of the differentially 
accessible regions determined with `edgeR` as well as for the entire dataset were constructed using the clustermap function 
in the seaborn python package, standardizing rows using z-scores.   

For analysis with `GREAT`, the consensus peaks for PCs and RACMs, as well as the top 500 differentially accessible
peaks, were used to generate `histograms for distance from TSS`. Separately, the 500 differentially accessible peaks were 
uploaded for gene ontology (GO) analysis and gene-region assignments included the closest 2 genes within 1 mb of flanking 
sequence on either end of the region.   

Data visualization of ATAC-seq peaks at specific loci was performed using `pyGenomeTracks` 37
. BAM files were 
first converted to bigwig files using `bamCoverage` in the deeptools package, and y-values were based on the average 
counts per sample
. To generate heatmaps for comparison of Chip-seq datasets with our ATAC-seq data, relevant 
files were downloaded from GEO11 or the ENCODE18,
database as bigwig files. If bigwig files could not be 
downloaded, bam files were downloaded and converted to bigwig files using bamCoverage in the deeptools package.   


Datasets previously aligned to mm10 were converted to mm9 using `crossmap` to facilitate comparison with our ATAC-seq dataset.   


Transcription factor binding site enrichment analysis was performed using Hypergeometric Optimization of Motif 
Enrichment `(HOMER)`  using 50,000 randomly selected background sequences matched for distance from TSS and GC 
content. 
