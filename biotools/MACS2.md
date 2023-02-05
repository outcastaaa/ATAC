[参考文章:包括原理、参数](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/blob/master/sessionV/lessons/04_peak_calling_macs.md)
## 详细参数
There are seven [major functions](https://github.com/taoliu/MACS#usage-of-macs2) available in MACS2 serving as sub-commands. We will only cover `callpeak` in this lesson, but you can use `macs2 COMMAND -h` to find out more, if you are interested.

`callpeak` is the main function in MACS2 and can be invoked by typing `macs2 callpeak`. If you type this command without parameters, you will see a full description of commandline options. Here is a shorter list of the commonly used ones: 

**Input file options**

* `-t`: The IP data file (this is the only REQUIRED parameter for MACS)
* `-c`: The control or mock data file
* `-f`: format of input file; Default is "AUTO" which will allow MACS to decide the format automatically.
* `-g`: mappable genome size which is defined as the genome size which can be sequenced; some precompiled values provided.

**Output arguments**

* `--outdir`: MACS2 will save all output files into speficied folder for this option
* `-n`: The prefix string for output files
* `-B/--bdg`: store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files

**Shifting model arguments**

* `-s`: size of sequencing tags. Default, MACS will use the first 10 sequences from your input treatment file to determine it
* `--bw`: The bandwidth which is used to scan the genome ONLY for model building. Can be set to the expected sonication fragment size.
* `--mfold`: upper and lower limit for model building

**Peak calling arguments**

* `-q`: q-value (minimum FDR) cutoff
* `-p`: p-value cutoff (instead of q-value cutoff)
* `--nolambda`: do not consider the local bias/lambda at peak candidate regions
* `--broad`: broad peak calling

> **NOTE:** Relaxing the q-value does not behave as expected in this case since it is partially tied to peak widths. Ideally, if you relaxed the thresholds, you would simply get more peaks but with MACS2 relaxing thresholds also results in wider peaks.

Now that we have a feel for the different ways we can tweak our command, let's set up the command for our run on Nanog-rep1:

```
$ macs2 callpeak -t bowtie2/H1hesc_Nanog_Rep1_aln.bam \
	-c bowtie2/H1hesc_Input_Rep1_aln.bam \
 	-f BAM -g 1.3e+8 \
	-n Nanog-rep1 \
	--outdir macs2
```

The tool is quite verbose so you should see lines of text being printed to the terminal, describing each step that is being carried out. If that runs successfully, go ahead and **re-run the same command but this time let's capture that information into a log file using `2>` to re-direct the stadard error to file**:

```
$ macs2 callpeak -t bowtie2/H1hesc_Nanog_Rep1_aln.bam \
	-c bowtie2/H1hesc_Input_Rep1_aln.bam \
 	-f BAM -g 1.3e+8 \
	-n Nanog-rep1 \
	--outdir macs2 2> macs2/Nanog-rep1-macs2.log
```

Ok, now let's do the same peak calling for the rest of our samples:

```bash
macs2 callpeak -t bowtie2/H1hesc_Nanog_Rep2_aln.bam -c bowtie2/H1hesc_Input_Rep2_aln.bam -f BAM -g 1.3e+8 --outdir macs2 -n Nanog-rep2 2> macs2/Nanog-rep2-macs2.log
	 
macs2 callpeak -t bowtie2/H1hesc_Pou5f1_Rep1_aln.bam -c bowtie2/H1hesc_Input_Rep1_aln.bam -f BAM -g 1.3e+8 --outdir macs2 -n Pou5f1-rep1 2> macs2/Pou5f1-rep1-macs2.log
	 
macs2 callpeak -t bowtie2/H1hesc_Pou5f1_Rep2_aln.bam -c bowtie2/H1hesc_Input_Rep2_aln.bam -f BAM -g 1.3e+8 --outdir macs2 -n Pou5f1-rep2 2> macs2/Pou5f1-rep2-macs2.log

```

## MACS2 Output files

### File formats
Before we start exploring the output of MACS2, we'll briefly talk about the new file formats you will encounter.

**narrowPeak:**

A narrowPeak (.narrowPeak) file is used by the ENCODE project to provide called peaks of signal enrichment based on pooled, normalized (interpreted) data. It is a BED 6+4 format, which means the first 6 columns of a standard BED file  with **4 additional fields**:
[bed文件](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/blob/master/sessionV/img/narrowPeak.png)  
```
1.chromosome 
2.start coordinate
3.end coordinate 
4.name 
5.score 
6.strand
==========Standard BED file fields========
7.signalValue-Measurement of overall enrichment for the region 
8.pValue-Statistical significance(-log10)
9.qValue-Statistical significance using false discovery rate(-log10)
10.peak-Point-source called for thispeak; 0-based offset from chromStart
==========narrowPeak specific fields=========
```



**WIG format:**

Wiggle format (WIG) allows the display of continuous-valued data in a track format. Wiggle format is line-oriented. It is composed of declaration lines and data lines, and require a separate wiggle track definition line. There are two options for formatting wiggle data: variableStep and fixedStep. These formats were developed to allow the file to be written as compactly as possible.

**BedGraph format:**

The BedGraph format also allows display of continuous-valued data in track format. This display type is useful for probability scores and transcriptome data. This track type is similar to the wiggle (WIG) format, but unlike the wiggle format, data exported in the bedGraph format are preserved in their original state. For the purposes of visualization, these can be interchangeable.

### MACS2 output files

	$ cd macs2/
	
	$ ls -lh
	
Let's first move the log files to the `log` directory:

	$ mv *.log ../../logs/
	
Now, there should be 6 files output to the results directory for each of the 4 samples, so a total of 24 files:

* `_peaks.narrowPeak`: BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue
* `_peaks.xls`: a tabular file which contains information about called peaks. Additional information includes pileup and fold enrichment
* `_summits.bed`: peak summits locations for every peak. To find the motifs at the binding sites, this file is recommended
* `_model.R`: an R script which you can use to produce a PDF image about the model based on your data and cross-correlation plot
* `_control_lambda.bdg`: bedGraph format for input sample
* `_treat_pileup.bdg`: bedGraph format for treatment sample

Let's first obtain a summary of how many peaks were called in each sample. We can do this by counting the lines in the `.narrowPeak` files:

	$ wc -l *.narrowPeak

We can also generate plots using the R script file that was output by MACS2. There is a `_model.R` script in the directory. Let's load the R module and run the R script in the command line using the `Rscript` command as demonstrated below:


	$ module load gcc/6.2.0 R/3.4.1
	
	$ Rscript Nanog-rep1_model.r
	
> **NOTE:** We need to load the `gcc/6.2.0` before loading R. You can find out which modules need to be loaded first by using module spider R/3.4.1` 
	
Now you should see a pdf file in your current directory by the same name. Create the plots for each of the samples and move them over to your laptop using `Filezilla`. 

Open up the pdf file for Nanog-rep1. The first plot illustrates **the distance between the modes from which the shift size was determined**. 

[picture of peak-modle](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/blob/master/sessionV/img/model-macs.png)    




The second plot is the  **cross-correlation plot**. This is a graphical representation of the Pearson correlation of positive- and negative- strand tag densities, shifting the strands relative to each other by increasing distance. We will talk about this in more detail in the next lesson.

> **NOTE:** [SPP](http://www.nature.com.ezp-prod1.hul.harvard.edu/nbt/journal/v26/n12/full/nbt.1508.html) is another very commonly used tool for *narrow* peak calling. While we will not be going through the steps for this peak caller in this workshop, we do have [a lesson on SPP](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/blob/may2017/sessionV/lessons/peak_calling_spp.md) that we encourage you to browse through if you are interested in learning more.

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

