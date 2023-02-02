# PICARD
功能非常强大  
## 详细参数  
```bash
picard --version
USAGE: PicardCommandLine <program name> [-h]

Available Programs:
--------------------------------------------------------------------------------------
Base Calling:                                    Tools that process sequencing machine data, e.g. Illumina base calls, and detect sequencing level attributes, e.g. adapters
    CheckIlluminaDirectory                       Asserts the validity for specified Illumina basecalling data.
    CollectIlluminaBasecallingMetrics            Collects Illumina Basecalling metrics for a sequencing run.
    CollectIlluminaLaneMetrics                   Collects Illumina lane metrics for the given BaseCalling analysis directory.
    ExtractIlluminaBarcodes                      Tool determines the barcode for each read in an Illumina lane.
    IlluminaBasecallsToFastq                     Generate FASTQ file(s) from Illumina basecall read data.
    IlluminaBasecallsToSam                       Transforms raw Illumina sequencing data into an unmapped SAM, BAM or CRAM file.
    MarkIlluminaAdapters                         Reads a SAM/BAM/CRAM file and rewrites it with new adapter-trimming tags.


--------------------------------------------------------------------------------------
Diagnostics and Quality Control:                 Tools that collect sequencing quality related and comparative metrics
    AccumulateQualityYieldMetrics                Combines multiple QualityYieldMetrics files into a single file.
    AccumulateVariantCallingMetrics              Combines multiple Variant Calling Metrics files into a single file
    BamIndexStats                                Generate index statistics from a BAM file
    CalculateFingerprintMetrics                  Calculate statistics on fingerprints, checking their viability
    CalculateReadGroupChecksum                   Creates a hash code based on the read groups (RG).
    CheckDuplicateMarking                        Checks the consistency of duplicate markings.
    CheckFingerprint                             Computes a fingerprint from the supplied input (SAM/BAM/CRAM or VCF) file and compares it to the provided genotypes
    CheckTerminatorBlock                         Asserts the provided gzip file's (e.g., BAM) last block is well-formed; RC 100 otherwise
    ClusterCrosscheckMetrics                     Clusters the results of a CrosscheckFingerprints run by LOD score
    CollectAlignmentSummaryMetrics               <b>Produces a summary of alignment metrics from a SAM or BAM file.</b>
    CollectArraysVariantCallingMetrics           Collects summary and per-sample from the provided arrays VCF file
    CollectBaseDistributionByCycle               Chart the nucleotide distribution per cycle in a SAM or BAM file
    CollectGcBiasMetrics                         Collect metrics regarding GC bias.
    CollectHiSeqXPfFailMetrics                   Classify PF-Failing reads in a HiSeqX Illumina Basecalling directory into various categories.
    CollectHsMetrics                             Collects hybrid-selection (HS) metrics for a SAM or BAM file.
    CollectIndependentReplicateMetrics           **EXPERIMENTAL - USE AT YOUR OWN RISK** Estimates the rate of independent replication rate of reads within a bam.

    CollectInsertSizeMetrics                     Collect metrics about the insert size distribution of a paired-end library.
    CollectJumpingLibraryMetrics                 Collect jumping library metrics.
    CollectMultipleMetrics                       Collect multiple classes of metrics.
    CollectOxoGMetrics                           Collect metrics to assess oxidative artifacts.
    CollectQualityYieldMetrics                   Collect metrics about reads that pass quality thresholds and Illumina-specific filters.
    CollectRawWgsMetrics                         Collect whole genome sequencing-related metrics.
    CollectRnaSeqMetrics                         Produces RNA alignment metrics for a SAM or BAM file.
    CollectRrbsMetrics                           <b>Collects metrics from reduced representation bisulfite sequencing (Rrbs) data.</b>
    CollectSamErrorMetrics                       Program to collect error metrics on bases stratified in various ways.
    CollectSequencingArtifactMetrics             Collect metrics to quantify single-base sequencing artifacts.
    CollectTargetedPcrMetrics                    Calculate PCR-related metrics from targeted sequencing data.
    CollectVariantCallingMetrics                 Collects per-sample and aggregate (spanning all samples) metrics from the provided VCF file
    CollectWgsMetrics                            Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments.
    CollectWgsMetricsWithNonZeroCoverage         **EXPERIMENTAL - USE AT YOUR OWN RISK** Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments.
    CompareMetrics                               Compare two metrics files.
    CompareSAMs                                  Compare two input SAM/BAM/CRAM files.
    ConvertHaplotypeDatabaseToVcf                Convert Haplotype database file to vcf
    ConvertSequencingArtifactToOxoG              Extract OxoG metrics from generalized artifacts metrics.
    CrosscheckFingerprints                       Checks that all data in the input files appear to have come from the same individual
    CrosscheckReadGroupFingerprints              DEPRECATED: USE CrosscheckFingerprints.
    EstimateLibraryComplexity                    Estimates the numbers of unique molecules in a sequencing library.
    ExtractFingerprint                           Computes a fingerprint from the input file.
    IdentifyContaminant                          Computes a fingerprint from the supplied SAM/BAM file, given a contamination estimate.
    LiftOverHaplotypeMap                         Lifts over a haplotype database from one reference to another
    MeanQualityByCycle                           Collect mean quality by cycle.
    QualityScoreDistribution                     Chart the distribution of quality scores.
    ValidateSamFile                              Validates a SAM/BAM/CRAM file.
    ViewSam                                      Prints a SAM or BAM file to the screen

--------------------------------------------------------------------------------------
Genotyping Arrays Manipulation:                  Tools that manipulate data generated by Genotyping arrays
    BpmToNormalizationManifestCsv                Program to convert an Illumina bpm file into a bpm.csv file.
    CombineGenotypingArrayVcfs                   Program to combine multiple genotyping array VCF files into one VCF.
    CompareGtcFiles                              Compares two GTC files.
    CreateBafRegressMetricsFile                  Program to generate a picard metrics file from the output of the bafRegress tool.
    CreateExtendedIlluminaManifest               Create an Extended Illumina Manifest for usage by the Picard tool GtcToVcf
    CreateVerifyIDIntensityContaminationMetricsFile    Program to generate a picard metrics file from the output of the VerifyIDIntensity tool.
    GtcToVcf                                     Program to convert an Illumina GTC file to a VCF
    MergePedIntoVcf                              Program to merge a single-sample ped file from zCall into a single-sample VCF.
    VcfToAdpc                                    Program to convert an Arrays VCF to an ADPC file.

--------------------------------------------------------------------------------------
Intervals Manipulation:                          Tools that process genomic intervals in various formats
    BedToIntervalList                            Converts a BED file to a Picard Interval List.
    IntervalListToBed                            Converts an Picard IntervalList file to a BED file.
    IntervalListTools                            A tool for performing various IntervalList manipulations
    LiftOverIntervalList                         Lifts over an interval list from one reference build to another.

--------------------------------------------------------------------------------------
Other:                                           Miscellaneous tools, e.g. those that aid in data streaming
    FifoBuffer                                   Provides a large, FIFO buffer that can be used to buffer input and output streams between programs.
    SortGff                                      Sorts a gff3 file, and adds flush directives

--------------------------------------------------------------------------------------
Read Data Manipulation:                          Tools that manipulate read data in SAM, BAM or CRAM format
    AddCommentsToBam                             Adds comments to the header of a BAM file.
    AddOATag                                     Record current alignment information to OA tag.
    AddOrReplaceReadGroups                       Assigns all the reads in a file to a single new read-group.
    BamToBfq                                     Converts a BAM file into a BFQ (binary fastq formatted) file
    BuildBamIndex                                Generates a BAM index ".bai" file.
    CleanSam                                     Cleans a SAM/BAM/CRAM files, soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads
    CollectDuplicateMetrics                      Collect Duplicate metrics from marked file.
    DownsampleSam                                Downsample a SAM or BAM file.
    FastqToSam                                   Converts a FASTQ file to an unaligned BAM or SAM file
    FilterSamReads                               Subsets reads from a SAM/BAM/CRAM file by applying one of several filters.
    FixMateInformation                           Verify mate-pair information between mates and fix if needed.
    GatherBamFiles                               Concatenate efficiently BAM files that resulted from a scattered parallel analysis
    MarkDuplicates                               Identifies duplicate reads.
    MarkDuplicatesWithMateCigar                  Identifies duplicate reads, accounting for mate CIGAR.
    MergeBamAlignment                            Merge alignment data from a SAM or BAM with data in an unmapped BAM file.

    MergeSamFiles                                Merges multiple SAM/BAM/CRAM (and/or) files into a single file.
    PositionBasedDownsampleSam                   Downsample a SAM or BAM file to retain a subset of the reads based on the reads location in each tile in the flowcell.
    ReorderSam                                   Reorders reads in a SAM or BAM file to match ordering in a second reference file.
    ReplaceSamHeader                             Replaces the SAMFileHeader in a SAM/BAM/CRAM file.
    RevertOriginalBaseQualitiesAndAddMateCigar   Reverts the original base qualities and adds the mate cigar tag to read-group files
    RevertSam                                    Reverts SAM/BAM/CRAM files to a previous state.
    SamFormatConverter                           Convert a BAM file to a SAM file, or a SAM to a BAM
    SamToFastq                                   Converts a SAM/BAM/CRAM file to FASTQ.
    SamToFastqWithTags                           Converts a SAM or BAM file to FASTQ alongside FASTQs created from tags.
    SetNmAndUqTags                               DEPRECATED: Use SetNmMdAndUqTags instead.
    SetNmMdAndUqTags                             Fixes the NM, MD, and UQ tags in a SAM/BAM/CRAM file
    SimpleMarkDuplicatesWithMateCigar            **EXPERIMENTAL - USE AT YOUR OWN RISK** Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules.
    SortSam                                      Sorts a SAM, BAM or CRAM file.
    SplitSamByLibrary                            Splits a SAM/BAM/CRAM file into individual files by library
    SplitSamByNumberOfReads                      Splits a SAM/BAM/CRAM file to multiple files.
    UmiAwareMarkDuplicatesWithMateCigar          **EXPERIMENTAL - USE AT YOUR OWN RISK** Identifies duplicate reads using information from read positions and UMIs.

--------------------------------------------------------------------------------------
Reference:                                       Tools that analyze and manipulate FASTA format references
    BaitDesigner                                 Designs oligonucleotide baits for hybrid selection reactions.
    CreateSequenceDictionary                     Creates a sequence dictionary for a reference sequence.
    ExtractSequences                             Subsets intervals from a reference sequence to a new FASTA file.
    NonNFastaSize                                Counts the number of non-N bases in a fasta file.
    NormalizeFasta                               Normalizes lines of sequence in a FASTA file to be of the same length.
    ScatterIntervalsByNs                         Writes an interval list created by splitting a reference at Ns.

--------------------------------------------------------------------------------------
Variant Evaluation and Refinement:               Tools that evaluate and refine variant calls, e.g. with annotations not offered by the engine
    FindMendelianViolations                      Finds mendelian violations of all types within a VCF
    GenotypeConcordance                          Calculates the concordance between genotype data of one sample in each of two VCFs - truth (or reference) vs. calls.

--------------------------------------------------------------------------------------
Variant Filtering:                               Tools that filter variants by annotating the FILTER column
    FilterVcf                                    Hard filters a VCF.

--------------------------------------------------------------------------------------
Variant Manipulation:                            Tools that manipulate variant call format (VCF) data
    FixVcfHeader                                 Replaces or fixes a VCF header.
    GatherVcfs                                   Gathers multiple VCF files from a scatter operation into a single VCF file
    LiftoverVcf                                  Lifts over a VCF file from one reference build to another.
    MakeSitesOnlyVcf                             Creates a VCF that contains all the site-level information for all records in the input VCF but no genotype information.
    MakeVcfSampleNameMap                         Creates a TSV from sample name to VCF/GVCF path, with one line per input.
    MergeVcfs                                    Combines multiple variant files into a single variant file
    RenameSampleInVcf                            Renames a sample within a VCF or BCF.
    SortVcf                                      Sorts one or more VCF files.
    SplitVcfs                                    Splits SNPs and INDELs into separate files.
    UpdateVcfSequenceDictionary                  Takes a VCF and a second file that contains a sequence dictionary and updates the VCF with the new sequence dictionary.
    VcfFormatConverter                           Converts VCF to BCF or BCF to VCF.
    VcfToIntervalList                            Converts a VCF or BCF file to a Picard Interval List

--------------------------------------------------------------------------------------
```

## 比对去重参数使用方法
MarkDuplicates——Identifies duplicate reads.
## 去重原理  

该工具的MarkDuplicates方法也可以识别duplicates。但是与samtools不同的是，该工具仅仅是对duplicates做一个标记，只在需要的时候对reads进行去重。  

它不仅考虑reads的比对位置，还会考虑其中的插入错配等情况（即会利用sam/bam文件中的CIGAR值），甚至reads的tail、lane以及flowcell。Picard主要考虑reads的5'端的比对位置，以及每个reads比对上的方向。  

因此我们可以从一定程度上认为，5' 端的位置、方向、以及碱基比对情况相同，Picard就将这些reads中碱基比对值Q>15的看作是best pair而其他的reads则当作是duplicate reads。甚至当reads的长度不同时，Picard依然利用上述原理进行去重。  

对Picard来说，reads的5' 端信息更为重要。若duplicates是PCR重复，那么它们的序列不一定完全相同。但是由于PCR扩增时，酶的前进方向是5'->3'方向，PCR重复序列中5' 端的部分相似的可能性更高。  




