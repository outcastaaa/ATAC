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

## 详细参数
```bash
Version:4.2.0.0


Required Arguments:

--INPUT,-I <String>           One or more input SAM or BAM files to analyze. Must be coordinate sorted.  This argument
                              must be specified at least once. Required.

--METRICS_FILE,-M <File>      File to write duplication metrics to  Required.

--OUTPUT,-O <File>            The output file to write marked records to  Required.


Optional Arguments:

--ADD_PG_TAG_TO_READS <Boolean>
                              Add PG tag to each read in a SAM or BAM  Default value: true. Possible values: {true,
                              false}

--arguments_file <File>       read one or more arguments files and add them to the command line  This argument may be
                              specified 0 or more times. Default value: null.

--ASSUME_SORT_ORDER,-ASO <SortOrder>
                              If not null, assume that the input file has this order even if the header says otherwise.
                              Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate,
                              unknown}  Cannot be used in conjunction with argument(s) ASSUME_SORTED (AS)

--ASSUME_SORTED,-AS <Boolean> If true, assume that the input file is coordinate sorted even if the header says
                              otherwise. Deprecated, used ASSUME_SORT_ORDER=coordinate instead.  Default value: false.
                              Possible values: {true, false}  Cannot be used in conjunction with argument(s)
                              ASSUME_SORT_ORDER (ASO)

--BARCODE_TAG <String>        Barcode SAM tag (ex. BC for 10X Genomics)  Default value: null.

--CLEAR_DT <Boolean>          Clear DT tag from input SAM records. Should be set to false if input SAM doesn't have this
                              tag.  Default true  Default value: true. Possible values: {true, false}

--COMMENT,-CO <String>        Comment(s) to include in the output file's header.  This argument may be specified 0 or
                              more times. Default value: null.

--COMPRESSION_LEVEL <Integer> Compression level for all compressed files created (e.g. BAM and VCF).  Default value: 2.

--CREATE_INDEX <Boolean>      Whether to create an index when writing VCF or coordinate sorted BAM output.  Default
                              value: false. Possible values: {true, false}

--CREATE_MD5_FILE <Boolean>   Whether to create an MD5 digest for any BAM or FASTQ files created.    Default value:
                              false. Possible values: {true, false}

--DUPLEX_UMI <Boolean>        Treat UMIs as being duplex stranded.  This option requires that the UMI consist of two
                              equal length strings that are separated by a hyphen (e.g. 'ATC-GTC'). Reads are considered
                              duplicates if, in addition to standard definition, have identical normalized UMIs.  A UMI
                              from the 'bottom' strand is normalized by swapping its content around the hyphen (eg.
                              ATC-GTC becomes GTC-ATC).  A UMI from the 'top' strand is already normalized as it is.
                              Both reads from a read pair considered top strand if the read 1 unclipped 5' coordinate is
                              less than the read 2 unclipped 5' coordinate. All chimeric reads and read fragments are
                              treated as having come from the top strand. With this option is it required that the
                              BARCODE_TAG hold non-normalized UMIs. Default false.  Default value: false. Possible
                              values: {true, false}

--DUPLICATE_SCORING_STRATEGY,-DS <ScoringStrategy>
                              The scoring strategy for choosing the non-duplicate among candidates.  Default value:
                              SUM_OF_BASE_QUALITIES. Possible values: {SUM_OF_BASE_QUALITIES,
                              TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM}

--GA4GH_CLIENT_SECRETS <String>
                              Google Genomics API client_secrets.json file path.  Default value: client_secrets.json.

--help,-h <Boolean>           display the help message  Default value: false. Possible values: {true, false}

--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP,-MAX_FILE_HANDLES <Integer>
                              Maximum number of file handles to keep open when spilling read ends to disk. Set this
                              number a little lower than the per-process maximum number of file that may be open. This
                              number can be found by executing the 'ulimit -n' command on a Unix system.  Default value:
                              8000.

--MAX_OPTICAL_DUPLICATE_SET_SIZE <Long>
                              This number is the maximum size of a set of duplicate reads for which we will attempt to
                              determine which are optical duplicates.  Please be aware that if you raise this value too
                              high and do encounter a very large set of duplicate reads, it will severely affect the
                              runtime of this tool.  To completely disable this check, set the value to -1.  Default
                              value: 300000.

--MAX_RECORDS_IN_RAM <Integer>When writing files that need to be sorted, this will specify the number of records stored
                              in RAM before spilling to disk. Increasing this number reduces the number of file handles
                              needed to sort the file, and increases the amount of RAM needed.  Default value: 500000.

--MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP,-MAX_SEQS <Integer>
                              This option is obsolete. ReadEnds will always be spilled to disk.  Default value: 50000.

--MOLECULAR_IDENTIFIER_TAG <String>
                              SAM tag to uniquely identify the molecule from which a read was derived.  Use of this
                              option requires that the BARCODE_TAG option be set to a non null value.  Default null.
                              Default value: null.

--OPTICAL_DUPLICATE_PIXEL_DISTANCE <Integer>
                              The maximum offset between two duplicate clusters in order to consider them optical
                              duplicates. The default is appropriate for unpatterned versions of the Illumina platform.
                              For the patterned flowcell models, 2500 is moreappropriate. For other platforms and
                              models, users should experiment to find what works best.  Default value: 100.

--PROGRAM_GROUP_COMMAND_LINE,-PG_COMMAND <String>
                              Value of CL tag of PG record to be created. If not supplied the command line will be
                              detected automatically.  Default value: null.

--PROGRAM_GROUP_NAME,-PG_NAME <String>
                              Value of PN tag of PG record to be created.  Default value: MarkDuplicates.

--PROGRAM_GROUP_VERSION,-PG_VERSION <String>
                              Value of VN tag of PG record to be created. If not specified, the version will be detected
                              automatically.  Default value: null.

--PROGRAM_RECORD_ID,-PG <String>
                              The program record ID for the @PG record(s) created by this program. Set to null to
                              disable PG record creation.  This string may have a suffix appended to avoid collision
                              with other program record IDs.  Default value: MarkDuplicates.

--QUIET <Boolean>             Whether to suppress job-summary info on System.err.  Default value: false. Possible
                              values: {true, false}

--READ_NAME_REGEX <String>    MarkDuplicates can use the tile and cluster positions to estimate the rate of optical
                              duplication in addition to the dominant source of duplication, PCR, to provide a more
                              accurate estimation of library size. By default (with no READ_NAME_REGEX specified),
                              MarkDuplicates will attempt to extract coordinates using a split on ':' (see Note below).
                              Set READ_NAME_REGEX to 'null' to disable optical duplicate detection. Note that without
                              optical duplicate counts, library size estimation will be less accurate. If the read name
                              does not follow a standard Illumina colon-separation convention, but does contain tile and
                              x,y coordinates, a regular expression can be specified to extract three variables:
                              tile/region, x coordinate and y coordinate from a read name. The regular expression must
                              contain three capture groups for the three variables, in order. It must match the entire
                              read name.   e.g. if field names were separated by semi-colon (';') this example regex
                              could be specified      (?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*$ Note that if no
                              READ_NAME_REGEX is specified, the read name is split on ':'.   For 5 element names, the
                              3rd, 4th and 5th elements are assumed to be tile, x and y values.   For 7 element names
                              (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values.
                              Default value: <optimized capture of last three ':' separated fields as numeric values>.

--READ_ONE_BARCODE_TAG <String>
                              Read one barcode SAM tag (ex. BX for 10X Genomics)  Default value: null.

--READ_TWO_BARCODE_TAG <String>
                              Read two barcode SAM tag (ex. BX for 10X Genomics)  Default value: null.

--REFERENCE_SEQUENCE,-R <File>Reference sequence file.  Default value: null.

--REMOVE_DUPLICATES <Boolean> If true do not write duplicates to the output file instead of writing them with
                              appropriate flags set.  Default value: false. Possible values: {true, false}

--REMOVE_SEQUENCING_DUPLICATES <Boolean>
                              If true remove 'optical' duplicates and other duplicates that appear to have arisen from
                              the sequencing process instead of the library preparation process, even if
                              REMOVE_DUPLICATES is false. If REMOVE_DUPLICATES is true, all duplicates are removed and
                              this option is ignored.  Default value: false. Possible values: {true, false}

--SORTING_COLLECTION_SIZE_RATIO <Double>
                              This number, plus the maximum RAM available to the JVM, determine the memory footprint
                              used by some of the sorting collections.  If you are running out of memory, try reducing
                              this number.  Default value: 0.25.

--TAG_DUPLICATE_SET_MEMBERS <Boolean>
                              If a read appears in a duplicate set, add two tags. The first tag, DUPLICATE_SET_SIZE_TAG
                              (DS), indicates the size of the duplicate set. The smallest possible DS value is 2 which
                              occurs when two reads map to the same portion of the reference only one of which is marked
                              as duplicate. The second tag, DUPLICATE_SET_INDEX_TAG (DI), represents a unique identifier
                              for the duplicate set to which the record belongs. This identifier is the index-in-file of
                              the representative read that was selected out of the duplicate set.  Default value: false.
                              Possible values: {true, false}

--TAGGING_POLICY <DuplicateTaggingPolicy>
                              Determines how duplicate types are recorded in the DT optional attribute.  Default value:
                              DontTag. Possible values: {DontTag, OpticalOnly, All}

--TMP_DIR <File>              One or more directories with space available to be used by this program for temporary
                              storage of working files  This argument may be specified 0 or more times. Default value:
                              null.

--USE_JDK_DEFLATER,-use_jdk_deflater <Boolean>
                              Use the JDK Deflater instead of the Intel Deflater for writing compressed output  Default
                              value: false. Possible values: {true, false}

--USE_JDK_INFLATER,-use_jdk_inflater <Boolean>
                              Use the JDK Inflater instead of the Intel Inflater for reading compressed input  Default
                              value: false. Possible values: {true, false}

--VALIDATION_STRINGENCY <ValidationStringency>
                              Validation stringency for all SAM files read by this program.  Setting stringency to
                              SILENT can improve performance when processing a BAM file in which variable-length data
                              (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT.
                              Possible values: {STRICT, LENIENT, SILENT}

--VERBOSITY <LogLevel>        Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING,
                              INFO, DEBUG}

--version <Boolean>           display the version number for this tool  Default value: false. Possible values: {true,
                              false}


Advanced Arguments:

--showHidden,-showHidden <Boolean>
                              display hidden arguments  Default value: false. Possible values: {true, false}


Argument INPUT was missing: Argument 'INPUT' is required
```


