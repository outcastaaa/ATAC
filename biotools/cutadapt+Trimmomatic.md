# cutadapt
[cutadapt使用参考文章](https://www.cnblogs.com/xudongliang/p/6404958.html#:~:text=%E7%94%A8%E6%B3%95%EF%BC%9A%20cutadapt%20-q%2010%20-o%20output.fastq%20input.fastq%20%E9%BB%98%E8%AE%A4%E5%8F%AA%E8%BF%87%E6%BB%A43%E7%AB%AF%E7%9A%84%E4%BD%8E%E8%B4%A8%E9%87%8F%E5%BA%8F%E5%88%97%EF%BC%8C,%E5%A6%82%E6%9E%9C%E6%83%B3%E8%A6%81%E8%BF%87%E6%BB%A45%E7%AB%AF%E4%BD%8E%E8%B4%A8%E9%87%8F%E5%BA%8F%E5%88%97%EF%BC%8C%E9%9C%80%E8%A6%81%E7%94%A8%E9%80%97%E5%8F%B7%E9%9A%94%E5%BC%80%20cutadapt%20-q%2015%2C10%20-o%20output.fastq%20input.fastq%205%E7%AB%AF%E7%94%A815%E8%BF%9B%E8%A1%8C%E8%BF%87%E6%BB%A4%EF%BC%8C3%E7%AB%AF%E7%94%A810%E8%BF%9B%E8%A1%8C%E8%BF%87%E6%BB%A4)  
[参考](https://www.jianshu.com/p/4ee2f4d2292f)  
## 常用参数
```bash
-g: #剪切reads 5'端adapter(双端测序第一条read)，加$表示adapter锚定在reads 5'端
-a: #剪切reads 3'端adapter(双端测序第一条read)，加$表示adapter锚定在reads3'端
-O MINLENGTH, --overlap=MINLENGTH #adapter与reads最小overlap,才算成功识别; Default: 3
-m LENGTH, --minimum-length=LENGTH: 根据最短长度筛选reads;Default: 0
--discard-untrimmed, --trimmed-only #丢掉不包含adapter的reads
 -e ERROR_RATE, --error-rate=ERROR_RATE  #adapter匹配允许的最大错配率（错配/匹配片段长度)；Default: 0.1
--no-trim: 不修剪adapter，直接输出满足跳进啊的reads
-u LENGTH, --cut=LENGTH:  #修剪reads 5'/3'端碱基,正数：从开始除移除碱基；负数：从末尾处移除碱基；
-q [5'CUTOFF,]3'CUTOFF, --quality-cutoff=[5'CUTOFF,]3'CUTOFF: #修剪低质量碱基
-l LENGTH, --length=LENGTH: #将reads修剪的最终长度
--trim-n: #修剪reads末端的'N'
-o FILE, --output=FILE: #输出文件
--info-file=FILE：每条reads和与reads匹配的adapter的信息
--too-short-output=FILE: #为reads长度最小值设定阈值筛选reads后，要丢弃的部分输出到文件；长度依据m值设定；   
--too-long-output=FILE：#为reads长度最大值设定阈值筛选reads后，要丢弃的部分输出到文件；长度依据M值设定； 
--untrimmed-output=FILE: #将没有adapter未做修剪的reads输出到一个文件;默认输出到trimmed reads结果文件
--max-n=COUNT：#reads中N的数量，设定整数或小数(N的占比)

双端测序参数
-A ADAPTER：  #第二条reads 3'adapter
-G ADAPTER：#第二条reads 5'adapter
-U LENGTH： #从第二条reads上修剪的长度
-p FILE, --paired-output=FILE： #第二条reads的输出结果
--untrimmed-paired-output=FILE：#第一条reads没有adapter时，将第二条reads输出到文件；默认输出到trimmed reads结果文件   
```
## 参数详解
```bash 
$ cutadapt --help
Options:
  --version           
  -h, --help            
  --debug               Print debugging information.
  -f FORMAT, --format=FORMAT  #文件格式
                        Input file format：'fasta', 'fastq' or 'sra-fastq'.
                        Ignored when reading csfasta/qual files.
                        Default: auto-detect from file name extension.

  Finding adapters::
    Parameters -a, -g, -b specify adapters to be removed from each read
    (or from the first read in a pair if data is paired). If specified
    multiple times, only the best matching adapter is trimmed (but see the
    --times option). When the special notation 'file:FILE' is used,
    adapter sequences are read from the given FASTA file.

    -a ADAPTER, --adapter=ADAPTER #剪切reads 3'端adapter，加$表示congreads 3'端第一个碱基（无法adapter部分匹配识别）
                        Sequence of an adapter ligated to the 3' end (paired
                        data: of the first read). The adapter and subsequent
                        bases are trimmed. If a '$' character is appended
                        ('anchoring'), the adapter is only found if it is a
                        suffix of the read.匹配adapter
    -g ADAPTER, --front=ADAPTER #剪切reads 5'端adapter,加$表示congreads 5'端第一个碱基匹配adapter（无法adapter部分匹配识别）
                        Sequence of an adapter ligated to the 5' end (paired
                        data: of the first read). The adapter and any
                        preceding bases are trimmed. Partial matches at the 5'
                        end are allowed. If a '^' character is prepended
                        ('anchoring'), the adapter is only found if it is a
                        prefix of the read.
    -b ADAPTER, --anywhere=ADAPTER #adapter在 5'和3'端都可能出现时使用，慎用
                        Sequence of an adapter that may be ligated to the 5'
                        or 3' end (paired data: of the first read). Both types
                        of matches as described under -a und -g are allowed.
                        If the first base of the read is part of the match,
                        the behavior is as with -g, otherwise as with -a. This
                        option is mostly for rescuing failed library preparations
                        -do not use if you know which end your adapter was ligated to!       
    -e ERROR_RATE, --error-rate=ERROR_RATE  #adapter匹配允许的最大错配率（错配/匹配片段长度）
                        Maximum allowed error rate (no. of errors divided by
                        the length of the matching region). Default: 0.1
    --no-indels         Allow only mismatches in alignments. Default: allow
                        both mismatches and indels 禁止adapter发生Insertions和deletions
    -n COUNT, --times=COUNT #从reads行修剪adapter次数
                        Remove up to COUNT adapters from each read. Default: 1
    -O MINLENGTH, --overlap=MINLENGTH #adapter与reads最小overlap,才算成功识别;Default: 3
                        If the overlap between the read and the adapter is
                        shorter than MINLENGTH, the read is not modified.
                        Reduces the no. of bases trimmed due to random adapter
                        matches. Default: 3
    --match-read-wildcards
                        Interpret IUPAC wildcards in reads. Default: False
    -N, --no-match-adapter-wildcards
                        Do not interpret IUPAC wildcards in adapters.
    --no-trim           Match and redirect reads to output/untrimmed-output as
                        usual, but do not remove adapters. #不修剪adapter输出reads
    --mask-adapter      Mask adapters with 'N' characters instead of trimming
                        them. #识别adapter后，用'N'代替adapter

  Additional read modifications:
    -u LENGTH, --cut=LENGTH  #修剪reads 5'/3'端碱基,正数：从开始除移除碱基；负数：从末尾处移除碱基；
                        Remove bases from each read (first read only if
                        paired). If LENGTH is positive, remove bases from the
                        beginning. If LENGTH is negative, remove bases from
                        the end. Can be used twice if LENGTHs have different
                        signs.
    --nextseq-trim=3'CUTOFF
                        NextSeq-specific quality trimming (each read). Trims
                        also dark cycles appearing as high-quality G bases
                        (EXPERIMENTAL).
    -q [5'CUTOFF,]3'CUTOFF, --quality-cutoff=[5'CUTOFF,]3'CUTOFF #修剪低质量碱基
                        Trim low-quality bases from 5' and/or 3' ends of each
                        read before adapter removal. Applied to both reads if
                        data is paired. If one value is given, only the 3' end
                        is trimmed. If two comma-separated cutoffs are given,
                        the 5' end is trimmed with the first cutoff, the 3'
                        end with the second.
    --quality-base=QUALITY_BASE
                        Assume that quality values in FASTQ are encoded as
                        ascii(quality + QUALITY_BASE). This needs to be set to
                        64 for some old Illumina FASTQ files. Default: 33
    -l LENGTH, --length=LENGTH #将reads修剪的最终长度
                        Shorten reads to LENGTH. This and the following
                        modificationsare applied after adapter trimming. 
    --trim-n            Trim N's on ends of reads. #修剪reads末端的'N'
    --length-tag=TAG    Search for TAG followed by a decimal number in the
                        description field of the read. Replace the decimal
                        number with the correct length of the trimmed read.
                        For example, use --length-tag 'length=' to correct
                        fields like 'length=123'. #修剪reads后，用此参数修改reads的长度
    --strip-suffix=STRIP_SUFFIX
                        Remove this suffix from read names if present. Can be
                        given multiple times.
    -x PREFIX, --prefix=PREFIX #为reads名添加前缀，使用{name}添加adapter名字
                        Add this prefix to read names. Use {name} to insert
                        the name of the matching adapter.
    -y SUFFIX, --suffix=SUFFIX #为reads名添加后缀，使用{name}添加adapter名字
                        Add this suffix to read names; can also include {name}

  Filtering of processed reads:
    -m LENGTH, --minimum-length=LENGTH #修建前后，reads短于最短长度时，丢弃这对reads
                        Discard trimmed reads that are shorter than LENGTH.
                        Reads that are too short even before adapter removal
                        are also discarded. In colorspace, an initial primer
                        is not counted. Default: 0
    -M LENGTH, --maximum-length=LENGTH #修建前后，reads长于最大长度时，丢弃这对reads
                        Discard trimmed reads that are longer than LENGTH.
                        Reads that are too long even before adapter removal
                        are also discarded. In colorspace, an initial primer
                        is not counted. Default: no limit
     --max-n=COUNT       Discard reads with too many N bases. If COUNT is an
                        integer, it is treated as the absolute number of N
                        bases. If it is between 0 and 1, it is treated as the
                        proportion of N's allowed in a read. #reads中N的数量，设定整数或小数(N的占比)
    --discard-trimmed, --discard #丢弃只有一个adapter的reads
                        Discard reads that contain an adapter. Also use -O to
                        avoid discarding too many randomly matching reads!
    --discard-untrimmed, --trimmed-only #丢掉不包含adapter的reads
                        Discard reads that do not contain the adapter.  

  Output:
    --quiet             Print only error messages.
    -o FILE, --output=FILE #输出文件
                        Write trimmed reads to FILE. FASTQ or FASTA format is
                        chosen depending on input. The summary report is sent
                        to standard output. Use '{name}' in FILE to
                        demultiplex reads into multiple files. Default: write
                        to standard output
    --info-file=FILE    Write information about each read and its adapter
                        matches into FILE. #每条reads与adapter的匹配信息
    -r FILE, --rest-file=FILE
                        When the adapter matches in the middle of a read,
                        write the rest (after the adapter) into FILE.
    --wildcard-file=FILE
                        When the adapter has N bases (wildcards), write
                        adapter bases matching wildcard positions to FILE.
                        When there are indels in the alignment, this will
                        often not be accurate.
    --too-short-output=FILE #为reads长度最大值设定阈值筛选reads后，要丢弃的部分输出到文件；
                        Write reads that are too short (according to length
                        specified by -m) to FILE. Default: discard reads
    --too-long-output=FILE #将依据最长长度筛选reads后，要丢弃的部分输出到文件；
                        Write reads that are too long (according to length
                        specified by -M) to FILE. Default: discard reads
    --untrimmed-output=FILE #将没有adapter未做修剪的reads输出到一个文件，而不是输出到trimmed reads结果文件
                        Write reads that do not contain any adapter to FILE.
                        Default: output to same file as trimmed reads

  Colorspace options:
    -c, --colorspace    Enable colorspace mode: Also trim the color that is
                        adjacent to the found adapter.
    -d, --double-encode
                        Double-encode colors (map 0,1,2,3,4 to A,C,G,T,N).
    -t, --trim-primer   Trim primer base and the first color (which is the
                        transition to the first nucleotide)
    --strip-f3          Strip the _F3 suffix of read names
    --maq, --bwa        MAQ- and BWA-compatible colorspace output. This
                        enables -c, -d, -t, --strip-f3 and -y '/1'.
    --no-zero-cap       Do not change negative quality values to zero in
                        colorspace data. By default, they are since many tools
                        have problems with negative qualities.
    -z, --zero-cap      Change negative quality values to zero. This is
                        enabled by default when -c/--colorspace is also
                        enabled. Use the above option to disable it.

  Paired-end options:  #双端测序参数
    The -A/-G/-B/-U options work like their -a/-b/-g/-u counterparts, but
    are applied to the second read in each pair.

    -A ADAPTER          3' adapter to be removed from second read in a pair. #第二条reads 3'adapter
    -G ADAPTER          5' adapter to be removed from second read in a pair. #第二条reads 5'adapter
    -B ADAPTER          5'/3 adapter to be removed from second read in a pair.#第二条reads 5'/3'adapter
    -U LENGTH           Remove LENGTH bases from second read in a pair #从第二条reads上修剪的长度
    -p FILE, --paired-output=FILE #第二条reads的输出结果
                        Write second read in a pair to FILE. 
    --pair-filter=(any|both)
                        Which of the reads in a paired-end read have to match
                        the filtering criterion in order for it to be
                        filtered. Default: any
    --interleaved       Read and write interleaved paired-end reads.
    --untrimmed-paired-output=FILE #第一条reads没有adapter时，将第二条reads输出到文件；默认输出到trimmed reads结果文件 
                        Write second read in a pair to this FILE when no
                        adapter was found in the first read. Use this option
                        together with --untrimmed-output when trimming paired-
                        end reads. Default: output to same file as trimmed
                        reads
    --too-short-paired-output=FILE #将reads2中太短的reads输出到文件
                        Write second read in a pair to this file if pair is
                        too short. Use together with --too-short-output.
    --too-long-paired-output=FILE #将reads2中太长的reads输出到文件
                        Write second read in a pair to this file if pair is
                        too long. Use together with --too-long-output.

```
# Trimmomatic
## [使用方法](https://zhuanlan.zhihu.com/p/91691632)

