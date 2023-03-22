# Nature code
##  Calculate Cross-correlation QC scores: function _xcor()
```
CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.cc.qc"
CC_PLOT_FILE="${SUBSAMPLED_TA_FILE}.cc.plot.pdf"

# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag

Rscript $(which run_spp.R) -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE}

sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
mv temp ${CC_SCORES_FILE}
```  
![qccode](./pictures/QC_code.png)  



##  Generate self-pseudoreplicates for each replicate (SE datasets): function _spr()

```bash
# ========================
# Create pseudoReplicates
# =======================
PR_PREFIX="${OFPREFIX}.filt.nodup"
PR1_TA_FILE="${PR_PREFIX}.SE.pr1.tagAlign.gz"
PR2_TA_FILE="${PR_PREFIX}.SE.pr2.tagAlign.gz"

# Get total number of read pairs
nlines=$( zcat ${FINAL_TA_FILE} | wc -l )
nlines=$(( (nlines + 1) / 2 ))

# Shuffle and split BED file into 2 equal parts
zcat ${FINAL_TA_FILE} | shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f ${FINAL_TA_FILE} | wc -c) -nosalt </dev/zero 2>/dev/null) | split -d -l ${nlines} - ${PR_PREFIX} # Will produce ${PR_PREFIX}00 and ${PR_PREFIX}01

# Convert reads into standard tagAlign file
gzip -nc “${PR_PREFIX}00" > ${PR1_TA_FILE}
rm "${PR_PREFIX}00"
gzip -nc “${PR_PREFIX}01" > ${PR2_TA_FILE}
rm "${PR_PREFIX}01"
```

## Generate self-pseudoreplicates for each replicate (PE datasets): function _spr_PE()

```bash
# ========================
# Create pseudoReplicates
# =======================
PR_PREFIX="${OFPREFIX}.filt.nodup"
PR1_TA_FILE="${PR_PREFIX}.PE2SE.pr1.tagAlign.gz"
PR2_TA_FILE="${PR_PREFIX}.PE2SE.pr2.tagAlign.gz"
joined=”temp.bedpe”

# Make temporary fake BEDPE file from FINAL_TA_FILE
zcat ${FINAL_TA_FILE} | sed 'N;s/\n/\t/' | gzip -nc > $joined
# Get total number of read pairs
nlines=$( zcat ${joined} | wc -l )
nlines=$(( (nlines + 1) / 2 ))

# Shuffle and split BEDPE file into 2 equal parts
zcat -f ${joined} | shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f ${FINAL_TA_FILE} | wc -c) -nosalt </dev/zero 2>/dev/null) | split -d -l ${nlines} - ${PR_PREFIX} # Will produce ${PR_PREFIX}00 and ${PR_PREFIX}01

# Convert fake BEDPE  to reads into standard tagAlign file
awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "${PR_PREFIX}00" | gzip -nc > ${PR1_TA_FILE}
rm "${PR_PREFIX}00"
awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "${PR_PREFIX}01" | gzip -nc > ${PR2_TA_FILE}
rm "${PR_PREFIX}01"
rm $joined
```  

![self-pseudoreplicates](./pictures/self-pseudoreplicates.png)  

## Generate pooled dataset and pooled-pseudoreplicates : function _ppr()  

```bash
# ========================
# Create pooled datasets
# =======================
REP1_TA_FILE=”${DATASET_PREFIX}.Rep1.tagAlign.gz”
REP2_TA_FILE=”${DATASET_PREFIX}.Rep2.tagAlign.gz”
POOLED_TA_FILE=”${DATASET_PREFIX}.Rep0.tagAlign.gz”

zcat ${REP1_TA_FILE} ${REP2_TA_FILE} | gzip -nc > ${POOLED_TA_FILE}

# ========================
# Create pooled pseudoreplicates
# =======================
REP1_PR1_TA_FILE=”${DATASET_PREFIX}.Rep1.pr1.tagAlign.gz”
REP1_PR2_TA_FILE=”${DATASET_PREFIX}.Rep1.pr2.tagAlign.gz”

REP2_PR1_TA_FILE=”${DATASET_PREFIX}.Rep2.pr1.tagAlign.gz”
REP2_PR2_TA_FILE=”${DATASET_PREFIX}.Rep2.pr2.tagAlign.gz”

PPR1_TA_FILE=”${DATASET_PREFIX}.Rep0.pr1.tagAlign.gz”
PPR2_TA_FILE=”${DATASET_PREFIX}.Rep0.pr1.tagAlign.gz”

zcat ${REP1_PR1_TA_FILE} ${REP2_PR1_TA_FILE} | gzip -nc > ${PPR1_TA_FILE}
zcat ${REP1_PR2_TA_FILE} ${REP2_PR2_TA_FILE} | gzip -nc > ${PPR2_TA_FILE}
#先形成了各自的假重复，再合并为一个pool
```
![pooled-pseudoreplicates](./pictures/pooled-pseudoreplicates.png)  

## Tn5_shift