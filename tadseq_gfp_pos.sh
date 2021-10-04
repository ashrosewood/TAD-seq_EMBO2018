#!/bin/bash
set -o errexit
set -o pipefail

################################################################################
# Requirements
################################################################################

module load samtools/0.1.18
module load bedtools/2.19.1

################################################################################
# Set default values
################################################################################
ASSEMBLY="TFbb_12072016"
POS=""
NEG=""
outDir=$( pwd | awk '{print $1"/tables"}' )
dataDir="data"
tempDir="/tmp/"
p="1e-5"
TFtables=0

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
    echo >&2 "
$(basename $0) - Calculate p-values and enrichments for TAD-seq positively selected samples over negatively selected samples
               - Specific for usage with TAD-seq data
                
USAGE: $(basename $0) -P <GFP_pos> -N <GFP_neg> [OPTIONS]
 -P     Experiment name for GFP_pos sample (should end with _pos<#GFP facs sorts>) [required]
 -N     Experiment name for GFP_neg sample (should end with _neg) [required]
 -p     p-value threshold for enrichment files                   [default: $p]
 -g     Genome assembly (e.g. dm3, hg19)                         [default: $ASSEMBLY]
 -d     Directory with GFP_pos and GFP_neg samples               [default: $dataDir]
 -o     Directory where output files will be put                 [default: $outDir]
 -T     Save a separate enrichment and p-value table for each TF [default: $TFtables]
NOTES:
The program will use 2 cores (set by SGE). 
"
    exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "P:N:p:g:d:o:T:" o
do
    case "$o" in
        P) POS="$OPTARG";;
        N) NEG="$OPTARG";;
        p) p="$OPTARG";;
        g) ASSEMBLY="$OPTARG";;
	d) dataDir="$OPTARG";;
        o) outDir="$OPTARG";;
	T) TFtables="$OPTARG";;
       \?) exit 1;;
    esac
done

if [ -z "$POS" -o -z "$NEG" ]; then
    echo >&2 "ERROR: -P -N are required!"; exit 1
fi


## set temp dir
TMP=$(mktemp -d -p "${tempDir}")

################################################################################
# Set processor usage
# Either assigned by SGE or half of the available cores
# (at least 4 as we anyway need 4 later on)
################################################################################

if [ -z "$NSLOTS" ]; then
    USE=$(grep processor /proc/cpuinfo | wc -l | awk '{if($1>6){print int($1/2)}else{print 4}}')
else
    USE=$NSLOTS
fi


################################################################################
# Set chrom sizes $CHROM
# Throw error message if either index or chromosome sizes file does not exist
################################################################################

CHROM=/groups/stark/genomes/chrom/${ASSEMBLY}.chrom.sizes


# does file exist
[ -e "$CHROM" ] || echo >&2 "ERROR: No chromosome size file found for genome assembly ${ASSEMBLY}!"

################################################################################
# run hyper geometric test for GFP+ and ++ vs GFP-
# All values saved in $TMP/${POS}_"$r$s.all.txt"
# P-value thresholded enrichments in $TMP/${POS}_"$r$s.enr.txt"
# -log10(p-value) tables in $TMP/${POS}_"$r$s.pvl.txt"
################################################################################

for r in 0 1 2; do for s in '+' '-'; do 
       ### get the total fragment count for pos and neg
       pos=$(cat $dataDir/${POS}_"$r$s.bed" | wc -l)
       neg=$(cat $dataDir/${NEG}_"$r$s.bed" | wc -l)
       unionBedGraphs -g $CHROM -i $dataDir/${POS}_"$r$s.bg" $dataDir/${NEG}_"$r$s.bg"  | \
            awk -vP=$pos -vN=$neg -vOFS='\t' 'BEGIN{F=N/P} {print $1,$2,$3,$4,$5,P,N,log(F*($4+1/F)/($5+1))/log(2),$4,$5+$4,P,N+P}' | \
            hyper.R -i - -m 9 -n 10 -M 11 -N 12 > $TMP/${POS}_"$r$s.all.txt"
       ### enr=$8, p-enr=$13, p-deplete=$14 make enrichment table where all enrichments are set to zero if with insignificant p-values 
       cat $TMP/${POS}_"$r$s.all.txt" | awk -vp=$p -vOFS='\t' '{if($8>0 && $13>p){r=0}else if($8<0 && $14>p){r=0}else{r=$8} for(i=$2;i<$3;i++){print $1,i,i+1,r}}' > $TMP/${POS}_"$r$s.enr.txt"
       ## make -log10(p-value) file and set lp to 100 if p = 0 to avoid inf 
       cat $TMP/${POS}_"$r$s.all.txt" | awk -vOFS='\t' '{if($13<=0){lp=100}else{lp=-log($13)/log(10)}; for(i=$2;i<$3;i++){print $1,i,i+1,lp}}' > $TMP/${POS}_"$r$s.pvl.txt"
    done
done

### make enrich and pval table per TF

if [ "$TFtables" -eq 1 ]; then
      mkdir -p $TMP/$POS
      ### grep by individual TFs and get the 
      awk -F'\t' '($8!=0){print $1}' $TMP/${POS}_0+.enr.txt | uniq | while read TF; do
          for r in 0 1 2; do for s in '+' '-'; do grep "^$TF" $TMP/${POS}_${r}${s}.enr.txt > $TMP/tf_$r$s.tmp; done; done
          paste $TMP/tf_{0,1,2}+.tmp $TMP/tf_{0,1,2}-.tmp | awk -vOFS='\t' '{print $1,$3,$4,$8,$12,$16,$20,$24}' > $TMP/$POS/$TF.enr.txt
      done


      awk -F'\t' '($8!=0){print $1}' $TMP/${POS}_0+.enr.txt | uniq | while read TF; do
      for r in 0 1 2; do for s in '+' '-'; do grep "^$TF" $TMP/${POS}_${r}${s}.pvl.txt > $TMP/tf_$r$s.tmp; done; done
           paste $TMP/tf_{0,1,2}+.tmp $TMP/tf_{0,1,2}-.tmp | awk -vOFS='\t' '{print $1,$3,$4,$8,$12,$16,$20,$24}' > $TMP/$POS/$TF.pvl.txt
      done
fi
      
##############################################################################################
# move all relavent files to desired output directory 
# keep everything except fa files from bowtie for the table
##############################################################################################

if [ ! -d "$outDir" ]; then
    mkdir $outDir
fi

mv $TMP/* $outDir


##############################################################################################
# exit
##############################################################################################

rm -rf $TMP

exit 0
