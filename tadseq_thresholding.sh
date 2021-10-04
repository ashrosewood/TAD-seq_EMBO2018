#!/bin/bash
set -o errexit
set -o pipefail

################################################################################
# Requirements
################################################################################

module load kent-ucsc/3.8f6f5e0a1cb75
module load bedtools/2.19.1

################################################################################
# Set default values
################################################################################
tempDir=/tmp               # change the tmp directory
ASSEMBLY="TFbb_10102015"   # back bone assembly
tableDir="tables"          # path to all.txt with p-values and log2FCs
peakDir="bigWigs"          # where to write tad bigBeds in linear coords for UCSC browser viewing
POS=""                     # name of the pos sample to evaluate ex Gal4-TFmix_12000ng_4xUAS_Lib14_merged_pos2_0+
FLANK=2100                 # length of flanking bb seq = start of TF
MINCOV=10                  # min coverage in GFP+ cells
MINLEN=60                  # min length of region that passes stringent threshold (region confined to TF sequence)
STRINGENT=1E-5             # stringent cutoff for TAD calling
LENIENT=1E-3               # lenient cutoff that determines TAD width
FC=1.5                     # fold-enrichment cutoff for TAD calling


################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
            echo >&2 "
$(basename $0) 
- Call TADs based from hypergeometric p-values, already calculated
for each position, requiring a lenient pval cut-off to overlap
a strict p-value cut-off that is at least as long as the minimal 
fragment length that overlaps the strict region.

- Specific for usage with TAD-seq data

AUTHOR: Ashley Woodfin
                
USAGE: $(basename $0) -i <BAM input file> -o <output files prefix> -b <barcodes> [OPTIONS]
 -P     Experiment name for GFP_pos sample (should end with _pos#) [required]
 -l     lenient p-value cutoff                                     [default: $p1]
 -s     strict p-value cutoff                                      [default: $p2]
 -f     fold change cutoff                                         [default: $FC]
 -c     min coverage in GFP+ cells                                 [default: $MINCOV]
 -L     min length of region that passes stringent threshold       [default: $MINLEN]
          (region confined to TF sequence) 
 -g     Genome assembly (e.g. which TAD-seq TF-backbone index)     [default: $ASSEMBLY]
 -t     Table dirctory with p-values and enrichments (all.txt)     [default: $tableDir]
 -o     Output directory for tads for browser (peaks.bb)           [default: $peakDir]
 -d     Mapped reads bed directoy                                  [default: $dataDir]

"
	    exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "P:l:s:f:c:L:g:t:o:d:" o
do
    case "$o" in
	P) POS="$OPTARG";;
	l) p1="$OPTARG";;
	s) p2="$OPTARG";;
	f) FC="$OPTARG";;
	c) MINCOV="$OPTARG";;
	L) MINLEN="$OPTARG";;
	g) ASSEMBLY="$OPTARG";;
	t) tableDir="$OPTARG";;
	o) peakDir="$OPTARG";;
	d) dataDir="$OPTARG";;
	\?) exit 1;;
    esac
done

if [ -z "$POS" ]; then
    echo >&2 "ERROR: -P is required!"; exit 1
fi


## set temp dir
TMP=$(mktemp -d -p "${tempDir}")


################################################################################
# Set chrom sizes $CHROM
# Throw error message if either index or chromosome sizes file does not exist
################################################################################

offset=/groups/stark/woodfin/analysis/tadseq/annotation/${ASSEMBLY}_linearOffset.txt

#######################################
# Call TADs and correct for FDR
#######################################

## get raw data table and multiple-testing correct it (use correct reading frame 0+)
cat $tableDir/${POS}.all.txt | /groups/stark/woodfin/starklab_repos/tad-seq/adjustPval.R -i - -c 13 > $tableDir/${POS}.all.bh.txt

## convert FC to log2
FC=$(awk -vFC=$FC 'BEGIN{print log(FC)/log(2)}')

## get stringent TADs with MINLEN and MINLEN within TF sequence (after FLANK)
cat $tableDir/${POS}.all.bh.txt | awk -v OFS='\t' -vPV=$STRINGENT -vFC=$FC -vMC=$MINCOV '($13<=PV && $8>=FC && $4>=MC){print $1, $2, $3}' | \
    bedtools merge | awk -vM=$MINLEN -vFL=$FLANK '($3-$2>=M && $3-M-FL>=0)' > $TMP/stringentTADcandidates

## get lenient TADs for merging
cat $tableDir/${POS}.all.bh.txt | awk -v OFS='\t' -vPV=$LENIENT '($13<=PV && $8>0){print $1, $2, $3}' | \
    bedtools merge > $TMP/lenientTADcandidates

## get lenient (longer?) TADs that overlap with stringent TADs
bedtools intersect -u -a $TMP/lenientTADcandidates -b $TMP/stringentTADcandidates > $TMP/tads

## get p-val and enr for these tads
bedtools intersect -wo -a $TMP/tads -b $tableDir/${POS}.all.bh.txt | \
    sort -k16,16g -k11,11gr | awk '!x[$1":"$2":"$3]++' > $TMP/tad_calls.txt

## call out the number of TADs
echo -ne "Number of TADs found: "; cat $TMP/tad_calls.txt | wc -l
## look at table and add TAD lengths
cat $TMP/tad_calls.txt | awk -vOFS='\t' '{print $1,$2,$3,$3-$2+1,$11,$16,$18}' | column -t > $tableDir/${POS}_tad_calls.txt

## make bigBed for linear genome
cat $tableDir/${POS}_tad_calls.txt | \
            awk -vFIX=$offset -vOFS='\t' ' 
           BEGIN{while((getline<FIX)>0)
            {
             TFlen[$1]=$2
            }}
            {if($1 in TFlen)
               {print "pvlSort", $2+TFlen[$1], $3+TFlen[$1], $1, 0, ".", $2+TFlen[$1], $3+TFlen[$1], "50,50,50"
               }
            }' | bedSort stdin stdout > $TMP/linear.bed

## save the tads on linear genome with bigWigs for hub
bedToBigBed $TMP/linear.bed /groups/stark/genomes/chrom/${ASSEMBLY}_linear.chrom.sizes $peakDir/${POS}.peaks.bb

##############################################################################################
# exit
##############################################################################################

rm -rf $TMP

exit 0
