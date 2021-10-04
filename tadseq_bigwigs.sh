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
POS=""
NEG=""
ASSEMBLY="TFbb_12072016_linear"
tablesDir="tables"
bedDir="data"
outDir=$( pwd | awk '{print $1"/bigWigs"}' )
tempDir="/tmp/"
p="1e-5"


################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
    echo >&2 "
$(basename $0) - Make bigWig files for TAD-seq postive and negative selected samples bedgraph files,
                    - positve over negative enrichments at all positions, 
                    - and hyper geometic p-values at all positions 
                      log10(depleted) - log10(enriched), enriched p-values (+ value) and depleted (- value)
               - Specific for usage with TAD-seq data with a custom 2bit UCSC genome

AUTHOR: Ashley Woodfin
                
USAGE: $(basename $0) -P <GFP_pos> -N <GFP_neg> [OPTIONS]
 -P     Experiment name for GFP_pos sample (should end with _pos#) [required]
 -N     Experiment name for GFP_neg sample                         [required]
 -p     p-value threshold for enrichment files                     [default: $p]
 -g     Linear Backbone Assembly (12072016 or 10102016)            [default: $ASSEMBLY]
 -t     Directory with enrichments and pvalue tables               [default: $tablesDir]
 -b     Directory with TF bedgraph files                           [default: $bedDir]
 -o     Output files directory                                     [default: $outDir]

NOTES:
The program will use 2 cores (set by SGE). 
"
    exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "P:N:p:g:t:b:o:" o
do
    case "$o" in
        P) POS="$OPTARG";;
	N) NEG="$OPTARG";;
        p) p="$OPTARG";;
        g) ASSEMBLY="$OPTARG";;
	t) tablesDir="$OPTARG";;
	b) bedDir="$OPTARG";;
        o) outDir="$OPTARG";;
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

CHROM=/groups/stark/genomes/chrom/${ASSEMBLY}.chrom.sizes
offset=/groups/stark/woodfin/analysis/tadseq/annotation/${ASSEMBLY}Offset.txt

# does file exist
[ -e "$CHROM" ] || echo >&2 "ERROR: No chromosome size file found for genome assembly ${ASSEMBLY}!"


################################################################################
# Make bigWigs files for UCSC session
################################################################################
for r in 0 1 2; do for s in '+' '-'; do
   cat $tablesDir/${POS}_"$r$s.all.txt" | awk -vOFS='\t' -vp=$p '{if($8>0 && $13<p || $14<p && $8<0){print $1,$2,$3,$8}else{print $1,$2,$3,0}}' |\
     awk -vFIX=$offset -vOFS='\t' ' 
     BEGIN{while((getline<FIX)>0)
      {
       TFlen[$1]=$2
      }}
      {if($1 in TFlen)
         {print "pvlSort", $2+TFlen[$1], $3+TFlen[$1], $4
         }
      }' | sort -k1,1 -k2,2n | wigToBigWig -clip stdin $CHROM $TMP/${POS}_$r$s.enr.bw
   
   cat $tablesDir/${POS}_"$r$s.all.txt" | awk -vOFS='\t' -vp=$p '{if($13<=0){lpEn=-300}else{lpEn=log($13)/log(10)}; if($14<=0){lpDe=-300}else{lpDe=log($14)/log(10)}; {print $1,$2,$3,lpEn, lpDe}}' |\
     awk -vOFS='\t' '{print $1, $2, $3, $5-$4}' |\
     awk -vFIX=$offset -vOFS='\t' ' 
     BEGIN{while((getline<FIX)>0)
      {
       TFlen[$1]=$2
      }}
      {if($1 in TFlen)
         {print "pvlSort", $2+TFlen[$1], $3+TFlen[$1], $4
         }
      }'| sort -k1,1 -k2,2n | wigToBigWig -clip stdin $CHROM $TMP/${POS}_$r$s.DNminusUpPvl.bw
   
   cat $bedDir/${POS}_"$r$s.bg" | awk -vFIX=$offset -vOFS='\t' ' 
     BEGIN{while((getline<FIX)>0)
      {
       TFlen[$1]=$2
      }}
      {if($1 in TFlen)
         {print "pvlSort", $2+TFlen[$1], $3+TFlen[$1], $4
         }
      }'| sort -k1,1 -k2,2n | wigToBigWig -clip stdin $CHROM $TMP/${POS}_$r$s.bw

   cat $bedDir/${NEG}_"$r$s.bg" | awk -vFIX=$offset -vOFS='\t' ' 
     BEGIN{while((getline<FIX)>0)
      {
       TFlen[$1]=$2
      }}
      {if($1 in TFlen)
         {print "pvlSort", $2+TFlen[$1], $3+TFlen[$1], $4
         }
      }'| sort -k1,1 -k2,2n | wigToBigWig -clip stdin $CHROM $TMP/${NEG}_$r$s.bw
  done
done


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
