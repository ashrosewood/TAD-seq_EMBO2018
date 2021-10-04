#!/bin/bash
set -o errexit
set -o pipefail

################################################################################
# Requirements
################################################################################

module load R/3.2.2
module load samtools/0.1.18
module load bedtools/2.19.1
module load gridengine
bsub=/groups/stark/software-all/shell/bsub_gridengine

################################################################################
# Set default values
################################################################################
ASSEMBLY="TFbb_12072016_linear"
POS=""
NEG=""
CDS="2100"
p="0.001"
scriptdir="/groups/stark/woodfin/analysis/tadseq/pipeline"
tempDir="/clustertmp/stark/woodfin/"
addAA="9"

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
    echo >&2 "
$(basename $0) - Peak calling and sequence extraction
                - Specific for usage with TAD-seq data
                
USAGE: $(basename $0) -i <BAM input file> -o <output files prefix> -b <barcodes> [OPTIONS]
 -P     Experiment name for GFP_pos sample (should end with _pos#) [required]
 -N     Experiment name for GFP_neg sample [required]
 -p     p-value threshold for enrichment files [default: $p]
 -g     Genome assembly (e.g. dm3, hg19) [default: $ASSEMBLY]
 -C     Length of the backbone flanking the fragment [default: $CDS]
 -A     Number of amino acids to add to the ends of the peaks (e.g. 9 is 3AA) [default: $addAA]

NOTES:
The program will use 7 cores for Rscript (set by SGE). 
"
    exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "P:N:p:g:C:" o
do
    case "$o" in
        P) POS="$OPTARG";;
	N) NEG="$OPTARG";;
        p) p="$OPTARG";;
        g) ASSEMBLY="$OPTARG";;
	C) CDS="$OPTARG";;
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
# Convert beds to bigBeds in linear genome coordinates for peak calling
################################################################################
for EXP in $POS $NEG; do
    for f in 0 1 2; do for s in '+' '-'; do 
       cat data/${EXP}_"$f$s.bed" |
       awk -vFIX=$offset -vOFS='\t' ' 
       BEGIN{while((getline<FIX)>0)
        {
         TFlen[$1]=$2
        }}
        {if($1 in TFlen)
           {print "pvlSort", $2+TFlen[$1], $3+TFlen[$1], $4, 0, $6
           }
        }' | bedSort stdin stdout > $TMP/${EXP}_"$f$s.bed"
       bedToBigBed $TMP/${EXP}_"$f$s.bed" $CHROM $TMP/${EXP}_"$f$s.bb"     
      done;done
done

################################################################################
# Call peaks for tad-seq data
################################################################################
if [ ! -d peaks ]; then
    mkdir peaks
fi

for f in "0+" "0-" "1+" "1-" "2+" "2-"; do
    $bsub -o log_peaks -C 2-4 "/groups/stark/software-all/shell/starr.0.0.6/call_peaks.sh -e ${TMP}/${POS}_${f}.bb -b ${TMP}/${NEG}_${f}.bb -g $ASSEMBLY -o peaks/${POS}_${f} -P $p"
done

# Wait for peak jobs to be done
mytemp=$(mktemp)

ls log_peaks/tmp*  | xargs -n1 basename | grep -v '\([A-Za-z0-9]*\.\)\{2\}' | awk '{print substr($1,1,10)}' | sort > $mytemp
qstat | sed '1,2d' | awk '{print $3}' | sort > ${mytemp}_2

while grep -q -F -f ${mytemp} ${mytemp}_2; do
    qstat | sed '1,2d' | awk '{print $3}' | sort > ${mytemp}_2
    sleep 2
done

rm ${mytemp}*

################################################################################
# calculate 0+ fdrs based on other frames 
################################################################################
Rscript ${scriptdir}/tadseq_peaks_FDRs_fromOtherFrames.R --Cores=7 --Sample=$POS

################################################################################
# annotate peaks with TF name and position of peak in TF
################################################################################

tf_anno=/groups/stark/woodfin/analysis/tadseq/annotation/${ASSEMBLY}_anno.bb

# get actual TF and postion in linear genome for overlaps
bigBedToBed $tf_anno stdout | awk -vOFS='\t' -vC=$CDS '{$2=$2-C; $3=$3+C; print $0}' | bedSort stdin stdout > $TMP/tf.tmp


## get the TF name and position per peak
cat peaks/fdrPvals/${POS}_0+_peaksfdr.txt | bedSort stdin stdout | intersectBed -wao -a stdin  -b $TMP/tf.tmp | cut --fields=16-18,20-22 --complement |\
    awk -vFIX=$offset -vOFS='\t' ' 
         BEGIN{while((getline<FIX)>0)
          {
           TFlen[$1]=$2
          }}
          {if($16 in TFlen)
             {print $0, $2-TFlen[$16], $3-TFlen[$16], $5-TFlen[$16]
             }
          }' | \
     awk -vOFS='\t' -vPF=$CDS -vA=$addAA '
          {
           start=($17-PF)%3 
           end=($18-PF)%3
           $17=$17-start-A
           $18=$18-end+A
           print $0
           }' | sort -k11,11gr > peaks/fdrPvals/${POS}_"0+_peaksfdr_anno.txt"

##############################################################################################
# get the sequence for the peaks in frame
##############################################################################################
if [ ! -d peaks/sequence ]; then
    mkdir peaks/sequence
fi

nonLinear=`echo $ASSEMBLY | sed 's/_linear//g'`

FASTA=/groups/stark/woodfin/analysis/tadseq/annotation/${nonLinear}.fa

awk -vOFS='\t' '{print $16,$17,$18,$16}' peaks/fdrPvals/${POS}_0+_peaksfdr_anno.txt | \
    sort -k1,1 -k2,3n |
     bedtools getfasta -fi $FASTA -bed - -fo peaks/sequence/${POS}_0+_peaks_coordinates.fa -name


##############################################################################################
# exit
##############################################################################################

rm -rf $TMP

exit 0
