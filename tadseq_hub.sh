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
ASSEMBLY="TFbb_12072016_linear"
bwDir=""
DIR=$(basename $(pwd))

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
    echo >&2 "
$(basename $0) - make a hub for 0+ samples for GFP positive samples
               - Specific for usage with TAD-seq data
                
USAGE: $(basename $0) -i <BAM input file> -o <output files prefix> -b <barcodes> [OPTIONS]
 -b     BigWig directory [required]
 -g     Genome assembly (e.g. dm3, hg19) [default: $ASSEMBLY]
 -h     Hub prefix [default: $DIR]

NOTES:
The program will use 1 core. 
"
    exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "b:g:h:" o
do
    case "$o" in
        b) bwDir="$OPTARG";;
        g) ASSEMBLY="$OPTARG";;
        h) DIR="$OPTARG";;
       \?) exit 1;;
    esac
done

if [ -z "$bwDir" -o -z "$DIR" ]; then
    echo >&2 "ERROR: -b is required!"; exit 1
fi

################################################################################
# Set up hub directory
################################################################################
hub_directory=/groups/stark/www/woodfin/${ASSEMBLY}"Hub"_${DIR}
track_directory=/groups/stark/www/woodfin/${ASSEMBLY}"Hub"_${DIR}/${ASSEMBLY}

if [ ! -d "$track_directory" ]; then
    mkdir -p $track_directory
fi

# add genomes.txt to hub directory
echo -e "genome ${ASSEMBLY}
trackDb ${ASSEMBLY}/trackDb.txt
description Sebastian TF Repressor mix linear assembly
twoBitPath ${ASSEMBLY}/${ASSEMBLY}.2bit
organism fly
defaultPos pvlSort:2500-10000
orderKey 4800
scientificName D. Melanogaster TF_REmix
htmlPath ${ASSEMBLY}/description.html" > $hub_directory/genomes.txt


## add hub.txt to the hub directory
echo -e "hub ${ASSEMBLY}"Hub"_${DIR}
shortLabel ${ASSEMBLY}"Hub"_${DIR}
longLabel TF and Repressor mix TADseq data
genomesFile genomes.txt
email ashley.woodfin@imp.ac.at" > $hub_directory/hub.txt

## add 2bit 
if [ ! -f "$track_directory/${ASSEMBLY}.2bit" ]; then
    cp /groups/stark/woodfin/analysis/tadseq/annotation/${ASSEMBLY}.2bit $track_directory/${ASSEMBLY}.2bit
fi

## add html
if [ ! -f "$track_directory/description.html" ]; then
    cp /groups/stark/woodfin/analysis/tadseq/annotation/${ASSEMBLY}/description.html $track_directory/description.html
fi

################################################################################
# sym link the bigwigs for 0+
################################################################################

Files=`ls -1 $bwDir/*bw | grep -E "enr|minus" | grep 0+`

if [ -f "${track_directory}/trackDb.txt" ]; then
    rm ${track_directory}/trackDb.txt
fi

if [ -f "${track_directory}/*bw" ]; then
    rm ${track_directory}/*bw
fi

###########################
### trackDb for enr and pvl
###########################
for File in $Files; do
    base=`basename $File .bw`
    BASE=`basename $File`
    color="50,0,0"
    Views="-5:5"
    
    ln -s -f ${File} ${track_directory}/${BASE}
    
    FRAME=$(echo $base | awk '{split($0,a,"."); print a[2]}')
    
    if [[ $FRAME =~ "minus" ]]
    then
        color="0,0,255"
        Views="-10:10"
    fi
    
    if [[ $FRAME =~ "0+" ]]
    then
        Vis="full"
    else
        Vis="hide"
    fi
    
    echo "processing:" ${File} $current_species "as" "color" $color
   
    echo -e "track ${base}.bw \nbigDataUrl ${base}.bw \nshortLabel $base \nlongLabel ${base} \ntype bigWig \ngroup tad \nvisibility $Vis \nautoScale off \nviewLimits $Views \nmaxHeightPixels 128:64:32 \ncolor ${color}\n" >> ${track_directory}/trackDb.txt
done

if [ -f "${track_directory}/peaks.txt" ]; then
    rm -f ${track_directory}/peaks.txt
fi

rm -f ${track_directory}/*bb

if [ ! -f "peaks/TF_RE_anno.bb" ]; then
    cp /groups/stark/woodfin/analysis/tadseq/annotation/${ASSEMBLY}_anno.bb peaks
fi
    
bigbeds=$( ls -1 peaks/*.bb )
for File in $bigbeds
do
    dir=${ASSEMBLY}"Hub"_${DIR}/${ASSEMBLY}
    FULL=`readlink -f $File`
    name=`basename $File .bb`
    BASE=`basename $File`
    ln -s -f ${FULL} ${track_directory}/${BASE}
    if [[ $FULL =~ "anno.bb" ]]
    then
        echo "processing:" ${File}
        echo -e "track type=\"bigBed 6\" name=$name description="${name}" itemRGB=On bigDataUrl=http://stark.imp.ac.at/woodfin/${dir}/${BASE}\n" >> ${track_directory}/peaks.txt
    else
	dir=${ASSEMBLY}"Hub"_${DIR}/${ASSEMBLY}
        echo "processing:" ${File}
        echo -e "track type=\"bigBed 9\" name=$name description="${name}" itemRGB=On bigDataUrl=http://stark.imp.ac.at/woodfin/${dir}/${BASE}\n" >> ${track_directory}/peaks.txt
    fi
done


echo -e "Hub directory is $hub_directory"

##############################################################################################
# exit
##############################################################################################

exit 0
