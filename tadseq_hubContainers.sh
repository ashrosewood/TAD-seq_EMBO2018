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
bigbed_dir=$( pwd | awk '{print $1"/peaks"}' )
OUTFRAME=0

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
 -p     Peak directory [default: $bigbed_dir]
 -o     Make out of frame containers (0/1) [default: $OUTFRAME]

NOTES:
The program will use 1 core. 
"
    exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "b:g:h:p:o:" o
do
    case "$o" in
        b) bwDir="$OPTARG";;
        g) ASSEMBLY="$OPTARG";;
        h) DIR="$OPTARG";;
	p) bigbed_dir="$OPTARG";;
	o) OUTFRAME="$OPTARG";;
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

if [ -d "$track_directory" ]; then
    rm -rf $track_directory
fi

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
# clean up directory if already exists
################################################################################

if [ -f "${track_directory}/trackDb.txt" ]; then
    rm ${track_directory}/trackDb.txt
fi

if [ -f "${track_directory}/*bw" ]; then
    rm ${track_directory}/*bw
fi

if [ -f "${track_directory}/*bb" ]; then
    rm ${track_directory}/*bb
fi

##################################################################
### trackDb for bc-collapsed frags, enrichments, pvals, and peaks
##################################################################

#set priority to 0
priority=0
EXPs=`ls -1 $bwDir/*pos*_0+.bw | sed 's/_0+.bw//g' | sed "s/.*\///"`

for EXP in $EXPs; do
    ## write track file first for overlayed 
    
    priority=$(echo $priority | awk '{print $1+1}')
    if [ $OUTFRAME -eq 1 ]; then
	################
	## bc-col frags
	################
	(   echo "track container_${EXP}_outframe"
            echo "container multiWig"
            echo "shortLabel ${EXP}_outframe"
            echo "longLabel  ${EXP}_outframe"
            echo "type bigWig"
            echo "configurable on"
            echo "visibility hide"
            echo "autoScale on"
            echo "viewLimits 0:1000"
            echo "maxHeightPixels 128:64:32"
            echo "color 153,51,255"
            echo "aggregate transparentOverlay"
            echo "showSubtrackColorOnUi on"
            echo "windowingFunction mean"
            echo "smoothingWindow off"
            echo "priority $priority"
            echo "") >> ${track_directory}/trackDb.txt
	
	## add the other frames to overlay if flag set
        for x in 1 2; do
            if [ $x -eq 1 ]; then
               colorPos="178,102,255"
	       colorNeg="255,178,102"
            else
	       colorPos="204,153,255"
	       colorNeg="255,204,153"
            fi
	    ( echo "track ${EXP}_${x}+.bw"
                    echo "bigDataUrl ${EXP}_${x}+.bw"
                    echo "shortLabel ${EXP}_${x}+"
                    echo "longLabel ${EXP}_${x}+"
                    echo "type bigWig"
                    echo "parent container_${EXP}_outframe"
                    echo "color $colorPos"
                    echo ""
                    echo "track ${EXP}_${x}-.bw"
                    echo "bigDataUrl ${EXP}_${x}-.bw"
                    echo "shortLabel ${EXP}_${x}-"
                    echo "longLabel ${EXP}_${x}-"
                    echo "type bigWig"
                    echo "parent container_${EXP}_outframe"
                    echo "color $colorNeg"
                    echo "" ) >> ${track_directory}/trackDb.txt
        done
    ##print separator
    echo -e "###\n" >> ${track_directory}/trackDb.txt

    ## add links to data
    for x in 1 2; do
	for s in + -; do
	    ln -s ${bwDir}/${EXP}_${x}${s}.bw ${track_directory}/${EXP}_${x}${s}.bw
	done
    done

    fi

    ### link the inframe
    for s in + -; do
       ln -s ${bwDir}/${EXP}_0${s}.bw ${track_directory}/${EXP}_0${s}.bw
    done
    
    
    #########
    ## peaks
    #########
    priority=$(echo $priority | awk '{print $1+1}')
    if [ -f ${bigbed_dir}/${EXP}_0+.peaks.bb ]; then
	(   echo "track ${EXP}_0+.peaks"
                    echo "bigDataUrl ${EXP}_0+.peaks.bb"
                    echo "shortLabel ${EXP}_0+.peaks"
                    echo "longLabel ${EXP}_0+.peaks"
                    echo "type bigBed 9"
                    echo "visibility dense"
                    echo "itemRgb on"
                    echo "useScore 0"
                    echo "priority $priority"
                    echo
                    echo -e "###\n"   ) >> ${track_directory}/trackDb.txt
                    
                #generate symbolic links to original data
                ln -s ${bigbed_dir}/${EXP}_0+.peaks.bb ${track_directory}/${EXP}_0+.peaks.bb
                
    fi
    
    #########################
    ## P-vals and enrichments
    #########################
    for x in enr DNminusUpPvl; do
	if [[ $x =~ "enr" ]]; then 
            color="50,0,0"
            Views="-3:3"
        else
            color="0,0,255"
            Views="-10:10"
        fi    
        priority=$(echo $priority | awk '{print $1+1}')
        ( echo "track ${EXP}_0+.${x}.bw"
          echo "bigDataUrl ${EXP}_0+.${x}.bw"
          echo "shortLabel ${EXP}_0+.${x}"
          echo "longLabel ${EXP}_0+.${x}"
          echo "type bigWig"
          echo "group tad"
          echo "visibility full"
          echo "autoScale off"
          echo "viewLimits $Views"
          echo "maxHeightPixels 128:64:32"
          echo "color $color"
          echo "priority $priority"
          echo
          echo -e "###\n") >> ${track_directory}/trackDb.txt

	#generate symbolic links to original data
        ln -s ${bwDir}/${EXP}_0+.${x}.bw ${track_directory}/${EXP}_0+.${x}.bw
    done
done

#####################
## add the anno file
#####################
(   echo "track ${ASSEMBLY}_anno"
                    echo "bigDataUrl ${ASSEMBLY}_anno.bb"
                    echo "shortLabel ${ASSEMBLY}_anno"
                    echo "longLabel ${ASSEMBLY}_anno"
                    echo "type bigBed 6"
                    echo "visibility pack"
                    echo "itemRgb on"
                    echo "priority $priority"
                    echo
                    echo -e "###\n"   ) >> ${track_directory}/trackDb.txt

ln -s /groups/stark/woodfin/analysis/tadseq/annotation/${ASSEMBLY}_anno.bb ${track_directory}/${ASSEMBLY}_anno.bb

############################
## add neg samples overlayed
############################
if [ $OUTFRAME -eq 1 ]; then
    NEGs=`ls -1 $bwDir/*neg*_0+.bw | sed 's/_0+.bw//g' | sed "s/.*\///"`
    
    for EXP in $NEGs; do
	priority=$(echo $priority | awk '{print $1+1}')
	################
	## bc-col frags
	################
	(   echo "track container_${EXP}_outframe"
            echo "container multiWig"
            echo "shortLabel ${EXP}"
            echo "longLabel  ${EXP}"
            echo "type bigWig"
            echo "configurable on"
            echo "visibility hide"
            echo "autoScale on"
            echo "viewLimits 0:1000"
            echo "maxHeightPixels 128:64:32"
            echo "color 153,51,255"
            echo "aggregate transparentOverlay"
            echo "showSubtrackColorOnUi on"
            echo "windowingFunction mean"
            echo "smoothingWindow off"
            echo "priority $priority"
            echo "") >> ${track_directory}/trackDb.txt
	
	## add the other frames to overlay
	for x in 1 2; do
	    if [ $x -eq 1 ]; then
		colorPos="178,102,255"
		colorNeg="255,178,102"
            else
		colorPos="204,153,255"
		colorNeg="255,204,153"
            fi
	    ( echo "track ${EXP}_${x}+.bw"
              echo "bigDataUrl ${EXP}_${x}+.bw"
              echo "shortLabel ${EXP}_${x}+"
              echo "longLabel ${EXP}_${x}+"
              echo "type bigWig"
              echo "parent container_${EXP}_outframe"
              echo "color $colorPos"
              echo ""
              echo "track ${EXP}_${x}-.bw"
              echo "bigDataUrl ${EXP}_${x}-.bw"
              echo "shortLabel ${EXP}_${x}-"
              echo "longLabel ${EXP}_${x}-"
              echo "type bigWig"
              echo "parent container_${EXP}_outframe"
              echo "color $colorNeg"
              echo "" ) >> ${track_directory}/trackDb.txt
	done
	for x in 1 2; do
	    for s in + -; do
		ln -s ${bwDir}/${EXP}_${x}${s}.bw ${track_directory}/${EXP}_${x}${s}.bw
	    done
	done
    done
fi

############################
## add sample container with -, + and ++
############################
SAMPLEs=`ls -1 $bwDir/*neg*_0+.bw | sed 's/_neg_0+.bw//g' | sed "s/.*\///"`

for EXP in $SAMPLEs; do
    priority=$(echo $priority | awk '{print $1+1}')
    ################
    ## bc-col frags
    ################
    (   echo "track container_${EXP}_0+"
        echo "container multiWig"
        echo "shortLabel ${EXP}"
        echo "longLabel  ${EXP}"
        echo "type bigWig"
        echo "configurable on"
        echo "visibility full"
        echo "autoScale on"
        echo "viewLimits 0:1000"
        echo "maxHeightPixels 128:64:32"
        echo "color 153,51,255"
        echo "aggregate transparentOverlay"
        echo "showSubtrackColorOnUi on"
        echo "windowingFunction mean"
        echo "smoothingWindow off"
        echo "priority $priority"
        echo ""
	if [ -d "tables/${EXP}_pos3" ]; then
	    echo "track ${EXP}_pos3_0+.bw"
            echo "bigDataUrl ${EXP}_pos3_0+.bw"
            echo "shortLabel ${EXP}_pos3_0+"
            echo "longLabel ${EXP}_pos3_0+"
            echo "type bigWig"
            echo "parent container_${EXP}_0+"
            echo "color 117,0,117"
            echo ""
	fi
        echo "track ${EXP}_pos2_0+.bw"
        echo "bigDataUrl ${EXP}_pos2_0+.bw"
        echo "shortLabel ${EXP}_pos2_0+"
        echo "longLabel ${EXP}_pos2_0+"
        echo "type bigWig"
        echo "parent container_${EXP}_0+"
        echo "color 139,0,139"
        echo ""
        echo "track ${EXP}_pos_0+.bw"
        echo "bigDataUrl ${EXP}_pos_0+.bw"
        echo "shortLabel ${EXP}_pos_0+"
        echo "longLabel ${EXP}_pos_0+"
        echo "type bigWig"
        echo "parent container_${EXP}_0+"
        echo "color 186,85,211"
        echo ""
        echo "track ${EXP}_neg_0+.bw"
        echo "bigDataUrl ${EXP}_neg_0+.bw"
        echo "shortLabel ${EXP}_neg_0+"
        echo "longLabel ${EXP}_neg_0+"
        echo "type bigWig"
        echo "parent container_${EXP}_0+"
        echo "color 216,191,216"
        echo "") >> ${track_directory}/trackDb.txt
    for s in + -; do
	ln -s ${bwDir}/${EXP}_neg_0${s}.bw ${track_directory}/${EXP}_neg_0${s}.bw
    done
done

echo "here are the hubname address(es):"
echo "http://stark.imp.ac.at/woodfin/${ASSEMBLY}"Hub"_${DIR}/hub.txt"


##############################################################################################
# exit
##############################################################################################

exit 0
