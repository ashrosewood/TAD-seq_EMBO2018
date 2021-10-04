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
bwDir=$( pwd | awk '{print $1"/bigWigs"}' )
bigbed_dir=$( pwd | awk '{print $1"/bigWigs"}' )
ASSEMBLY="TFbb_12072016_linear"
DIR=$(basename $(pwd))
OUTFRAME=0

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
    echo >&2 "
$(basename $0) 
- Make a hub for 0+ samples for TAD-seq postive and negative selected samples bedgraph files,
- Specific for usage with TAD-seq data with a custom 2bit UCSC genome

AUTHOR: Ashley Woodfin
                
USAGE: $(basename $0) -b <BigWig directory> [OPTIONS]
 -b     BigWig directory                                         [default: $bwDir]
 -g     Linear Backbone Assembly (12072016 or 10102016)          [default: $ASSEMBLY]
 -h     Hub prefix                                               [default: $DIR]
 -p     Peak directory                                           [default: $bigbed_dir]
 -o     Make other out of frame open reading frame tracks (0/1)  [default: $OUTFRAME]

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
hub_directory=/groups/stark/www/data/ucsc/tadseq/${ASSEMBLY}"Hub"_${DIR}
track_directory=/groups/stark/www/data/ucsc/tadseq/${ASSEMBLY}"Hub"_${DIR}/${ASSEMBLY}

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

priority=1

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

#####################
## add luciferase files
#####################
for x in 4xUAS Nhe2; do
    priority=$(echo $priority | awk '{print $1+1}')
    (   echo "track Luciferase_${x}"
                    echo "bigDataUrl Luciferase_${x}.bb"
                    echo "shortLabel Luciferase_${x}"
                    echo "longLabel Luciferase_${x}"
                    echo "type bigBed 9"
                    echo "visibility pack"
                    echo "itemRgb on"
                    echo "priority $priority"
                    echo
                    echo -e "###\n"   ) >> ${track_directory}/trackDb.txt

ln -s /groups/stark/woodfin/analysis/tadseq/annotation/Luciferase_${x}.bb ${track_directory}/Luciferase_${x}.bb
done

##########################
## Make tracks_db.txt file
##########################

## set the number of postive selections based on the total number
## since the samples are named pos pos2 historically

EXPs=`ls -1 $bwDir/*neg_0+.bw | sed 's/_neg_0+.bw//g' | sed "s/.*\///"`
for EXP in $EXPs; do
    PosLimit=`grep $EXP experiments.txt | wc -l | awk '{print $1-1}'`
    if [ $PosLimit -gt 1 ];then
       Pos=$( seq 2 1 $PosLimit | awk '{print "pos", "pos"$1}' | tr '\n' ' ' )
    elif [ $PosLimit -gt 0 ];then
	Pos="pos"
    else
	Pos=""
    fi

    #############################
    ## in-frame neg, pos samples
    #############################
    for x in neg $Pos; do
	priority=$(echo $priority | awk '{print $1+1}')
	if [ $x == "neg" ]; then
	    color="216,191,216"
	elif [ $x == "pos" ]; then
	    color="186,85,211"
	elif [ $x == "pos2" ]; then
	    color="139,0,139"
	fi
	  (  echo "track ${EXP}_${x}_0+.bw"
             echo "bigDataUrl ${EXP}_${x}_0+.bw"
             echo "shortLabel ${EXP}_${x}_0+"
             echo "longLabel ${EXP}_${x}_0+"
             echo "type bigWig"
	     echo "group tad"
             echo "visibility full"
             echo "autoScale on"
             echo "maxHeightPixels 128:64:32"
             echo "color $color"
             echo "priority $priority"
             echo
             echo -e "###\n" )>> ${track_directory}/trackDb.txt
	### link the inframe
	  ln -s ${bwDir}/${EXP}_${x}_0+.bw ${track_directory}/${EXP}_${x}_0+.bw
    done
    
    ########################
    ## bc-col frags outframe
    ########################
    if [ $OUTFRAME -eq 1 ]; then
	priority=$(echo $priority | awk '{print $1+1}')
	(   echo "track container_${EXP}_outframe"
            echo "container multiWig"
            echo "shortLabel ${EXP}_outframe"
            echo "longLabel  ${EXP}_outframe"
            echo "type bigWig"
            echo "configurable on"
            echo "visibility hide"
            echo "autoScale on"
            echo "maxHeightPixels 128:64:32"
            echo "aggregate transparentOverlay"
            echo "showSubtrackColorOnUi on"
            echo "windowingFunction mean"
            echo "smoothingWindow off"
            echo "priority $priority"
            echo "" ) >> ${track_directory}/trackDb.txt
	    for x in neg $Pos; do
		for y in 0- 1+ 1- 2+ 2-; do
		    if [ $x == "neg" ] ; then
			color="255,204,153"
		    elif [ $x == "pos" ]; then
			color="255,178,102"
		    elif [ $x == "pos2" ]; then
			color="255,128,0"
		    fi
		    ( echo "track ${EXP}_${x}_${y}.bw"
                      echo "bigDataUrl ${EXP}_${x}_${y}.bw"
                      echo "shortLabel ${EXP}_${x}_${y}"
                      echo "longLabel ${EXP}_${x}_${y}"
                      echo "type bigWig"
		      echo "parent container_${EXP}_outframe"
                      echo "color $color"
		      echo "priority $priority"
                      echo ""
                    ) >> ${track_directory}/trackDb.txt
		    ln -s ${bwDir}/${EXP}_${x}_${y}.bw ${track_directory}/${EXP}_${x}_${y}.bw
		done
	    done
    ##print separator
    echo -e "###\n" >> ${track_directory}/trackDb.txt
    fi
    
    #########
    ## peaks
    #########
    for x in $Pos; do
	priority=$(echo $priority | awk '{print $1+1}')
	if [ -f ${bigbed_dir}/${EXP}_${x}_0+.peaks.bb ]; then
	    (   echo "track ${EXP}_${x}_0+.peaks"
                echo "bigDataUrl ${EXP}_${x}_0+.peaks.bb"
                echo "shortLabel ${EXP}_${x}_0+.peaks"
                echo "longLabel ${EXP}_${x}_0+.peaks"
                echo "type bigBed 9"
                echo "visibility dense"
                echo "itemRgb on"
                echo "useScore 0"
                echo "priority $priority"
                echo
                echo -e "###\n"   ) >> ${track_directory}/trackDb.txt
            
            #generate symbolic links to original data
            ln -s ${bigbed_dir}/${EXP}_${x}_0+.peaks.bb ${track_directory}/${EXP}_${x}_0+.peaks.bb
                
	fi
    done
    
    #########################
    ## P-vals and enrichments
    #########################
    for x in $Pos; do
	for y in enr DNminusUpPvl; do
	    if [[ $y =~ "enr" ]]; then 
		color="50,0,0"
		Views="-3:3"
            else
		color="0,0,255"
		Views="-5:5"
            fi    
            priority=$(echo $priority | awk '{print $1+1}')
            ( echo "track ${EXP}_${x}_0+.${y}.bw"
              echo "bigDataUrl ${EXP}_${x}_0+.${y}.bw"
              echo "shortLabel ${EXP}_${x}_0+.${y}"
              echo "longLabel ${EXP}_${x}_0+.${y}"
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
            ln -s ${bwDir}/${EXP}_${x}_0+.${y}.bw ${track_directory}/${EXP}_${x}_0+.${y}.bw
	done
    done
done

#####################
## add conservtion score
#####################
priority=$(echo $priority | awk '{print $1+1}')
( echo "track sequenceConservation.bw"
              echo "bigDataUrl sequenceConservation.bw"
              echo "shortLabel sequenceConservation"
              echo "longLabel sequenceConservation.bw"
              echo "type bigWig"
              echo "group tad"
              echo "visibility full"
              echo "autoScale off"
              echo "viewLimits 0:5.5"
              echo "maxHeightPixels 128:64:32"
              echo "color 128,128,128"
              echo "priority $priority"
              echo
              echo -e "###\n") >> ${track_directory}/trackDb.txt

ln -s /groups/stark/woodfin/analysis/tadseq/annotation/sequenceConservation.bw ${track_directory}/sequenceConservation.bw

echo "links are found here:"
echo $track_directory 

echo "here are the hubname address(es):"
echo "http://stark.imp.ac.at/data/ucsc/tadseq/${ASSEMBLY}"Hub"_${DIR}/hub.txt"


##############################################################################################
# exit
##############################################################################################

exit 0
