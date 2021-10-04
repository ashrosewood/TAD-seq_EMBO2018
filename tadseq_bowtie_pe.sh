#!/bin/bash
set -o errexit
set -o pipefail

################################################################################
# Requirements
################################################################################

module load samtools/0.1.18
module load kent-ucsc/3.8f6f5e0a1cb75
module load bowtie/0.12.9
module load bedtools/2.19.1

################################################################################
# Set default values
################################################################################
EXP=""
BAM=""
BARCODE=""
YLINK=""
ASSEMBLY="TFbb_12072016"
BACKBONE="tadSeqBackBone_12072016"
bcCol="12"                                                     # usually column 12 for MiSeq BZ2(col 12) 
BCMM="2"
BOWTIE_MM="3"
CDS="2100"
MULTofTHREE="N"
YMM="2"
R1_OFFSET="0"
TAGLEN="20"                                                    #this script MUST be modified if the molecular barcode and Ylinker are > 30 for paired-end 50
tempDir="/tmp/"
outDir=$( pwd | awk '{print $1"/data"}' )
FASTQ="0"
keepFilteredReads="0"
Y_FILT="1"
BARCODE_LEN=8 

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
    echo >&2 "
$(basename $0) - Extract barcoded reads from a paired-end BAM file and map the reads against a dm3 TF assembly
                - Specific for usage with TAD-seq data
                
USAGE: $(basename $0) -i <BAM input file> -o <output files prefix> -b <barcodes> [OPTIONS]
 -i     Input file (BAM) [required]
 -e     Experiment name for output files prefix [required]
 -b     Barcodes [required]
        Note: This script does NOT allow mismatches in BC
 -B     Number of the column in the bam file that contains the barcode [default: $bcCol]
 -g     Genome assembly (e.g. dm3, hg19) [default: $ASSEMBLY]
 -m     Maximum number of mismatches per indidual read during bowtie mapping [default: $BOWTIE_MM]
 -c     Maximum number of mismatches in Barcode for collapsing [default: $BCMM]
 -y     Y-linker sequence [required]
 -M     Multiple of three [required]
         * if Gal4_TF construct (N)
         * if TF_Gal4 construct (Y)         
 -Y     Maximum number of mismatches in Barcode for filtering [default: $YMM]
 -C     Length of the backbone flanking the fragment [default: $CDS]
 -r     Read 1 offset [default: $R1_OFFSET]
 -T     Directory to which temporary files will be written [default: $tempDir]
 -o     Directory where output files will be put [default: $outDir]
 -L     Read length (after Ylinker and Random BC excluded) to use for alignment [default: $TAGLEN]
 -F     Save gzipped fastq files and exit without mapping (0/1) [default: $FASTQ]
        Reads will be saved as EXP_R1.fastq.gz and EXP_R2.fastq.gz
 -f     Save gzipped fasta reads that contain ylinker and exit without mapping (0/1) [default: $keepFilteredReads]
        Reads will be saved as EXP_ylinker_read_1.fa.gz and EXP_ylinker_read_2.fa.gz
 -a    Backbone index prefix [default: $BACKBONE]
 -t    Filter reads for ylinker (0/1) [default: $Y_FILT]
 -l    Barcode length [default: $BARCODE_LEN] 

NOTES:
The program will use \$NSLOTS cores (set by SGE) or half of all cores on a machine
(minimum of 4 cores required). 
Note that if there are not 20 positions left in read 2 after the molecular barcode and Y-linker sequence the job will be killed."
    exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "i:e:b:B:g:m:c:y:M:Y:C:r:T:o:L:F:f:a:t:l:" o
do
    case "$o" in
        i) BAM="$OPTARG";;
        e) EXP="$OPTARG";;
        b) BARCODE="$OPTARG";;
        B) bcCol="$OPTARG";;
        g) ASSEMBLY="$OPTARG";;
        R) RAN_BARCODE_LEN="$OPTARG";;
        g) ASSEMBLY="$OPTARG";;
        m) BOWTIE_MM="$OPTARG";;
        c) BCMM="$OPTARG";;
        y) YLINK="$OPTARG";;
        M) MULTofTHREE="$OPTARG";;
        Y) YMM="$OPTARG";;
        C) CDS="$OPTARG";;
        r) R1_OFFSET="$OPTARG";;
        T) tempDir="$OPTARG";;
        o) outDir="$OPTARG";;
	L) TAGLEN="$OPTARG";;
	F) FASTQ="$OPTARG";;
	f) keepFilteredReads="$OPTARG";;
	a) BACKBONE="$OPTARG";;
	t) Y_FILT="$OPTARG";;
	l) BARCODE_LEN="$OPTARG";;
       \?) exit 1;;
    esac
done

if [ -z "$BAM" -o -z "$EXP" -o -z "$BARCODE" -o -z "$YLINK" -o -z "$MULTofTHREE" ]; then
    echo >&2 "ERROR: -i -e -b -y -M are required!"; exit 1
fi


## set params based on input
#BARCODE_LEN=$( echo $BARCODE | awk '{print length($1)}' )    #column alows up to 9, but we are not always using 9
TMP=$(mktemp -d -p "${tempDir}")
YLen=$( echo "$YLINK" | awk '{print length($1)}')

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
# Set $INDEX and $SIZES
# Throw error message if either index or chromosome sizes file does not exist
################################################################################

INDEX=/groups/stark/indices/bowtie/tadseq/${ASSEMBLY}
SIZES=/groups/stark/genomes/chrom/${ASSEMBLY}.chrom.sizes


# does file exist
[ -e "$SIZES" ] || echo >&2 "ERROR: No chromosome size file found for genome assembly ${ASSEMBLY}!"
[ -e "${INDEX}.1.ebwt" ] || echo >&2 "ERROR: No bowtie index files found for genome assembly ${ASSEMBLY}!"

################################################################################
# Extract barcoded reads
# Reads are saved in $TMP/$SAMPLE_R1.fq and $TMP/SAMPLE_R2.fq
################################################################################

samtools view $BAM | \
            awk -vT=${TMP} -vB=${BARCODE} -vS=${EXP} -vCOL=${bcCol} -vBCL=${BARCODE_LEN} -F'\t' '
            BEGIN{
                  n=split(B,Bc,"|");for(i=1;i<=n;i++){BCS[substr(Bc[i],1,BCL)]=1}
                  }
                 {
                  bc=substr($COL,6,BCL); if(bc in BCS){
                    {
                     print ($1"\n"$10"\n+\n"$11) > (T "/" S "_R1.fastq") 
                     getline
                     print ($1"\n"$10"\n+\n"$11) > (T "/" S "_R2.fastq")
                    }
                  }
                 }'

# Make outdir
if [ ! -d "$outDir" ]; then
    mkdir $outDir
fi

# Save FASTQ file
if [ "$FASTQ" = "1" ]
then
    cat $TMP/${EXP}_R1.fastq | gzip > $outDir/${EXP}_R1.fastq.gz
    cat $TMP/${EXP}_R2.fastq | gzip > $outDir/${EXP}_R2.fastq.gz
    rm -rf $TMP
    exit 0
else
    gzip $TMP/${EXP}_R1.fastq
    gzip $TMP/${EXP}_R2.fastq
fi


################################################################################
# Extract and trim tag sequences
# Filter for the presence of Y-linker sequence with up to 2 mismatches
# Reads are saved in $TMP/${EXP}_read_1.fa.gz and $TMP/${EXP}_read_1.fa.gz
################################################################################
if [ $Y_FILT = "1" ]
   then
       ## extract and trim tag sequences, filter for presence of Y-linker sequence
       paste <(zcat $TMP/${EXP}_R1.fastq.gz | awk -vS=${R1_OFFSET} '(NR%4==2){print substr($1,1+S)}') <(zcat $TMP/${EXP}_R2.fastq.gz | awk -vOFS='\t' -vL=$YLen '(NR%4==2){ST=11+L; print substr($1,ST), substr($1,1,10), substr($1,11,L)}') | \
	   awk -vMM=$YMM -vYL=$YLINK '{mm=0; for(i=1;i<=length(YL);i++){if(substr($4,i,1)!=substr(YL,i,1)){mm++}} if(mm<=MM){print $0}}' | \
	   awk -vL=$TAGLEN -vTMP=$TMP '{
            print (">" substr($1,1,L) "_" $3 "\n" substr($1,1,L)) > TMP"/read_1.fa"
            print (">" substr($2,1,L) "_" $3 "\n" substr($2,1,L)) > TMP"/read_2.fa"
        }'
else
    if [ $Y_FILT = "0" ]
    then
	paste <(zcat $TMP/${EXP}_R1.fastq.gz | awk -vS=${R1_OFFSET} '(NR%4==2){print substr($1,1+S)}') <(zcat $TMP/${EXP}_R2.fastq.gz | awk -vOFS='\t' -vL=$YLen '(NR%4==2){ST=11+L; print substr($1,ST), substr($1,1,10), substr($1,11,L)}') | \
            awk -vL=$TAGLEN -vTMP=$TMP '{
            print (">" substr($1,1,L) "_" $3 "\n" substr($1,1,L)) > TMP"/read_1.fa"
            print (">" substr($2,1,L) "_" $3 "\n" substr($2,1,L)) > TMP"/read_2.fa"
        }'
    fi
fi

## zip and save Ylinker filtered fa files
if [ "$keepFilteredReads" = "1" ]
then
    cat $TMP/read_1.fa | gzip > $outDir/${EXP}_ylinker_read_1.fa.gz
    cat $TMP/read_2.fa | gzip > $outDir/${EXP}_ylinker_read_2.fa.gz
    rm -rf $TMP
    exit 0
fi

##############################################################################################
# Align Y-linker filtered reads with bowtie
# All fragments saved in $TMP/${EXP}_frags.bed
# Collapse by position and molecular BC saved in $TMP/${EXP}_collapsed_frags.bed
# Further collapsed by molecular barcodes with mismatches > $BCMM saved in $TMP/${EXP}_all.bed
##############################################################################################

## save fa files for unaligned and multimapped for table
bowtie -v $BOWTIE_MM -m 1 --quiet --strata --best --al $TMP/${EXP}_aligned.fa --un $TMP/${EXP}_unaligned.fa --max $TMP/${EXP}_valid_alignments.fa -f -p $USE  $INDEX -X 2000 -1 $TMP/read_1.fa -2 $TMP/read_2.fa | \
    awk -F'\t' -vOFS='\t' '{
          if(NR%2==1){
             split($1,T,"/")
             R[T[2]]=T[1]
             S=((T[2]==1)?"+":"-")
             B=$4
          }else{
                split($1,T,"/")
                R[T[2]]=T[1]
                print $3,B,$4+length($5),R[1]"_"R[2],"1",S
          }}' | sort -k1,1 -k2,2n -k3,3nr -k6,6 -k5,5nr > $TMP/${EXP}_frags.bed


## position and molecular BC collapsing
## Awk grab molecular barcode
cat $TMP/${EXP}_frags.bed | \
    awk -F'\t' -vOFS='\t' '{BC=split($4,T,"_"); print $1, $2, $3, T[2], $5, $6}' | \
    sort -k1,1 -k2,2n -k3,3nr -k6,6 -k4,4 | \
    uniq -c | awk -vOFS='\t' '{print $2,$3,$4,$5,$1,$7}' | \
    sort -k1,1 -k2,2n -k3,3nr -k6,6 -k5,5nr > $TMP/${EXP}_collapsed_frags.bed


## similar barcodes collapsing (BCs that are have only up to <=BCMM mismatches are considered identical)
cat $TMP/${EXP}_collapsed_frags.bed | \
    awk -vMM=$BCMM '{key=($1":"$2":"$3":"$6); if(key!=old){delete(BC); old=key; BC[$4]=1; print}else{m=0;for(bc in BC){mm=0;for(i=1;i<=10;i++){if(substr(bc,i,1)!=substr($4,i,1)){mm++;if(mm>MM){break}}}if(mm<=MM){m=1;break}} if(m==0){print;BC[$4]=1}}}' > $TMP/${EXP}_all.bed


##############################################################################################
# split into the different reading frames, which depends on the length of the plasmid backbone
# saved in $TMP/${EXP}_<fragstrand>.bed
# calculate the coverage for each frame 
# saved in $TMP/${EXP}_<fragstrand>.bg
##############################################################################################

if [ $MULTofTHREE = "Y" ]; then
    cat $TMP/${EXP}_all.bed | awk -F'\t' -vOFS='\t' -vPF=$CDS -vT=$TMP -vE=$EXP '(($3-$2)%3==0){frame=($2-PF)%3; if(frame<0){frame=frame+3}; print $0 > (T "/" E "_" frame $6 ".bed")}'
else
    cat $TMP/${EXP}_all.bed | awk -F'\t' -vOFS='\t' -vPF=$CDS -vT=$TMP -vE=$EXP '{frame=($2-PF)%3; if(frame<0){frame=frame+3}; print $0 > (T "/" E "_" frame $6 ".bed")}'
fi

for f in 0 1 2; do for s in '+' '-'; do 
        cat $TMP/${EXP}_"$f$s.bed" | genomeCoverageBed -i stdin -bga -g $SIZES > $TMP/bedgraph_tmp
        ## fix genomeCoverageBed bug that does not report a chromosome (here: TF) at all if the TF does not have any read matching to it - even despite the use of "-bga"
        cat $SIZES | awk -vOFS='\t' -vBG=$TMP/bedgraph_tmp 'BEGIN{while(getline<BG){GOTIT[$1]=1; print}} (!($1 in GOTIT)){print $1, 0, $2, 0}' | sort -k1,1 -k2,2n > $TMP/${EXP}_"$f$s.bg"
done; done

#######################################################
# align unmapped reads to the backbone and contaminants
#######################################################

BBindex="/groups/stark/indices/bowtie/tadseq/${BACKBONE}"

bowtie -v $BOWTIE_MM -m 1 --quiet --strata --best --al $TMP/${EXP}_backBone_aligned.fa --un $TMP/${EXP}_backBone_unaligned.fa --max $TMP/${EXP}_backBone_valid_alignments.fa -f -p $USE  $BBindex -X 1000 -1 $TMP/${EXP}_unaligned_1.fa -2 $TMP/${EXP}_unaligned_2.fa| \
            awk -F'\t' -vOFS='\t' '{
              if(NR%2==1){
                 split($1,T,"/")
                 R[T[2]]=T[1]
                 S=((T[2]==1)?"+":"-")
                 B=$4
              }else{
                    split($1,T,"/")
                    R[T[2]]=T[1]
                    print $3,B,$4+length($5),R[1]"_"R[2],".",S
              }}' > $TMP/${EXP}_backBone_frags.bed

## contaminants
Cindex="/groups/stark/muerdter/genomes/contaminants/contaminants"

bowtie -v $BOWTIE_MM -m 1 --quiet --strata --best --al $TMP/${EXP}_cont_aligned.fa --un $TMP/${EXP}_cont_unaligned.fa --max $TMP/${EXP}_cont_valid_alignments.fa -f -p $USE  $Cindex -X 1000 -1 $TMP/${EXP}_backBone_unaligned_1.fa -2 $TMP/${EXP}_backBone_unaligned_2.fa| \
            awk -F'\t' -vOFS='\t' '{
              if(NR%2==1){
                 split($1,T,"/")
                 R[T[2]]=T[1]
                 S=((T[2]==1)?"+":"-")
                 B=$4
              }else{
                    split($1,T,"/")
                    R[T[2]]=T[1]
                    print $3,B,$4+length($5),R[1]"_"R[2],".",S
              }}' > $TMP/${EXP}_cont_frags.bed

##############################################################################################
# make a summary file for all reads, ylinker filtering, and alignment
##############################################################################################

## percentages relative to all reads
allreads=$( paste <(zcat $TMP/${EXP}_R1.fastq.gz | awk '(NR%4==2)') | wc -l )
Ylink=$( grep -v "^>" $TMP/read_2.fa | wc -l | awk -vA=$allreads '{print $1, $1/A*100}' )
YL=$( echo $Ylink | awk '{print $1}' )
noYlink=$( echo $Ylink | cut -f 1 | awk -vA=$allreads '{NoY=A-$1; print NoY, NoY/A*100}' )

aligned=$( cat $TMP/${EXP}_aligned_2.fa | grep -v "^>" | wc -l )
al=$( echo $aligned | awk -vA=$YL '{print $1, $1/A*100}'  )
mult=$( cat $TMP/${EXP}_valid_alignments_2.fa | grep -v "^>" | awk -vA=$YL 'BEGIN{i=0}{i++;}END{print i+0, (i+0)/A*100}' )
un=$( cat $TMP/${EXP}_unaligned_2.fa | grep -v "^>" | awk -vA=$YL 'BEGIN{i=0}{i++;}END{print i+0, (i+0)/A*100}' )


BC_col=$( cat $TMP/${EXP}_collapsed_frags.bed | wc -l | awk -vA=$YL '{print $1, $1/A*100}')
coor_col=$( cat $TMP/${EXP}_all.bed | wc -l | awk -vA=$YL '{print $1, $1/A*100}' )
inframe_col=$( cat $TMP/${EXP}_0+.bed | wc -l | awk -vA=$YL '{print $1, $1/A*100}' )

if [ -f $TMP/${EXP}_cont_valid_alignments_1.fa ]
then
    Cont_map=$( cat $TMP/${EXP}_cont_frags.bed $TMP/${EXP}_cont_valid_alignments_1.fa | grep -v "^>" | awk -vA=$YL 'BEGIN{i=0}{i++;}END{print i+0, (i+0)/A*100}' )
else
    Cont_map=$( cat $TMP/${EXP}_cont_frags.bed | awk -vA=$YL 'BEGIN{i=0}{i++;}END{print i+0, (i+0)/A*100}' )
fi

if [ -f $TMP/${EXP}_backBone_valid_alignments_1.fa ]
then
    BB_map=$( cat $TMP/${EXP}_backBone_frags.bed $TMP/${EXP}_backBone_valid_alignments_1.fa | grep -v "^>" | awk -vA=$YL 'BEGIN{i=0}{i++;}END{print i+0, (i+0)/A*100}' )
else
    BB_map=$( cat $TMP/${EXP}_backBone_frags.bed | awk -vA=$YL 'BEGIN{i=0}{i++;}END{print i+0, (i+0)/A*100}' )
fi

if [ -f $TMP/${EXP}_cont_unaligned_1.fa ]
then
    NONE=$( cat $TMP/${EXP}_cont_unaligned_1.fa | grep -v "^>" | awk -vA=$YL 'BEGIN{i=0}{i++;}END{print i+0, (i+0)/A*100}' )
else
    NONE=$( echo "0 0")
fi



(   echo "All_reads $allreads 100"
    echo "No_Ylinker $noYlink"
    echo "Ylinker $Ylink"
    echo "TF_align $al"
    echo "TF_multimapped $mult"
    echo "TF_unaligned $un"
    echo "TF_BC_col $BC_col"
    echo "TF_Coord_col $coor_col"
    echo "TF_inFramePos $inframe_col"
    echo "BB_align $BB_map"
    echo "Cont_align $Cont_map"
    echo "no_mapping $NONE"
) | tr ' ' '\t' >> $TMP/${EXP}_read_report.txt


##############################################################################################
# move all relavent files to desired output directory 
# keep everything except fa files from bowtie for the table
##############################################################################################

rm $TMP/*fa
rm $TMP/*gz
rm $TMP/*backBone*
rm $TMP/*cont*

mv $TMP/${EXP}_* $outDir

##############################################################################################
# exit
##############################################################################################

rm -rf $TMP

exit 0
