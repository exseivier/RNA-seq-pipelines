#!/usr/bin/env bash

<<HEADER_COMMENT
Pipeline Trapnell_et_al_2012

Author: Javier Montalvo
Contact: javier.montalvo@cinvestav.mx
RNA-Lab (Laboratory 12)
Langebio Cinvestav
Irapuato Guanajuato

This script trims and maps reads to a reference genome using hisat2
program, then it assembles the reads into a transcriptome and
measures the mRNA abundance by transcript with sammtools and cufflinks.
This version (2.1) of the script supports the batch processing of fastq
files

Pipeline summary
For every pair of reads: Fwd and Rev;
do
	hisat2 -> samtools[ SAM<=>BAM; sorting; indexing] -> cufflinks
done
HEADER_COMMENT

echo ""
echo "OUTPUT DATA"
echo ""
# User input mapping results path
<<ERASE
echo -n "Mapping results path, [Press ENTER if default]?: "; read MAPOUT
if [ "$MAPOUT" == "" ]; then
	MAPOUT="hisat2_out"
	export MAPOUT=$MAPOUT
	[ -d $MAPOUT ] || mkdir $MAPOUT
fi
echo "I got $MAPOUT"

# User input assembly results path
echo -n "Assembly results path, [Press ENTER if default]?: "; read ASSEMOUT
if [ "$ASSEMOUT" == "" ]; then
	ASSEMOUT="cl_out"
	export ASSEMOUT=$ASSEMOUT
	[ -d $ASSEMOUT ] || mkdir $ASSEMOUT
fi
echo "I got $ASSEMOUT"
ERASE

#User input for question, is there a reference genome index?
echo -n "Is there a reference genome index [yes|no]?: "; read IS_THERE_A_INDEX
if [ "$IS_THERE_A_INDEX" == "" ]; then
	IS_THERE_A_INDEX="no"
fi
echo "$IS_THERE_A_INDEX, I have a reference genome"

echo -n "Labels?: "; read LABELS
[ $LABELS ] || LABELS="q1,q2"
echo $LABELS
#####################
#      HISAT2       #
#        &&         #
#     CUFFLINKS     #
#####################
cd $PRWD
job=""

if [ "$IS_THERE_A_INDEX" == "no" ]; then
	job=$(echo "cd $PRWD; module load hisat2/2.0.4; hisat2-build $GENOMES/$REFGENOME $INDEX/$REFGENOME; sleep 10" | qsub -N indexing -l nodes=1:ppn=8,vmem=10gb,mem=10gb -V -q default)
fi
if [ "$job" == "" ]; then
	DEPENDENCY=""
else
	DEPENDENCY="-W depend=afterok:${job}"
fi

while IFS='' read -r LINE || [[ -n "$LINE" ]];
do
arrLINE=(${LINE//"|"/ })
SAMFILE_NAME=${arrLINE[0]}
if [ "$SAMFILE_NAME" == "" ]; then
	SAMFILE_NAME="output.sam"
fi
READS_F=${arrLINE[1]}
READS_R=${arrLINE[2]}
HISAT2_OPTIONAL_PARAMETERS=$(echo ${arrLINE[3]} | sed 's/;/ /g')
echo $HISAT2_OPTIONAL_PARAMETERS
CUFFLINKS_OPTIONAL_PARAMETERS=$(echo ${arrLINE[4]} | sed 's/;/ /g')
echo $HISAT2_OPTIONAL_PARAMETERS $CUFFLINKS_OPTIONAL_PARAMETERS "Hello"
if [ "$READS_F" == "" ] || [ "$READS_R" == "" ]; then
	echo "You have to specify the forward and reverse reads file"
	echo "Try again!"
	exit 1
fi
echo "
hisat2 -q -x $INDEX/$REFGENOME \
$HISAT2_OPTIONAL_PARAMETERS \
-1 $READS/$READS_F -2 $READS/$READS_R \
-S $MAPOUT/$SAMFILE_NAME"

echo "Name $SAMFILE_NAME; reads forward $READS_F; reads reverse $READS_R"

#<<MAIN_PROGRAM
hisat2_job=""
tjob=$(echo "
# Set working directory to HOME
cd $PRWD
# Loading HISAT2
# Loading samttols
# Loading cufflinks, cuffmerge and cuffdiff
# It is needed to change these settings depending on the type of cluster administration
# or depending on the module version you want to use.
module load hisat2/2.0.4
echo "hisat2 module loaded"
module load samtools/1.3.1
echo "samtools module loaded"
module load cufflinks/2.2.1
echo "cufflinks module loaded"
# Creating reference genome index. If index is not created
# Mapping reads to reference genome 
hisat2 -q -x $INDEX/$REFGENOME \
$HISAT2_OPTIONAL_PARAMETERS \
-1 $READS/$READS_F -2 $READS/$READS_R \
-S $MAPOUT/$SAMFILE_NAME
# SAM <=> BAM
samtools view -bT $GENOMES/$REFGENOME $MAPOUT/$SAMFILE_NAME > $MAPOUT/${SAMFILE_NAME%.*}.bam
# Sorting BAM
samtools sort $MAPOUT/${SAMFILE_NAME%.*}.bam -o $MAPOUT/${SAMFILE_NAME%.*}.sort.bam
# Indexing BAM
samtools index $MAPOUT/${SAMFILE_NAME%.*}.sort.bam
# Assembling genome and estimating RNA abundances
# Creating specific output directory
mkdir $ASSEMOUT/${SAMFILE_NAME%.*}
cufflinks -L ${SAMFILE_NAME%.*} $CUFFLINKS_OPTIONAL_PARAMETERS  -o $ASSEMOUT/${SAMFILE_NAME%.*}  $MAPOUT/${SAMFILE_NAME%.*}.sort.bam
echo "$ASSEMOUT/${SAMFILE_NAME%.*}/transcripts.gtf" >> $ASSEMOUT/GTFs.txt" | qsub -N $SAMFILE_NAME -l nodes=1:ppn=8,vmem=10gb,mem=10gb $DEPENDENCY -V -q default)
#MAIN_PROGRAM
if [ "$hisat2_job" == "" ]; then
	hisat2_job=$tjob
else
	hisat2_job=${hisat2_job}:${tjob}
fi
done < $COMMANDS/$MAPPING_ASSEMBLY_COMMANDS

##################
#    CUFFMERGE   #
#        &&      #
#     CUFFDIFF   #
##################

[ $CUFFMERGE_OUT ] || export CUFFMERGE_OUT="cm_out"
[ $CUFFDIFF_OUT ] || export CUFFDIFF_OUT="cd_out"
[ -d $CUFFMERGE_OUT ] || mkdir $CUFFMERGE_OUT
[ -d $CUFFDIFF_OUT ] || mkdir $CUFFDIFF_OUT

while IFS='' read -r LINE || [[ -n "$LINE" ]];
do
	arrLINE=(${LINE//"|"/ })
	SAMFILE_NAME=${arrLINE[0]}
	if [ "$SAMFILE_NAME" == "" ]; then
		SAMFILE_NAME="Assembling_transc"
	fi
	MAPPED_READS=$(echo ${arrLINE[1]} | sed 's/;/ /g')
	CUFFMERGE_OPTIONAL_PARAMETERS=$(echo ${arrLINE[2]} | sed 's/;/ /g')
	CUFFDIFF_OPTIONAL_PARAMETERS=$(echo ${arrLINE[3]} | sed 's/;/ /g')
	LABELS=${arrLINE[4]}
	echo "
module load cufflinks/2.2.1
cd $PRWD
mkdir $PRWD/$CUFFDIFF_OUT/${SAMFILE_NAME%.*}
mkdir $PRWD/$CUFFMERGE_OUT/${SAMFILE_NAME%.*}
cuffmerge -o $PRWD/$CUFFMERGE_OUT/${SAMFILE_NAME%.*} \
$CUFFMERGE_OPTIONAL_PARAMETERS \
$ASSEMOUT/GTFs.txt
cd $PRWD/$MAPOUT
cuffdiff -o $PRWD/$CUFFDIFF_OUT/${SAMFILE_NAME%.*} \
-L $LABELS \
$CUFFDIFF_OPTIONAL_PARAMETERS \
$PRWD/$CUFFMERGE_OUT/${SAMFILE_NAME%.*}/merged.gtf \
$MAPPED_READS" | \
qsub -N ${SAMFILE_NAME%.*}_ensam -l nodes=1:ppn=8,vmem=10gb,mem=10gb -W depend=afterok:$hisat2_job -V -q default
done < $COMMANDS/$DIFFEXP_COMMANDS

#<<EXPERIMENTAL_CODE
# Phase 0

#######################
#         SNP         #
#          &          #
#   Variant Calling   #
#######################

cd $PRWD
[ $VCF_OUT ] || export VCF_OUT="vcf_out"
[ -d $VCF_OUT ] || mkdir $VCF_OUT
while IFS='' read -r LINE || [[ -n "$LINE" ]];
do
	arrLINE=(${LINE//"|"/ })
	NAME=${arrLINE[0]}
	FILES=$(echo ${arrLINE[1]} | sed 's/;/ /g')
	VCF_OPTIONAL_PARAMETERS=$(echo ${arrLINE[2]} | sed 's/;/ /g')
	SNPVC_OP_PARAM=$(echo ${arrLINE[3]} | sed 's/;/ /g')
	echo "
cd $PRWD/$MAPOUT
module load samtools/1.3.1
module load bcftools/1.2

samtools mpileup -g -f $PRWD/$GENOMES/$REFGENOME \
$VCF_OPTIONAL_PARAMETERS \
$FILES \
> $PRWD/$VCF_OUT/${NAME}.bcf

bcftools call -O z -m -A -v $PRWD/$VCF_OUT/${NAME}.bcf \
$SNPVC_OP_PARAM \
-o $PRWD/$VCF_OUT/${NAME}_VC.vcf.gz" | \
qsub -N $NAME -l nodes=1:ppn=8,vmem=10gb,mem=10gb -W depend=afterok:$hisat2_job -V -q default

done < $PRWD/$COMMANDS/$VCF_COMMANDS
#EXPERIMENTAL_CODE
