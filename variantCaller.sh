#!/usr/bin/env bash

###################@
#      HISAT2,     @
#     SAMTOOLS     @
#         &        @
#     BCFTOOLS     @
###################@


echo -n "Is there a reference genome index [yes|no] ?: "; read IS_THERE_AN_INDEX
if [ "$IS_THERE_AN_INDEX" == "" ]; then
	IS_THERE_AN_INDEX="no"
fi

echo -n "Are the reads mapped to a reference genome [yes|no] ?: "; read ARE_READS_MAPPED
[ $ARE_READS_MAPPED ] || ARE_READS_MAPPED="no"

cd $PRWD
job=""

# If there is not an index, run an hisat2-build
if [ "$IS_THERE_A_INDEX" == "no" ]; then
	job=$(echo "cd $PRWD; module load hisat2/2.0.4; hisat2-build $GENOMES/$REFGENOME $INDEX/$REFGENOME; sleep 10" | qsub -N indexing -l nodes=1:ppn=8,vmem=10gb,mem=10gb -V -q default)
fi

# $DEPENDENCY value allocating
if [ "$job" == "" ]; then
	DEPENDENCY=""
else
	DEPENDENCY="-W depend=afterok:${job}"
fi
# Processing files
if [ "$ARE_READS_MAPPED" == "no" ]; then

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

###########################################################################
<<DEBUG
echo "
hisat2 -q -x $INDEX/$REFGENOME \
$HISAT2_OPTIONAL_PARAMETERS \
-1 $READS/$READS_F -2 $READS/$READS_R \
-S $MAPOUT/$SAMFILE_NAME"

echo "Name $SAMFILE_NAME; reads forward $READS_F; reads reverse $READS_R"
DEBUG
###########################################################################

hisat2_job=""
tjob=$(echo "
# Set working directory to HOME
cd $PRWD
# Loading HISAT2
# Loading samttols
# It is needed to change these settings depending on the type of cluster administration
# or depending on the module version you want to use.
module load hisat2/2.0.4
echo "hisat2 module loaded"
module load samtools/1.3.1
echo "samtools module loaded"
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
samtools index $MAPOUT/${SAMFILE_NAME%.*}.sort.bam" | qsub -N $SAMFILE_NAME -l nodes=1:ppn=8,vmem=10gb,mem=10gb $DEPENDENCY -V -q default)

if [ "$hisat2_job" == "" ]; then
	hisat2_job=$tjob
else
	hisat2_job=${hisat2_job}:${tjob}
fi
done < $COMMANDS/$MAPPING_ASSEMBLY_COMMANDS

fi


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

