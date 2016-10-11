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


#User input for question, is there a reference genome index?
echo -n "Is there a reference genome index [yes|no]?: "; read IS_THERE_A_INDEX
if [ "$IS_THERE_A_INDEX" == "" ]; then
	IS_THERE_A_INDEX="no"
fi
echo "$IS_THERE_A_INDEX, I have a reference genome"

echo ""
echo "PBS SETTINGS"
echo ""
#User input number of nodes
echo -n "Number of nodes, [Press ENTER if default (1)]?: "; read NUMBER_OF_NODES
if [ "$NUMBER_OF_NODES" == "" ]; then
	NUMBER_OF_NODES=1
fi
echo "I got $NUMBER_OF_NODES nodes"

#User input processors per node
echo -n "Processors per node, [Press ENTER if default (1)]?: "; read PROCESSOR_PER_NODE
if [ "$PROCESSOR_PER_NODE" == "" ]; then
	PROCESSOR_PER_NODE=1
fi
echo "I got $PROCESSOR_PER_NODE processors per node"

#User input queue
echo -n "Queue name, [Press ENTER if (default)]?: "; read QUEUE
if [ "$QUEUE" == "" ]; then
	QUEUE="default"
fi
echo "$QUEUE is the name of the queue"

#User input virtual memory
echo -n "Virtual memory amount in gb, [Press ENTER if default (5gb)]?: "; read VIRTUAL_MEMORY
if [ "$VIRTUAL_MEMORY" == "" ]; then
	VIRTUAL_MEMORY=5
fi
echo "I got ${VIRTUAL_MEMORY}gb of virtual memory"

#User input memory
echo -n "Memory amount in gb, [Press ENTER if default (5gb)]?: "; read MEMORY
if [ "$MEMORY" == "" ]; then
	MEMORY=5
fi
echo "I got ${MEMORY}gb of memory"

#<<DEBUG2
echo "echo "qsub -N $PIPELINE_NAME -l nodes=$NUMBER_OF_NODES:ppn=$PROCESSOR_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb -V -q $QUEUE"" | bash
echo "qsub -N $PIPELINE_NAME -l nodes=$NUMBER_OF_NODES:ppn=$PROCESSOR_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb -V -q $QUEUE"
#DEBUG2


while IFS='' read -r LINE || [[ -n "$LINE" ]]; do
arrLINE=(${LINE//"|"/ })
SAMFILE_NAME=${arrLINE[0]}
if [ "$SAMFILE_NAME" == "" ]; then
	SAMFILE_NAME="output.sam"
fi
READS_F=${arrLINE[1]}
READS_R=${arrLINE[2]}
HISAT2_OPTIONAL_PARAMETERS=${arrLINE[3]}
CUFFLINKS_OPTIONAL_PARAMETERS=${arrLINE[4]}
if [ "$READS_F" == "" ] || [ "$READS_R" == "" ]; then
	echo "You have to specify the forward and reverse reads file"
	echo "Try again!"
	exit 1
fi
echo "Name $SAMFILE_NAME; reads forward $READS_F; reads reverse $READS_R"
cd $PRWD

#<<MAIN_PROGRAM
echo "
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
if [ "$IS_THERE_A_INDEX" == "no" ]; then
	hisat2-build $GENOMES/$REFGENOME $INDEX/$REFGENOME
fi

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
echo "$ASSEMOUT/${SAMFILE_NAME%.*}/transcripts.gtf" >> $ASSEMOUT/GTFs.txt" | qsub -N $SAMFILE_NAME -l nodes=$NUMBER_OF_NODES:ppn=$PROCESSOR_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb -V -q $QUEUE
#MAIN_PROGRAM
done < $COMMANDS/$MAPPING_ASSEMBLY_COMMANDS
#__EOF__//~~CAT
