#!/usr/bin/env bash

<<HEADER_COMMENT
Pipeline Trapnell_et_al_2012

Author: Javier Montalvo
Contact: javier.montalvo@cinvestav.mx
RNA-Lab (Laboratory 12)
Langebio Cinvestav
Irapuato Guanajuato

This script maps reads to a reference genome using hisat2 program,
then it assembles the reads into a transcriptome and measures the
level of gene expression using sammtools and cufflinks.

Pipeline summary
hisat2 -> samtools[ SAM<=>BAM; sorting; indexing] -> cufflinks
HEADER_COMMENT

<<PBS_SETTINGS
#PBS -N <str>
#PBS -l nodes=<int>:ppn=<int>,vmem=<int>gb,mem=<int>gb
#PBS -V
#PBS -q <str>
PBS_SETTINGS


#User input pipeline name
echo -n "Pipeline name?: "; read PIPELINE_NAME
if [ "$PIPELINE_NAME" == "" ]; then
	PIPELINE_NAME="default"
fi
echo "I got $PIPELINE_NAME"

echo ""
echo "INPUT DATA"
echo ""
# User input HOME path
echo -n "Home path?: "; read CASA
if [ "$CASA" == "" ]; then
	echo "It is important to sepcify the HOME path!"
	echo "Try again!"
	exit 1
else
	echo "I got $CASA"
fi

# User input GENOMES path
echo -n "Genomes path?: "; read GENOMAS
if [ "$GENOMAS" == "" ]; then
	echo "You have to specify the GENOMES path!"
	echo "Try again!"
	exit 1
else
	echo "I got $GENOMAS"
fi

# User input INDEX path
echo -n "Index path?: "; read INDICE
if [ "$INDICE" == "" ]; then
	echo "You have to specify the INDEX path!"
	echo "Try again!"
	exit 1
else
	echo "I got $INDICE"
fi

# User input genoma fasta file name
echo -n "Genoma fasta file name?: "; read GENOMA_FASTA_FILE
if [ "$GENOMA_FASTA_FILE" == "" ]; then
	echo "I do not know what is the genome fasta file name!"
	echo "Try again!"
	exit 1
else
	echo "I got $GENOMA_FASTA_FILE"
fi

# User input READS path
echo -n "Reads path?: "; read LECTURAS
if [ "$LECTURAS" == "" ]; then
	echo "I do not know where the reads are!"
	echo "Would you like to tell me?"
	exit 1
else
	echo "I got $LECTURAS"
fi

# User input trimmed forward reads file name
echo -n "Trimmed forward reads file name?: "; read READS_F
# User input reverse reads file name
echo -n "Trimmed reverse reads file name?: "; read READS_R
if [ "$READS_F" == "" ] || [ "$READS_R" == "" ]; then
	echo "It is needed to specify the name of the reads files!"
	echo "try again!"
	exit 1
else
	echo "I got $READS_F"
	echo "I got $READS_R"
fi

echo ""
echo "OUTPUT DATA"
echo ""
# User input mapping results path
echo -n "Mapping results path?: "; read MAPOUT
if [ "$MAPOUT" == "" ]; then
	MAPOUT="hisat2_out"
fi
echo "I got $MAPOUT"

# User input assembly results path
echo -n "Assembly results path?: "; read ASSEMOUT
if [ "$ASSEMOUT" == "" ]; then
	ASSEMOUT="cl_out"
fi
echo "I got $ASSEMOUT"

# User input sam file name
echo -n "SAM file name?: "; read SAM_FILE
if [ "$SAM_FILE" == "" ]; then
	SAM_FILE="output.sam"
fi
echo "I got $SAM_FILE"

echo ""
echo "OPTIONAL PARAMETERS"
echo ""

#User input hisat2 optional parameters
echo -n "Hisat2 optional parameters, [Press ENTER if default]?: "; read HISAT2_OPTIONAL_PARAMETERS
if [ "$HISAT2_OPTIONAL_PARAMETERS" == "" ]; then
	echo "Hisat2 optional parameters were set as default!"
else
	echo "Hisat2 optional parameters are:"
	echo $HISAT2_OPTIONAL_PARAMETERS
fi

#User input cufflinks optional parameters
echo -n "Cufflinks optional parameters, [Press ENTER if default]?: "; read CUFFLINKS_OPTIONAL_PARAMETERS
if [ "$CUFFLINKS_OPTIONAL_PARAMETERS" == "" ]; then
	echo "Cufflinks optional parameters were set as default!"
else
	echo "Cufflinks optional parameters are:"
	echo $CUFFLINKS_OPTIONAL_PARAMETERS
fi

#User input LABEL
echo -n "Label for this work?: "; read LABEL
if [ "$LABEL" == "" ]; then
	LABEL="treatment"
fi
echo "I got $LABEL"

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
echo -n "Number of nodes?: "; read NUMBER_OF_NODES
if [ "$NUMBER_OF_NODES" == "" ]; then
	NUMBER_OF_NODES=1
fi
echo "I got $NUMBER_OF_NODES nodes"

#User input processors per node
echo -n "Processors per node?: "; read PROCESSOR_PER_NODE
if [ "$PROCESSOR_PER_NODE" == "" ]; then
	PROCESSOR_PER_NODE=1
fi
echo "I got $PROCESSOR_PER_NODE processors per node"

#User input queue
echo -n "Queue name?: "; read QUEUE
if [ "$QUEUE" == "" ]; then
	QUEUE="default"
fi
echo "$QUEUE is the name of the queue"

#User input virtual memory
echo -n "Virtual memory amount in gb?: "; read VIRTUAL_MEMORY
if [ "$VIRTUAL_MEMORY" == "" ]; then
	VIRTUAL_MEMORY=5
fi
echo "I got ${VIRTUAL_MEMORY}gb of virtual memory"

#User input memory
echo -n "Memory amount in gb?: "; read MEMORY
if [ "$MEMORY" == "" ]; then
	MEMORY=5
fi
echo "I got ${MEMORY}gb of memory"

<<DEBUG2
echo "echo "qsub -N $PIPELINE_NAME -l nodes=$NUMBER_OF_NODES:ppn=$PROCESSOR_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb -V -q $QUEUE"" | bash
echo "qsub -N $PIPELINE_NAME -l nodes=$NUMBER_OF_NODES:ppn=$PROCESSOR_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb -V -q $QUEUE"
DEBUG2

<<DEBUG1
echo "
# This is a comment
echo "Hello world!"
echo "This is not a comment!"
#This is another comment!
cd $CASA
mkdir $MAPOUT
mkdir $ASSEMOUT" | bash
DEBUG1


#<<MAIN_PROGRAM
echo "
# Set working directory to HOME
cd $CASA
# Creating output directories
mkdir $MAPOUT
mkdir $ASSEMOUT
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
	hisat2-build $GENOMAS/$GENOMA_FASTA_FILE $INDICE/$GENOMA_FASTA_FILE
fi

# Mapping reads to reference genome 
hisat2 -q -x $INDICE/$GENOMA_FASTA_FILE \
$HISAT2_OPTIONAL_PARAMETERS \
-1 $LECTURAS/$READS_F -2 $LECTURAS/$READS_R \
-S $MAPOUT/$SAM_FILE

# SAM <=> BAM
samtools view -bT $GENOMAS/$GENOMA_FASTA_FILE $MAPOUT/$SAM_FILE > $MAPOUT/${SAM_FILE%.*}.bam
# Sorting BAM
samtools sort $MAPOUT/${SAM_FILE%.*}.bam -o $MAPOUT/${SAM_FILE%.*}.sort.bam
# Indexing BAM
samtools index $MAPOUT/${SAM_FILE%.*}.sort.bam

# Assembling genome and estimating RNA abundances
# Creating specific output directory
mkdir $ASSEMOUT/${SAM_FILE%.*}
cufflinks -L $LABEL $CUFFLINKS_OPTIONAL_PARAMETERS  -o $ASSEMOUT/${SAM_FILE%.*}  $MAPOUT/${SAM_FILE%.*}.sort.bam
echo "$ASSEMOUT/${SAM_FILE%.*}/transcripts.gtf" >> $ASSEMOUT/GTFs.txt" | qsub -N $PIPELINE_NAME -l nodes=$NUMBER_OF_NODES:ppn=$PROCESSOR_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb -V -q $QUEUE
#MAIN_PROGRAM
#__EOF__//~~CAT
