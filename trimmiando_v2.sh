#!/usr/bin/env bash
# Parsing arguments


	# User input data
echo ""
echo "User input data"
echo ""
echo -n "Adapters file path?: "; read ADAPTERS
echo -n "Trimmomatic program path?: "; read TRIMMO
if [ "$ADAPTERS" == "" ] || [ "$TRIMMO" == "" ]; then
	echo "You have to specify those paths"
	echo "Try again!"
	exit 1
fi
echo ""
echo "ILLUMINACLIP additional parameters"
echo ""
	# ILLUMINACLIP additional parameters
echo -n "Mismatch allowed?: [default (2)]"; read MISMATCH
echo -n "Palindrome clip thershold?: [default (15)]"; read PALINDROME_CLIP
echo -n "Simple clip threshold?: [default (10)]"; read SIMPLE_CLIP
if [ "$MISMATCH" == "" ]; then
	MISMATCH=2
fi
if [ "$PALINDROME_CLIP" == "" ]; then
	PALINDROME_CLIP=15
fi
if [ "$SIMPLE_CLIP" == "" ]; then
	SIMPLE_CLIP=10
fi
job=""
tjob=""
while IFS='' read -r LINE || [[ -n "$LINE" ]]; do
	arrLINE=(${LINE//"|"/ })
	NAME=${arrLINE[0]}
	READS_F=${arrLINE[1]}
	READS_R=${arrLINE[2]}
	TRIMMOMATIC_OPTIONS=$( echo ${arrLINE[3]} | sed 's/;/ /g')
	if [ "$READS_F" == "" ] || [ "$READS_R" == "" ]; then
		echo "You have to specify the filenames of the reads!"
		echo "Try again!"
		exit 1
	fi
	if [ "$TRIMMOMATIC_OPTIONS" == "" ]; then
		echo "Trimmomatic options are needed!"
		echo "Try again!"
		exit 1
	fi
	echo "I got $READS_F and $READS_R"
	echo "I got $TRIMMOMATIC_OPTIONS"
tjob=$(echo "
# Module Trimmomatic-0.32 loading
module load Trimmomatic/0.32
# Setting working directory
cd $PRWD
# Path to Trimmomatic
#TRIMMO=/data/software/Trimmomatic-0.32
# Running Trimmomatic
java -jar -Xmx1024m $TRIMMO/trimmomatic-0.32.jar \
PE \
$READS/$READS_F $READS/$READS_R \
$READS/trim${READS_F%.fastq.gz}_P.fastq.gz $READS/trim${READS_F%.fastq.gz}_U.fastq.gz \
$READS/trim${READS_R%.fastq.gz}_P.fastq.gz $READS/trim${READS_R%.fastq.gz}_U.fastq.gz \
ILLUMINACLIP:$ADAPTERS:$MISMATCH:$PALINDROME_CLIP:$SIMPLE_CLIP \
$TRIMMOMATIC_OPTIONS" | qsub -N ${NAME%.sam}_trimm \
-l nodes=1:ppn=8,vmem=10gb,mem=10gb \
-V -q default)
if [ "$job" == "" ]; then
	job=$tjob
else
	job=${job}:${tjob}
fi
done < "$COMMANDS/$COMMAND_FILE"

# FastQC analysis
cd $PRWD
while IFS='' read -r LINE || [[ -n "$LINE" ]];
do
	arrLINE=(${LINE//"|"/ })
	NAME=${arrLINE[0]}
	READS_F=${arrLINE[1]}
	READS_R=${arrLINE[2]}
	for i in {1..2};
	do
		READS_FoR=${arrLINE[$i]}
		echo "trim${READS_FoR%.fastq.gz}_P.fastq.gz"
echo "
# qsub feed begins
cd $PRWD/$READS
module load FastQC/0.11.2
READS_FoR=${arrLINE[$i]}
fastqc trim${READS_FoR%.fastq.gz}_P.fastq.gz" | qsub -N trim${NAME%.sam}_fqc.$i -l nodes=1:ppn=8,vmem=10gb,mem=10gb -W depend=afterok:$job -V -q default
	done
done < $COMMANDS/$COMMAND_FILE

