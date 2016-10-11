#!/usr/bin/env bash

	# Parsing arguments
BATCH_PROCESS=""
SINGLE_PROCESS=""
COMMFILE=""
while [[ $# -gt 0 ]];
do
	key=$1
	case $key in
		-b|--batch)
		BATCH_PROCESS="TRUE"
		SINGLE_PROCESS="FALSE"
		;;
		-s|--single)
		SINGLE_PROCESS="TRUE"
		BATCH_PROCESS="FALSE"
		;;
		-c|--commfile)
		COMMFILE=$2
		shift
		;;
	esac
	shift
done

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

	# qsub settings
echo ""
echo "QSUB settings"
echo ""
echo -n "Number of nodes?: [Press Enter if default (1)]"; read NUMBER_OF_NODES
echo -n "Processors per node?: [Press ENETR if default (1)]"; read PROCESSORS_PER_NODE
echo -n "Queue?: [Press ENTER if (default)]"; read QUEUE
echo -n "Virtual memory?: [Press ENTER if default (5gb)]"; read VIRTUAL_MEMORY
echo -n "Memory?: [Press ENTER if default (5gb)]"; read MEMORY
if [ "$NUMBER_OF_NODES" == "" ]; then
	NUMBER_OF_NODES=1
fi
if [ "$PROCESSORS_PER_NODE" == "" ]; then
	PROCESSORS_PER_NODE=1
fi
if [ "$QUEUE" == "" ]; then
	QUEUE="default"
fi
if [ "$VIRTUAL_MEMORY" == "" ]; then
	VIRTUAL_MEMORY=5
fi
if [ "$MEMORY" == "" ]; then
	MEMORY=5
fi

echo $BATCH_PROCESS $SINGLE_PROCESS $COMMFILE
if [ "$BATCH_PROCESS" == "" ] && [ "$SINGLE_PROCESS" == "" ]; then
	echo "Usage:"
	echo "trimiando.sh [options]: [-b|-s] -c  <command file name>"
	echo "Options:"
	echo "	-b|--batch: Processes in batch mode"
	echo "	-s|--single: Processes in single mode"
	echo "	-c|--commfile: The command file name (this field is needed)"
	exit 1
fi

	# Running Trimmomatic in SINGLE mode
if [ "$BATCH_PROCESS" == "FALSE" ] && [ "$SINGLE_PROCESS" == "TRUE" ]; then
	echo -n "Name process?: "; read NAME
	if [ "$NAME" == "" ]; then
		NAME="Unknown"
	fi
	echo -n "Forward reads?: "; read READS_F
	echo -n "Reverse reads?: "; read READS_R
	if [ "$READS_F" == "" ] || [ "$READS_R" == "" ]; then
		echo "Reads files name are needed!"
		echo "Try again!"
		exit 1
	else
		echo "I got $READS_F and $READS_R"
	fi

	echo -n "Trimmomatic parameters?: "; read TRIMMOMATIC_OPTIONS
	if [ "$TRIMMOMATIC_OPTIONS" == "" ]; then
		echo "You have to specify the minimal parameters to run Trimmomatic!"
		echo "Try again!"
		exit 1
	else
		echo "I got $TRIMMOMATIC_OPTIONS"
	fi
	echo "
# Module Trimmomatic-0.32 loading
module load Trimmomatic/0.32
# Setting working directory
cd $PRWD/$READS
# Path to Trimmomatic
#TRIMMO=/data/software/Trimmomatic-0.32
# Running Trimmomatic
java -jar -Xmx1024m $TRIMMO/trimmomatic-0.32.jar \
PE \
$READS_F $READS_R \
trim${READS_F%.fastq.gz}_P.fastq.gz trim${READS_F%.fastq.gz}_U.fastq.gz \
trim${READS_R%.fastq.gz}_P.fastq.gz trim${READS_R%.fastq.gz}_U.fastq.gz \
ILLUMINACLIP:$ADAPTERS:$MISMATCH:$PALINDROME_CLIP:$SIMPLE_CLIP
$TRIMMOMATIC_OPTIONS" | qsub -N $NAME \
-l nodes=$NUMBER_OF_NODES:ppn=$PROCESSORS_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb \
-V -q $QUEUE
	
	# Running Trimmomatic in BATCH mode
elif [ "$BATCH_PROCESS" == "TRUE" ] && [ "$SINGLE_PROCESS" == "FALSE" ]; then
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
		echo "
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
$TRIMMOMATIC_OPTIONS" | qsub -N ${NAME%.sam} \
-l nodes=$NUMBER_OF_NODES:ppn=$PROCESSORS_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb \
-V -q $QUEUE
	done < "$COMMANDS/$COMMAND_FILE"
fi

