#!/usr/bin/env bash

	# User input reads
	# Trimmomatic parameters
BATCH_PROCESS=""
SINGLE_PROCESS=""
COMMFILE=""
	# qsub settings
do
	key=$1
	echo $1 $2
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
echo $BATCH_PROCESS $SINGLE_PROCESS $COMMFILE
if [ "$BATCH_PROCESS" == "" ] && [ "$SINGLE_PROCESS" == "" ] && [ "$COMMFILE" == "" ]; then
	echo "Usage:"
	echo "trimiando.sh [options]: [-b|-s] -c  <command file name>"
	echo "Options:"
	echo "	-b|--batch: Processes in batch mode"
	echo "	-s|--single: Processes in single mode"
	echo "	-c|--commfile: The command file name (this field is needed)"
	exit 1
fi
if [ "$BATCH_PROCESS" == "TRUE" ] && [ "$COMMFILE" == "" ]; then
	echo "You have to specify the command file:"
	echo "trimmiando.sh -b -c <command file name>"
	echo "Try again!"
fi

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
		# Running Trimmomatic in SINGLE PROCESS
	echo "
# Module Trimmomatic-0.32 loading
module load Trimmomatic/0.32
# Setting working directory
cd /LUSTRE/usuario/montalvo/sra/
# Path to Trimmomatic
TRIMMO=/data/software/Trimmomatic-0.32
# Running Trimmomatic
java -jar -Xmx1024m $TRIMMO/trimmomatic-0.32.jar \
PE \
$READS_F $READS_R \
trim${READS_F%.fastq.gz}_P.fastq.gz trim${READS_F%.fastq.gz}_U.fastq.gz \
trim${READS_R%.fastq.gz}_P.fastq.gz trim${READS_R%.fastq.gz}_U.fastq.gz \
$TRIMMOMATIC_OPTIONS" | qsub -N $NAME \
-l nodes=$NUMBER_OF_NODES:ppn=$PROCESSOR_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb \
-V -q $QUEUE
	
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
			# Running Trimmomatic in BATCH PROCESS
		echo "
# Module Trimmomatic-0.32 loading
module load Trimmomatic/0.32
# Setting working directory
cd /LUSTRE/usuario/montalvo/sra/
# Path to Trimmomatic
TRIMMO=/data/software/Trimmomatic-0.32
# Running Trimmomatic
java -jar -Xmx1024m $TRIMMO/trimmomatic-0.32.jar \
PE \
$READS_F $READS_R \
trim${READS_F%.fastq.gz}_P.fastq.gz trim${READS_F%.fastq.gz}_U.fastq.gz \
trim${READS_R%.fastq.gz}_P.fastq.gz trim${READS_R%.fastq.gz}_U.fastq.gz \
$TRIMMOMATIC_OPTIONS" | qsub -N $NAME \
-l nodes=$NUMBER_OF_NODES:ppn=$PROCESSOR_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb \
-V -q $QUEUE
	done < "$COMMFILE"
fi

