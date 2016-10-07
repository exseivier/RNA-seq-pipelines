#!/usr/bin/env bash

while IFS='' read -r LINE || [[ -n $LINE ]];
do
	echo $LINE
	arrLINE=(${LINE//" "/ })
	key=${arrLINE[0]}
	value=${arrLINE[1]}
	case $key in 
		-pr)
		PRWD=$value
		;;
		-i)
		INDEX=$value
		;;
		-d)
		READS=$value
		;;
		-g)
		GENOMES=$value
		;;
		-t)
		TRANSCRIPTOMES=$value
		;;
		-r)
		REFGENOME=$value
		;;
		-c)
		COMMANDS=$value
		;;
		-s)
		SRAFILES=$value
		;;
		-tr)
		TRIMMING_COMMANDS=$value
		;;
		-a)
		PIPELINE_COMMANDS=$value
		;;
		*)
		echo "Bad option"
		;;
	esac
	shift
done < $1

echo $PRWD $INDEX $GENOMES $REFGENOME $TRANSCRIPTOMES $COMMANDS $READS $SRAFILES $TRIMMING_COMMANDS $PIPELINE_COMMANDS
echo "Redirecting..."
cd $PRWD
echo "You have been redirected to $PRWD"
echo ""
echo "Exporting variables"
export PRWD=$PRWD
export INDEX=$INDEX
export READS=$READS
export GENOMES=$GENOMES
export TRANSCRIPTOMES=$TRANSCRIPTOMES
export REFGENOME=$REFGENOME
export COMMANDS=$COMMANDS
export SRAFILES=$SRAFILES
export TRIMMING_COMMANDS=$TRIMMING_COMMANDS
export PIPELINE_COMMANDS=$PIPELINE_COMMANDS
echo ""
echo "Variables were exported! Have a nice work!..."
