#!/usr/bin/env bash

cd $PRWD
[ -e $COMMANDS/$TRIMMING_COMMANDS ] && rm $COMMANDS/$TRIMMING_COMMANDS
echo ""
echo "Welcome to fdump.sh!"
echo "Transforming SRA files into fsatq"
echo -n "Fastq_dump optional parameters, [Press ENTER if default]?: "
read OPTIONAL_PARAMETERS;
echo -n "Number of nodes?: "; read NODES; if [ "$NODES" == "" ]; then NODES=1; fi
echo -n "Processors per node?: "; read PROCESSOR_PER_NODE; if [ "$PROCESSOR_PER_NODE" == "" ]; then PROCESSOR_PER_NODE=1; fi
echo -n "Virtual memory?: "; read VIRTUAL_MEMORY; if [ "$VIRTUAL_MEMORY" == "" ]; then VIRTUAL_MEMORY=5; fi
echo -n "Memory?: "; read MEMORY; if [ "$MEMORY" == "" ]; then MEMORY=5; fi
echo -n "Queue?: "; read QUEUE; if [ "$QUEUE" == "" ]; then QUEUE="default"; fi

while IFS='' read -r LINE || [[ -n $LINE ]];
do
	arrLINE=(${LINE//"|"/ })

<<DEBUG_1					#...[OK]
	#echo ${arrLINE[1]}
DEBUG_1

<<DEBUG_2					#...[OK]
#	echo "--outdir $PRWD/$READS --gzip --split-files $OPTIONAL_PARAMETERS ${arrLINE[0]}"
#	echo "${arrLINE[0]%.*}"
#	echo "${arrLINE[0]%.*}_x.fastq.gz"
#	touch $READS/${arrLINE[0]%.*}_1.fastq.gz
#	touch $READS/${arrLINE[0]%.*}_2.fastq.gz
#	files=$(ls $READS/${arrLINE[0]%.sra}_*.* | tr "\n" "|" | sed 's/sra\///g')
#	echo "${arrLINE[0]%.sra}.sam|${files}" >> $COMMANDS/$COMMAND_FILE
#	echo "qsub -N $NAME -l nodes=$NODES:ppn=$PROCESSOR_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb -V -q $QUEUE"
DEBUG_2

#<<MAIN_PRPGRAM
	echo "
cd $PRWD
fastq-dump --outdir $PRWD/$READS \
--gzip \
--split-files \
$OPTIONAL_PARAMETERS \
${arrLINE[0]}" | qsub -N ${arrLINE[0]%.sra} -l nodes=$NODES:ppn=$PROCESSOR_PER_NODE,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb -V -q $QUEUE

#files=$(ls $READS/${arrLINE[0]%.sra}_*.* | tr "\n" "|" | sed 's/sra\///g')
#echo "${arrLINE[0]%.sra}.sam|${files}" >> $COMMANDS/$TRIMMING_COMMANDS" 
#MAIN_PRPGRAM

done < $READS/$SRAFILES
