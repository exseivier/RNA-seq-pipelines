#!/usr/bin/env bash

echo -n "Pattern?: "; read PATTERN
echo -n "Nodes?: "; read NODES; if [ "$NODES" == "" ]; then NODES=1; fi
echo -n "Processor per node?: "; read PPN; if [ "$PPN" == "" ]; then PPN=1; fi
echo -n "Virtual memory?: "; read VIRTUAL_MEMORY; if [ "$VIRTUAL_MEMORY" == "" ]; then VIRTUAL_MEMORY=5; fi
echo -n "Memory?: "; read MEMORY; if [ "$MEMORY" == "" ]; then MEMORY=5; fi
echo -n "Queue?: "; read QUEUE; if [ "$QUEUE" == "" ]; then QUEUE="default"; fi

for file in $PATTERN;
do
	echo "cd $PRWD/$READS; module load FastQC/0.11.2; fastqc $file" \
	| qsub -N ${file%.*.*} -l nodes=$NODES:ppn=$PPN,vmem=${VIRTUAL_MEMORY}gb,mem=${MEMORY}gb -V -q $QUEUE
done
