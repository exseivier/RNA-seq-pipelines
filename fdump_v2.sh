#!/usr/bin/env bash

cd $PRWD
echo ""
echo "Welcome to fdump.sh!"
echo "Transforming SRA files into fsatq"
echo -n "Fastq_dump optional parameters, [Press ENTER if default]?: "
read OPTIONAL_PARAMETERS;

<<COMMENT1
echo -n "Number of nodes?: "; read NODES; if [ "$NODES" == "" ]; then NODES=1; fi
echo -n "Processors per node?: "; read PROCESSOR_PER_NODE; if [ "$PROCESSOR_PER_NODE" == "" ]; then PROCESSOR_PER_NODE=1; fi
echo -n "Virtual memory?: "; read VIRTUAL_MEMORY; if [ "$VIRTUAL_MEMORY" == "" ]; then VIRTUAL_MEMORY=5; fi
echo -n "Memory?: "; read MEMORY; if [ "$MEMORY" == "" ]; then MEMORY=5; fi
echo -n "Queue?: "; read QUEUE; if [ "$QUEUE" == "" ]; then QUEUE="default"; fi
COMMENT1

while IFS='' read -r LINE || [[ -n $LINE ]];
do
	arrLINE=(${LINE//"|"/ })
	tjob=""
	job=""
#<<MAIN_PRPGRAM
tjob=$(echo "
cd $PRWD
fastq-dump --outdir $PRWD/$READS \
--gzip \
--split-files \
$OPTIONAL_PARAMETERS \
$READS/${arrLINE[0]}" | qsub -N ${arrLINE[0]%.sra}_fdump -l nodes=1:ppn=8,vmem=10gb,mem=10gb -V -q default)

if [ "$job" == "" ]; then
	job=$tjob
else
	job=${job}:${tjob}
fi
#MAIN_PRPGRAM
done < $READS/$SRAFILES
echo $job
cd $PRWD
while IFS='' read -r LINE || [[ -n "$LINE" ]];
do
arrLINE=(${LINE//"|"/ })
for i in {1..2};
do
echo "
cd $PRWD/$READS
module load FastQC/0.11.2
fastqc ${arrLINE[0]%.sra}_${i}.fastq.gz" | \
qsub -N ${arrLINE[0]%.sra}_${i}_fqc -l nodes=1:ppn=8,vmem=10gb,mem=10gb -W depend=afterok:$job -V -q default
done
done < $READS/$SRAFILES
