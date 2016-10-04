#!/usr/bin/env bash

echo -n "Pattern?: "; read PATTERN
echo -n "Home?: "; read HOME

for file in $PATTERN;
do
	echo "cd $HOME; module load FastQC/0.11.2; fastqc $file" \
	| qsub -N ${file%.*.*} -l nodes=1:ppn=8,vmem=5gb,mem=5gb -V -q default;
done
