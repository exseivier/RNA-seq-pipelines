#!/usr/bin/env bash

#PBS -N job3
#PBS -l nodes=1:ppn=1,vmem=2gb,mem=2gb
#PBS -V
#PBS -q default
cd ~/
sleep 30
if [ -e job2.out ]; then
	touch final_job.txt
else
	echo "Error no archivo job2.out" >> error.txt
	exit 1
fi

