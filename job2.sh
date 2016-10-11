#!/usr/bin/env bash

#PBS -N job2
#PBS -l nodes=1:ppn=1,vmem=2gb,mem=2gb
#PBS -V
#PBS -q default
cd ~/
sleep 30
if [ -e job1.out ]; then
	touch job2.out
else
	echo "Error no archivo job1.out" >> error.txt
	exit 1
fi
rm job2.out

