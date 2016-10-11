#!/usr/bin/env bash

#PBS -N job1
#PBS -l nodes=1:ppn=1,vmem=2gb,mem=2gb
#PBS -V
#PBS -q default
cd ~/
sleep 30
touch job1.out

