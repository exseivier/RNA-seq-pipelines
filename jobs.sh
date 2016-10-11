#!/usr/bin/env bash

cd ~/scripts
j1=$(qsub job1.sh)
j2=$(qsub -W depend=afterok:$j1 job2.sh)
qsub -W depend=afterok:$j2 job3.sh

