#!/usr/bin/env bash

cd $PRWD/$READS
[ -e ../$COMMANDS/$COMMAND_FILE ] && rm ../$COMMANDS/$COMMAND_FILE
for file in *.sra;
do
	files=$(ls ${file%.sra}_*.* | tr "\n" "|")
	echo "${file%.sra}.sam|${files}" >> ../$COMMANDS/$COMMAND_FILE
done;
#__EOF__//~~CAT
