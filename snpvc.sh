#!/usr/bin/env bash

while [ $# -gt 0 ];
do
	key=$1
	case $key in
		-nodes)
		NODES=$2
		shift
		;;
		-ppn)
		PPN=$2
		shift
		;;
		-vmem)
		VMEM=$2
		shift
		;;
		-mem)
		MEM=$2
		shift
		;;
		-q)
		QUEUE=$2
		shift
		;;
		-qop)
		OP=$2
		shift
		;;
	esac
	shift
done

while IFS='' read -r LINE || [[ -n "$LINE" ]];
do

	arrLINE=(${LINE//"|"/ })
	NAME=${arrLINE[0]}
echo "
cd $PRWD/$VCF_OUT
module load bcftools/1.2
bcftools call -O z -m -A -v -V ${NAME}.bcf \
$OPTIONAL_PARAMETERS \
-o ${NAME}.vcf.gz" | \
qsub -N $NAME -l nodes=$NODES:ppn=$PPN,vmem=${VMEM}gb,mem=${MEM}gb -V -q $QUEUE

done < $COMMANDS/$VCF_COMMANDS
