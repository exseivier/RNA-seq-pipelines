#!/usr/bin/env bash
while [[ $# -gt 0 ]];
do
	key=$1
	echo $key
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
		*)
		echo "WTF!"
		;;
	esac
	shift
done

cd $PRWD
[ $VCF_OUT ] || export VCF_OUT="vcf_out"
[ -d $VCF_OUT ] || mkdir $VCF_OUT
[ $MAPOUT ] || export MAPOUT="hisat2_out"
[ -d $MAPOUT ] || mkdir $MAPOUT
[ $LOCAL_LOG ] || export LOCAL_LOG="logs"
[ -d $LOCAL_LOG ] || mkdir $LOCAL_LOG
cd $PRWD/$LOCAL_LOG
while IFS='' read -r LINE || [[ -n "$LINE" ]];
do
	arrLINE=(${LINE//"|"/ })
	NAME=${arrLINE[0]}
	FILES=$(echo ${arrLINE[1]} | sed 's/;/ /g')
	VCF_OPTIONAL_PARAMETERS=$(echo ${arrLINE[2]} | sed 's/;/ /g')
	echo "
cd $PRWD/$MAPOUT
module load samtools/1.3.1
samtools mpileup -g -f $PRWD/$GENOMES/$REFGENOME \
$FILES \
$VCF_OPTIONAL_PARAMETERS \
> $PRWD/$VCF_OUT/${NAME}.bcf" |\
qsub -N $NAME -l nodes=$NODES:ppn=$PPN,vmem=${VMEM}gb,mem=${MEM}gb -V -q default

done < $PRWD/$COMMANDS/$VCF_COMMANDS
