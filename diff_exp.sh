#!/usr/bin/env bash
cd $PRWD
echo ""
echo "QSUB settings"
echo ""
echo -n "Nodes?: "; read NODES
[ $NODES ] || NODES=1
echo -n "Processors per node?: "; read PPN
[ $PPN ] || PPN=1
echo -n "Virtual memory?: "; read VIRMEM
[ $VIRMEM ] || VIRMEM=5
echo -n "Memory?: "; read MEM
[ $MEM ] || MEM=5
echo -n "Queue?: "; read QUEUE
[ $QUEUE ] || QUEUE="default"
echo -n "Cuffmerge output?: "; read CUFFMERGE_OUT
[ $CUFFMERGE_OUT ] || CUFFMERGE_OUT="cm_out"
echo -n "Cuffdiff output?: "; read CUFFDIFF_OUT
[ $CUFFDIFF_OUT ] || CUFFDIFF_OUT="cd_out"

[ -d $PRWD/$CUFFMERGE_OUT ] || mkdir $PRWD/$CUFFMERGE_OUT
export CUFFMERGE_OUT=$CUFFMERGE_OUT
[ -d $PRWD/$CUFFDIFF_OUT ] || mkdir $PRWD/$CUFFDIFF_OUT
export CUFFDIFF_OUT=$CUFFDIFF_OUT

while IFS='' read -r LINE || [[ -n "$LINE" ]];
do
arrLINE=(${LINE//"|"/ })
NAME=${arrLINE[0]}
MAPPED_READS_1=${arrLINE[1]}
MAPPED_READS_2=${arrLINE[2]}
CUFFMERGE_OPTIONAL_PARAMETERS=$( echo ${arrLINE[3]} | sed 's/;/ /g')
CUFFDIF_OPTIONAL_PARAMETERS=$( echo ${arrLINE[4]} | sed 's/;/ /g')
LABELS=${arrLINE[5]}

#<<MAIN_PROGRAM
echo "
module load cufflinks/2.2.1
cd $PRWD

cuffmerge -o $CUFFMERGE_OUT \
$CUFFMERGE_OPTIONAL_PARAMETERS \
$ASSEMOUT/GTFs.txt

cuffdiff -o $CUFFDIFF_OUT \
-L $LABELS \
$CUFFDIF_OPTIONAL_PARAMETERS \
$PRWD/$CUFFMERGE_OUT/merged.gtf \
$MAPOUT/$MAPPED_READS_1 $MAPOUT/$MAPPED_READS_2" \
| qsub -N $NAME -l nodes=$NODES:ppn=$PPN,vmem=${VIRMEM}gb,mem=${MEM}gb -V -q $QUEUE
#MAIN_PROGRAM
echo "Data..."
echo $NAME $MAPPED_READS_1 $MAPPED_READS_2 $CUFFMERGE_OPTIONAL_PARAMETERS $CUFFDIF_OPTIONAL_PARAMETERS $LABELS
echo $CUFFMERGE_OPTIONAL_PARAMETERS $CUFFDIF_OPTIONAL_PARAMETERS
done < $COMMANDS/$DIFFEXP_COMMANDS
