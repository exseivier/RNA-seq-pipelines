#!/usr/bin/env bash

while IFS='' read -r LINE || [[ -n $LINE ]];
do
	echo $LINE
	arrLINE=(${LINE//" "/ })
	key=${arrLINE[0]}
	value=${arrLINE[1]}
	case $key in 
		-pr)
		PRWD=$value
		;;
		-i)
		INDEX=$value
		;;
		-d)
		READS=$value
		;;
		-g)
		GENOMES=$value
		;;
		-t)
		TRANSCRIPTOMES=$value
		;;
		-r)
		REFGENOME=$value
		;;
		-mo)
		MAPOUT=$value
		;;
		-ao)
		ASSEMOUT=$value
		;;
		-c)
		COMMANDS=$value
		;;
		-tc)
		COMMAND_FILE=$value
		;;
		-mac)
		MAPPING_ASSEMBLY_COMMANDS=$value
		;;
		-dec)
		DIFFEXP_COMMANDS=$value
		;;
		-s)
		SRAFILES=$value
		;;
		-tf)
		TRANSCRIPTOME_FILE=$value
		;;
		-log)
		LOCAL_LOG=$value
		;;
		-vc)
		VCF_COMMANDS=$value
		;;
		-vo)
		VCF_OUT=$value
		;;
		--cm-out)
		CUFFMERGE_OUT=$value
		;;
		--cd-out)
		CUFFDIFF_OUT=$value
		;;
		*)
		echo "Bad option"
		;;
	esac
	shift
done < $1

export PRWD=$PRWD; [ -d $PRWD ] || mkdir $PRWD
export INDEX=$INDEX; [ -d $PRWD/$INDEX ] || mkdir $PRWD/$INDEX
export READS=$READS; [ -d $PRWD/$READS ] || mkdir $PRWD/$READS
export GENOMES=$GENOMES; [ -d $PRWD/$GENOMES ] || mkdir $PRWD/$GENOMES
export TRANSCRIPTOMES=$TRANSCRIPTOMES
export REFGENOME=$REFGENOME
export COMMANDS=$COMMANDS; [ -d $PRWD/$COMMANDS ] || mkdir $PRWD/$COMMANDS
export SRAFILES=$SRAFILES; [ -f $PRWD/$READS/$SRAFILES ] || touch $PRWD/$READS/$SRAFILES
export COMMAND_FILE=$COMMAND_FILE; [ -f $PRWD/$COMMANDS/$COMMAND_FILE ] || touch $PRWD/$COMMANDS/$COMMAND_FILE
export MAPPING_ASSEMBLY_COMMANDS=$MAPPING_ASSEMBLY_COMMANDS; [ -f $PRWD/$COMMANDS/$MAPPING_ASSEMBLY_COMMANDS ] || touch $PRWD/$COMMANDS/$MAPPING_ASSEMBLY_COMMANDS
export DIFFEXP_COMMANDS=$DIFFEXP_COMMANDS; [ -f $PRWD/$COMMANDS/$DIFFEXP_COMMANDS ] || touch $PRWD/$COMMANDS/$DIFFEXP_COMMANDS
export TRANSCRIPTOME_FILE=$TRANSCRIPTOME_FILE; [ -f $PRWD/$TRANSCRIPTOMES/$TRANSCRIPTOME_FILE ] || touch $PRWD/$TRANSCRIPTOMES/$TRANSCRIPTOME_FILE
export LOCAL_LOG=$LOCAL_LOG; [ -d $PRWD/$LOCAL_LOG ] || mkdir $PRWD/$LOCAL_LOG
export VCF_COMMANDS=$VCF_COMMANDS; [ -f $PRWD/$COMMANDS/$VCF_COMMANDS ] || touch $PRWD/$COMMANDS/$VCF_COMMANDS
export VCF_OUT=$VCF_OUT; [ -d $PRWD/$VCF_OUT ] || mkdir $PRWD/$VCF_OUT
export MAPOUT=$MAPOUT; [ -d $PRWD/$MAPOUT ] || mkdir $PRWD/$MAPOUT
export ASSEMOUT=$ASSEMOUT; [ -d $PRWD/$ASSEMOUT ] || mkdir $PRWD/$ASSEMOUT
export CUFFMERGE_OUT=$CUFFMERGE_OUT; [ -d $PRWD/$CUFFMERGE_OUT ] || mkdir $PRWD/$CUFFMERGE_OUT
export CUFFDIFF_OUT=$CUFFDIFF_OUT; [ -d $PRWD/$CUFFDIFF_OUT ] || mkdir $PRWD/$CUFFDIFF_OUT

echo $PRWD $INDEX $GENOMES $REFGENOME $TRANSCRIPTOMES $COMMANDS $READS $SRAFILES $TRIMMING_COMMANDS $PIPELINE_COMMANDS
echo "Redirecting..."
cd $PRWD
echo "You have been redirected to $PRWD"
echo ""
echo "Exporting variables"

echo ""
echo "Variables were exported! Have a nice work!..."

#__EOF__//~~CAT
