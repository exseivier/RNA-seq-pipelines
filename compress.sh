#PBS -N compress_fq.gz
#PBS -l nodes=1:ppn=1
#PBS -V
#PBS -q default

export CASA="/LUSTRE/usuario/montalvo"
export LECTURAS="$CASA/sra"

cd $LECTURAS/
gzip -9 SRR1598911_1.fastq
gzip -9 SRR1598911_2.fastq
