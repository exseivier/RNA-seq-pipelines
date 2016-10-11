#PBS -N dl_sra
#PBS -l nodes=1:ppn=1
#PBS -V
#PBS -q default

cd $PBS_O_WORKDIR

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR159/SRR1598911/SRR1598911.sra -P $HOME/sra/

