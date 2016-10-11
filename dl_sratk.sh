#PBS -N dl_sratk
#PBS -l nodes=1:ppn=1
#PBS -V
#PBS -q default

cd $PBS_O_WORKDIR

wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.7.0/sratoolkit.2.7.0-ubuntu64.tar.gz -P $HOME/bin/

