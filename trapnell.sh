#PBS -N Trapnell_et_al 
#PBS -l nodes=1:ppn=8,mem=20gb,vmem=20gb
#PBS -V
#PBS -q default

# Exporting environmental variables
# Home in LUSTRE
export CASA="/LUSTRE/usuario/montalvo"
# Genomes
export GENOMAS="$CASA/genomes"
# Genome indexes
export INDICE="$CASA/index"
# Reads
export LECTURAS="$CASA/sra"
# Mapping output
export MAPOUT="$CASA/hisat_out"
# Assembling
export ASSEMOUT="$CASA/cl_out"
# Merging transcripts

# Loading HISAT2
# Loading samttols
# Loading cufflinks, cuffmerge and cuffdiff
module load hisat2/2.0.4
module load samtools/1.3.1
module load cufflinks/2.2.1

# Creating reference genome index. If index is not created, uncomment this
# hisat2-build $GENOMAS/Amaranthus_hypochondriacus.faa $INDICE/Amaranthus_hypochondriacus.faa

# Mapping reads to reference genome
hisat2 -q -x $INDICE/Amaranthus_hypochondriacus.faa \
-1 $LECTURAS/trimSRR1598911_1_P.fastq.gz -2 $LECTURAS/trimSRR1598911_2_P.fastq.gz \
-S $MAPOUT/SRR1598911.sam
hisat2 -q -x $INDICE/Amaranthus_hypochondriacus.faa \
-1 $LECTURAS/trimSRR1598913_1_P.fastq.gz -2 $LECTURAS/trimSRR1598913_2_P.fastq.gz \
-S $MAPOUT/SRR1598913.sam

# SAM <=> BAM
# Changing to mapped reads directory
cd $MAPOUT/
samtools view -bT $GENOMAS/Amaranthus_hypochondriacus.faa SRR1598911.sam > SRR1598911.bam
samtools view -bT $GENOMAS/Amaranthus_hypochondriacus.faa SRR1598913.sam > SRR1598913.bam
# Sorting BAM
samtools sort SRR1598911.bam -o SRR1598911.sort.bam
samtools sort SRR1598913.bam -o SRR1598913.sort.bam
# Indexing BAM
samtools index SRR1598911.sort.bam
samtools index SRR1598913.sort.bam


# Assembling genome and estimating RNA abundances
# Creating directory to SRR1598911 assembly
mkdir $ASSEMOUT/SRR1598911
cufflinks -L Floral -o $ASSEMOUT/SRR1598911  SRR1598911.sort.bam
# Creating directory to SRR1598913 assembly
mkdir $ASSEMOUT/SRR1598913
cufflinks -L Root -o $ASSEMOUT/SRR1598913  SRR1598913.sort.bam
cd $ASSEMOUT/
echo "$ASSEMOUT/SRR1598911/transcripts.gtf" > $ASSEMOUT/GTFs.txt # This instruction will erase any data previously stored in GTFs.txt
echo "$ASSEMOUT/SRR1598913/transcripts.gtf" >> $ASSEMOUT/GTFs.txt

# Merging assembled transcripts
# Creating merged transcript directory
mkdir $CASA/cm_out
# Changing to merged transcript directory
cd $CASA/cm_out
cuffmerge -o $CASA/cm_out $ASSEMOUT/GTFs.txt

# Differential expression analysis
# Creating cuffdiff output directory
mkdir $CASA/cd_out
# Channging to the cuffdiff output directory
cd $CASA/cd_out
cuffdiff -o $CASA/cd_out -L Floral,Root $CASA/cm_out/merged.gtf $MAPOUT/SRR1598911.sort.bam $MAPOUT/SRR1598913.sort.bam

