#PBS -l nodes=1:ppn=28
#PBS -l walltime=160:00:00
#PBS -S /bin/bash
#PBS -N Canu
#PBS -o simple_job_script.out
#PBS -m a
#PBS -q long

#Go to working directory

cd /beegfs/work/workspace/ws/tu_bciqu01-PathToWorkingDirectory

#BAM to fasta file

conda activate samtools
samtools bam2fq input.bam | seqtk seq -A > input.fa
conda deactivate

#Load the module

module load bio/canu/2.2

#Run the command

canu -p Prefix_Name -d Directory_Name genomeSize=23m -untrimmed correctedErrorRate=0.12 maxInputCoverage=800 'batOptions=-eg 0.10 -sb 0.01 -dg 2 -db 1 -dr 3' -pacbio-hifi input.fasta useGrid=false