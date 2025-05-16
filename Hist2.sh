#Install Hisat2 and samtools through conda, then run the following.

#Indexing the genome
hisat2-build -p 6 -f genome.fasta genome_index

#Mapping of RNA-Seq from the sequencing company reads to the genome
hisat2 --phred33 -dta -x genome_index -1 forward_reads.fq -2 reverse_reads.fq -S output.sam

#Convert the sam file to bam file
samtools view -@ 28 -Sb -o output.sam output.bam

#Sort the bam file
samtools sort output.bam -o sorted_output.bam