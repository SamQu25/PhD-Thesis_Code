install.packages("Biostrings")
install.packages("rtracklayer")
install.packages("seqinr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("dplyr")

library(Biostrings)
library(rtracklayer)
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)

#Total bases in genome
fasta <- readDNAStringSet("genome.fasta")
gff <- import("augustus_output")
total_bases <- sum(width(fasta))
cat("Total bases in genome:", total_bases, "\n")
#==================================================================================

#Calculate the percentage of CDS regions
cds_region <- gff[gff$type == "CDS"]
cds_length <- sum(width(cds_region))
percentage_cds <- (cds_length / sum(width(fasta))) * 100
cat("Percentage CDS (%):", percentage_cds, "\n")
#==================================================================================

#Calcukate average gene size
gene_regions <- gff[gff$type == "gene"]
average_gene_size <- mean(width(gene_regions)) / 1e3
cat("Average gene size (kb):", average_gene_size, "\n")
#==================================================================================

#Calculate gene density
gene_count <- length(gene_regions)
total_length_kb <- sum(width(fasta)) / 1e3
gene_density <- gene_count / total_length_kb
cat("Average gene density (gene/kb):", gene_density, "\n")
#==================================================================================

#Calculate the number of genes
protein_coding_genes <- gene_regions[grepl("gene", gene_regions$type)]
cat("Number of protein-coding genes:", length(protein_coding_genes), "\n")
#==================================================================================

#Calculate the number of exons and introns
exon_regions <- gff[gff$type == "CDS"]
cat("Number of exons:", length(exon_regions), "\n")

intron_regions <- gff[gff$type == "intron"]
cat("Number of introns:", length(intron_regions), "\n")
#==================================================================================

#Calculate the intro/exon ratio
intron_exon_rate <- length(intron_regions) / length(exon_regions)
cat("Introns/Exons ratio:", intron_exon_rate, "\n")
#==================================================================================

#Calculate the number of introns per gene
introns_per_gene <- length(intron_regions) / length(protein_coding_genes)
cat("Introns per gene:", introns_per_gene, "\n")
#==================================================================================


#Calculate average intron length
exons_by_gene <- split(exon_regions, sapply(exon_regions$type, function(attr) {
  sub(".*gene_id \"([^\"]+)\".*", "\\1", attr)
}))
introns_by_gene <- lapply(exons_by_gene, gaps)
intron_lengths <- unlist(lapply(introns_by_gene, width))
average_intron_length <- mean(intron_lengths) / 1e3
cat("Average intron length (kb):", average_intron_length, "\n")
#==================================================================================

#Calculate the percentage of bases in intergene regions
gene_ranges <- gff[gff$type == "gene", ]
gene_lengths <- width(gene_ranges)
total_gene_bases <- sum(gene_lengths)
cat("Total number of bases in genes:", total_gene_bases, "\n")
intergenic_distances <- total_bases - total_gene_bases
cat("Total number of bases in intergenic regions:", intergenic_distances, "\n")
intergenic_distances_percent <- (intergenic_distances / total_bases) * 100
cat("Percentage of intergenic regions (%):", intergenic_distances_percent, "\n")
#==================================================================================

#Calculate N50, number of contigs, and GC content
contig_lengths <- width(fasta)

#Calculate N50 (in kb)
n50_value <- function(lengths) {
    sorted_lengths <- sort(lengths, decreasing = TRUE)
    cumsum_lengths <- cumsum(sorted_lengths)
    half_length <- sum(sorted_lengths) / 2
    n50 <- sorted_lengths[which(cumsum_lengths >= half_length)[1]]
    return(n50 / 1e3)
}

N50contig <- n50_value(contig_lengths)

#Number of Contigs
num_contigs <- length(contig_lengths)

#GC-Content
gc_content <- function(dna_sequences)  {
    gc_counts <- sum(alphabetFrequency(dna_sequences, baseOnly = TRUE)[, c("G", "C")])
    total_bases <- sum(width(dna_sequences))
    return((gc_counts / total_bases) * 100)
}

GC_content <- gc_content(fasta)

cat("N50 (kb):\t\t", N50contig, "\n")
cat("Number of contigs:\t", num_contigs, "\n")
cat("GC Content (%):\t\t", GC_content, "\n")
#==================================================================================