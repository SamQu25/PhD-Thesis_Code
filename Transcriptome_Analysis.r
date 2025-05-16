#Install and load required packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("genefilter")
library(genefilter)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")
library(Rsubread)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2", force = TRUE)
library(DESeq2)
install.packages("fastmap")
library(fastmap)
install.packages("dplyr", dependencies=TRUE, force = TRUE)
install.packages("ggplot2", dependencies=TRUE)
install.packages("gplots", dependencies=TRUE)
install.packages("RColorBrewer", dependencies=TRUE)
install.packages("zoom", dependencies=TRUE)
install.packages("shiny", dependencies=TRUE)
install.packages("plotly", dependencies=TRUE)
install.packages("readr")

library(dplyr)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(zoom)
library(shiny)
library(plotly)
library(readr)

#==================================================================================

#List bam files
bam.files <- list.files(path = "/Path/to/sorted_bams",
				pattern = "\\.bam$",
				full.names = TRUE,
				recursive = FALSE)

#==================================================================================

#Build Index
buildindex(basename = "Index", reference = "genome.fasta")

#==================================================================================

#Summary of the proportion of mapped reads (Takes a while - ca. 220min)
props <- propmapped(files = bam.files)
props

#==================================================================================

#Counting / Can be paralleled (e.g. snakemaker; Takes ca. 220min)
fc <- featureCounts(bam.files,
                    isPairedEnd=TRUE,
                    isGTFAnnotationFile=TRUE,
                    annot.ext="augustus_output.gff",
                    GTF.featureType="CDS",
                    GTF.attrType="gene_id"
                    )
names(fc)
fc$stat

#==================================================================================

#Number of genes
dim(fc$counts)

#==================================================================================

#Make counts file containing gene name, their counts and what sample they belong to
fc_counts <- fc$counts
fc_counts <- dplyr::as_tibble(fc_counts, rownames = "gene_id")
colnames(fc_counts) <- c("gene_id", "Control_Alignment_1_sorted", "Control_Alignment_2_sorted", "Control_Alignment_3_sorted", "Treatment_Alignment_1_sorted", "Treatment_Alignment_2_sorted", "Treatment_Alignment_3_sorted")
fc_counts <- as.data.frame(fc_counts)

#==================================================================================

#Creating metadata file
metadata <- data.frame(
    id = c("Control_Alignment_1_sorted", "Control_Alignment_2_sorted", "Control_Alignment_3_sorted", "Treatment_Alignment_1_sorted", "Treatment_Alignment_2_sorted", "Treatment_Alignment_3_sorted"),
    dex = c("control", "control", "control", "treated", "treated", "treated")
)
metadata$id <- factor(metadata$id, level = c("Control_Alignment_1_sorted", "Control_Alignment_2_sorted", "Control_Alignment_3_sorted", "Treatment_Alignment_1_sorted", "Treatment_Alignment_2_sorted", "Treatment_Alignment_3_sorted"))

#==================================================================================

#Construct DESeqDataSet object
dds <- DESeqDataSetFromMatrix(  countData = fc_counts,
                                colData = metadata,
                                design = ~dex, tidy = TRUE
)

#==================================================================================

#Run DDSeq function
dds <- DESeq(dds)

#==================================================================================

#Result table
res <- results(dds)

#==================================================================================

#Short summary list by p-values
res <- res[order(res$padj),]
head(res)

#==================================================================================

#Diagnostic MA plot
plotMA(res, ylim = c(-10,10))

#==================================================================================

#Independent filtering
#For weakly expressed genes we can't see differential expression cause low read counts suffer from high read noise
#Can be displayed by examining ratio of small p-values (<0.01) for genes binned by mean normalized count
#Create bins using normalize function
qs <- c(0, quantile(res$baseMean[res$baseMean > 0], 0:7/7))
#"Cut" the genes into the bins
bins <- cut(res$baseMean, qs)
#Rename the levels of the bins using the middle point
levels(bins) <- paste0("~", round(0.5*qs[-1] + 0.5*qs[-length(qs)]))
#Caluclate the ration of £p£ values less than 0.01 for each bin
ratios <- tapply(res$padj, bins, function(p) mean(p < 0.01, na.rm = TRUE))
#Plot these ratios
barplot(ratios, xlab = "Mean normalized counts", ylab = "Ratio of small $p$ values")

#==================================================================================

#Removing these genes gains benefits for multiple testing adjustment. Removing weakly expressed genes we can find more genes to be significant which we keep and improve the test.
#Independent filtering is automatically done by DESeq2. Can be controlled by the results function.
#Filter is independent because it is blind to the sample assignment ("Treatment" / "Control")

metadata(res)$alpha
metadata(res)$filterThreshold
plot(metadata(res)$filterNumRej,
    type = "b", ylab = "Number of rejections",
    xlab = "Quantiles of filter")

lines(metadata(res)$lo.fit, col = "red")
abline(v = metadata(res)$filterTheta)

#==================================================================================

#Euclidean distance between samples on non-rlog normalized data
sampleDist <- dist(t(assay(dds)))
sampleDist

#==================================================================================

#Visualize distance as heatmap
sampleDistMatrix <- as.matrix(sampleDist)
rownames(sampleDistMatrix) <- paste(dds$dex, sep = "-")
colnames(sampleDistMatrix) <- NULL
colours = colorRampPalette(rev(brewer.pal(9, "Blues"))) (255)
heatmap.2(sampleDistMatrix, trace = "none", col = colours)

#==================================================================================

#Remove genes with no expression
padj_Na <- res[!is.na(res$padj), ]
#Extract genes with 2 > log2FC > -2
sub_1 <- subset(padj_Na, abs(log2FoldChange) < 2)
#Extract genes with 2 <= log2FC <= -2
sub_2 <- subset(padj_Na, abs(log2FoldChange) >= 2 & padj > 0.01)
#Extract genes with log2FC >= 2 and padj <= 0.01 -> Significantly up-regulated genes
sub_up <- subset(padj_Na, log2FoldChange >= 2 & padj <= 0.01)
#Extract genes with log2FC <= -2 and padj <= 0.01 -> Significantly down-regulated genes
sub_down <- subset(padj_Na, log2FoldChange <= -2 & padj <= 0.01)

nrow(padj_Na)
nrow(sub_1)
nrow(sub_2)
nrow(sub_up)
nrow(sub_down)

#==================================================================================

#Convert DESeq2 data into data frames
sub_1_data <- data.frame(
  gene = rownames(sub_1),
  log2FoldChange = sub_1$log2FoldChange,
  padj = sub_1$padj
)
sub_1_data <- sub_1_data %>% arrange(desc(log2FoldChange))
sub_2_data <- data.frame(
  gene = rownames(sub_2),
  log2FoldChange = sub_2$log2FoldChange,
  padj = sub_2$padj
)
sub_2_data <- sub_2_data %>% arrange(desc(log2FoldChange))
sub_up_data <- data.frame(
  gene = rownames(sub_up),
  log2FoldChange = sub_up$log2FoldChange,
  padj = sub_up$padj
)
sub_up_data <- sub_up_data %>% arrange(desc(log2FoldChange))
sub_down_data <- data.frame(
  gene = rownames(sub_down),
  log2FoldChange = sub_down$log2FoldChange,
  padj = sub_down$padj
)
sub_down_data <- sub_down_data %>% arrange(log2FoldChange)

#==================================================================================

#Blast search of each gene was conducted using the NCBI BLAST+ command line tool

install.packages("tibble")
library(tibble)

#Read in BLAST hits
blast_hits <- read_tsv("BLAST_Hits.tsv", col_names = TRUE, na = "", progress = FALSE, show_col_types = FALSE)

blast_hits <- column_to_rownames(blast_hits, var = "Gene")

#Order by number
numeric_order <- as.numeric(gsub("g", "", rownames(blast_hits)))
blast_hits <- blast_hits[order(numeric_order), ]

head(blast_hits, 35)
dim(blast_hits)

#==================================================================================

res_with_blast <- as.data.frame(res_with_blast)
#Order by number
numeric_order <- as.numeric(gsub("g", "", rownames(res_with_blast)))
res_with_blast <- res_with_blast[order(numeric_order), ]

head(res_with_blast, 20)
dim(res_with_blast)

#==================================================================================

#Merge res_with_blast with blast_hits
blast_hits$Gene <- rownames(blast_hits)
res_with_blast$Gene <- rownames(res_with_blast)

merged_blast_hits <- merge(res_with_blast, blast_hits, by = "Gene", all.x = TRUE)

rownames(merged_blast_hits) <- merged_blast_hits$Gene

merged_blast_hits$Gene <- NULL

#Order by number
numeric_order <- as.numeric(gsub("g", "", rownames(merged_blast_hits)))
merged_blast_hits <- merged_blast_hits[order(numeric_order), ]

head(merged_blast_hits)
dim(merged_blast_hits)

#==================================================================================

# Export the data frame to a CSV file
write.csv(merged_blast_hits, file = "res_with_blast.csv", row.names = TRUE)

#Remove genes with no expression
merged_blast_hits_NoNA <- merged_blast_hits[!is.na(merged_blast_hits$padj), ]

#Extract genes with 2 > log2FC > -2
merged_blast_hits_2log2FC2 <- subset(merged_blast_hits_NoNA, abs(log2FoldChange) < 2)

#Extract genes with 2 <= log2FC <= -2
merged_blast_hits_non_signif <- subset(merged_blast_hits_NoNA, abs(log2FoldChange) >= 2 & padj > 0.01)

#Extract genes with log2FC >= 2 and padj <= 0.01 -> Significantly up-regulated genes
merged_blast_hits_up <- subset(merged_blast_hits_NoNA, log2FoldChange >= 2 & padj <= 0.01)

#Extract genes with log2FC <= -2 and padj <= 0.01 -> Significantly down-regulated genes
merged_blast_hits_down <- subset(merged_blast_hits_NoNA, log2FoldChange <= -2 & padj <= 0.01)

#==================================================================================

# Export the data frame to a CSV file
write.csv(merged_blast_hits_NoNA, file = "merged_blast_hits_NoNA.csv", row.names = TRUE)
write.csv(merged_blast_hits_2log2FC2, file = "merged_blast_hits_2log2FC2.csv", row.names = TRUE)
write.csv(merged_blast_hits_non_signif, file = "merged_blast_hits_non_signif.csv", row.names = TRUE)
write.csv(merged_blast_hits_up, file = "merged_blast_hits_up.csv", row.names = TRUE)
write.csv(merged_blast_hits_down, file = "merged_blast_hits_down.csv", row.names = TRUE)

head(merged_blast_hits)

merged_blast_hits_rownames <- merged_blast_hits
merged_blast_hits_rownames$Gene_ID <- rownames(merged_blast_hits_rownames)
merged_blast_hits_rownames <- merged_blast_hits_rownames[, c("Gene_ID", setdiff(names(merged_blast_hits_rownames), "Gene_ID"))]
head(merged_blast_hits_rownames)

#===================================================================================

#KO mapping conducted using the KAAS server
#https://www.genome.jp/kaas-bin/kaas_main?mode=est_b

kaas_mapping <- read.table("KAAS_Annotation.txt", header = FALSE, sep = "\t", fill = TRUE, col.names = c("Gene_ID", "KO"))
head(kaas_mapping, 10)

#===================================================================================

#Merge KAAS result with merged_blast_hits
merged_blast_hits_kaas <- merge(merged_blast_hits_rownames, kaas_mapping, by = "Gene_ID", all.x = TRUE)
merged_blast_hits_kaas <- merged_blast_hits_kaas %>%
    mutate(num = as.numeric(str_extract(Gene_ID, "\\d+"))) %>%
    arrange(num) %>%
    dplyr::select(-num)
saveRDS(merged_blast_hits_kaas, file = "merged_blast_hits_kaas.rds")
head(merged_blast_hits_kaas, 10)

#====================================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("topGO")
install.packages("readxl")

library(topGO)
library(readxl)

#====================================================================================

#GO_result.tsv from earlier InterproScan analysis

GO_results <- read.delim("GO_result.tsv", header = FALSE, sep = "\t")
colnames(GO_results) <- c("ID", "Gene_ID", "GO_term")
GO_results <- GO_results[, -1]
GO_results$GO_term <- gsub("\\(.*?\\)", "", GO_results$GO_term)
head(GO_results)
print(length(unique(GO_results$Gene_ID)))

GO_results_pivot <- GO_results %>%
    group_by(Gene_ID) %>%
    summarise(GO_term = paste(GO_term, collapse = ","), .groups = "drop")
head(GO_results_pivot)
dim(GO_results_pivot)

GO_results_pivot <- GO_results_pivot %>%
    tidyr::separate(GO_term, into = paste0("GO_term", 1:10), sep = ",", remove = TRUE)
GO_results_pivot[is.na(GO_results_pivot)] <- ""

GO_results_pivot <- GO_results_pivot %>%
    mutate(Gene_number = as.numeric(gsub("g", "", Gene_ID))) %>%
    arrange(Gene_number) %>%
    dplyr::select(-Gene_number)

#====================================================================================

#head(merged_blast_hits_kaas)
#dim(merged_blast_hits_kaas)

all_genes <- merged_blast_hits_kaas[, 1:7]
genes_up <- merged_blast_hits_up[, 1:7]
genes_up$Gene_ID <- rownames(genes_up)
row.names(genes_up) <- NULL
genes_up <- genes_up[, c(ncol(genes_up), 1:(ncol(genes_up) - 1))]
genes_up <- genes_up[, -ncol(genes_up)]
genes_down <- merged_blast_hits_down[, 1:7]
genes_down$Gene_ID <- rownames(genes_down)
row.names(genes_down) <- NULL
genes_down <- genes_down[, c(ncol(genes_down), 1:(ncol(genes_down) - 1))]
genes_down <- genes_down[, -ncol(genes_down)]

head(all_genes)
dim(all_genes)
all_genes_list <- all_genes[, 1, drop = FALSE]
head(all_genes_list)
dim(all_genes_list)

head(genes_up)
dim(genes_up)
genes_up_list <- genes_up[, 1, drop = FALSE]
head(genes_up_list)
dim(genes_up_list)

head(genes_down)
dim(genes_down)
genes_down_list <- genes_down[, 1, drop = FALSE]
head(genes_down_list)
dim(genes_down_list)

GO_results_pivot$GO_terms <- apply(GO_results_pivot[, 2:ncol(GO_results_pivot)], 1, function(x) paste(x[!is.na(x) & x != ""], collapse = ","))

geneID2GO <- function(GO_results_pivot) {
    result <- strsplit(GO_results_pivot$GO_terms, split = ",")
    names(result) <- GO_results_pivot$Gene_ID
    return(result)
}
geneID2GO_map <- geneID2GO(GO_results_pivot)
str(geneID2GO_map)

all_genes_list <- as.vector(all_genes_list$Gene_ID)
genes_up_list <- as.vector(genes_up_list$Gene_ID)

gene_list_up <- factor(as.integer(all_genes_list %in% genes_up_list))
names(gene_list_up) <- all_genes_list
print(levels(gene_list_up))

topGOdata_up <- new(
    "topGOdata",
    ontology = "BP",
    allGenes = gene_list_up,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO_map,
    nodeSize = 10
)

result <- runTest(topGOdata_up, algorithm = "classic", statistic = "fisher")
result_table <- GenTable(topGOdata_up, classicFisher = result, orderBy = "classicFisher", topNodes = 10)
print(result_table)
pdf("Images/topGO_upregulated.pdf", width = 8, height = 12)
showSigOfNodes(topGOdata_up, score(result), firstSigNodes = 5, useInfo = "all")
dev.off()

go_df <- as.data.frame(result_table)
go_df$classicFisher <- as.numeric(go_df$classicFisher)
go_df <- go_df[order(go_df$classicFisher), ]
go_df$Term <- factor(go_df$Term, levels = go_df$Term[order(go_df$classicFisher)])
p <- ggplot(go_df, aes(x = reorder(Term, classicFisher), y = -log10(classicFisher))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    xlab("GO Terms") +
    ylab("-log10(p-value)") +
    ggtitle("Top Upregulated GO Terms") +
    theme_classic() +
    theme(
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        plot.title = element_text(face = "bold", size = 15)
    ) 

print(p)
ggsave("Images/topGO_upregulated_barplot.png", plot = p, width = 15, height = 12)

#====================================================================================

genes_down_list <- as.vector(genes_down_list$Gene_ID)

gene_list_down <- factor(as.integer(all_genes_list %in% genes_down_list))
names(gene_list_down) <- all_genes_list
print(levels(gene_list_down))

topGOdata_down <- new(
    "topGOdata",
    ontology = "BP",
    allGenes = gene_list_down,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO_map,
    nodeSize = 10
)

result <- runTest(topGOdata_down, algorithm = "classic", statistic = "fisher")
result_table <- GenTable(topGOdata_down, classicFisher = result, orderBy = "classicFisher", topNodes = 10)
print(result_table)
pdf("Images/topGO_downregulated.pdf", width = 8, height = 12)
showSigOfNodes(topGOdata_down, score(result), firstSigNodes = 5, useInfo = "all")
dev.off()

go_df <- as.data.frame(result_table)
go_df$classicFisher <- as.numeric(go_df$classicFisher)
go_df <- go_df[order(go_df$classicFisher), ]
go_df$Term <- factor(go_df$Term, levels = go_df$Term[order(go_df$classicFisher)])
g <- ggplot(go_df, aes(x = reorder(Term, classicFisher), y = -log10(classicFisher))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    xlab("GO Terms") +
    ylab("-log10(p-value)") +
    ggtitle("Top Downregulated GO Terms") +
    theme_classic() +
    theme(
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        plot.title = element_text(face = "bold", size = 15)
    )

print(g)
ggsave("Images/topGO_downregulated_barplot.png", plot = g, width = 15, height = 12)

#====================================================================================