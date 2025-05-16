#Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rtracklayer")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("topGO")
BiocManager::install("Rgraphviz")
install.packages("readxl")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("writexl")
install.packages("tidyverse")
install.packages("svglite")

library(topGO)
library(readxl)
library(writexl)
library(ggplot2)
library(dplyr)
library(rtracklayer)
library(tidyverse)
library(svglite)
library(Rgraphviz)

#==================================================================================

# Replace "path_to_file.gff3" with the actual path to your GFF3 file
gff3_file <- "augustus_output.gff3"

# Import the GFF3 file
gff3_data <- import(gff3_file, format = "gff3")

# Convert to a data frame
gff3_df <- as.data.frame(gff3_data)

gff3_df$phase <- NULL
gff3_df$strand <- NULL
gff3_df$type <- NULL
gff3_df$source <- NULL
colnames(gff3_df)[colnames(gff3_df) == "seqnames"] <- "gene"

#==================================================================================

df_up <- read.csv("merged_blast_hits_up.csv", header = TRUE, stringsAsFactors = FALSE)
df_down <- read.csv("merged_blast_hits_down.csv", header = TRUE, stringsAsFactors = FALSE)

colnames(df_up)[colnames(df_up) == "X"] <- "gene"
colnames(df_down)[colnames(df_down) == "X"] <- "gene"

#==================================================================================

gff3_df <- gff3_df %>%
    mutate(Regulation = ifelse(gene %in% df_up$gene, "Upregulated", 
                               ifelse(gene %in% df_down$gene, "Downregulated", "Not regulated")))

df_up_selected <- df_up %>%
    select(gene, Description)

df_down_selected <- df_down %>%
    select(gene, Description) %>%
    rename(Description_down = Description)

gff3_df <- gff3_df %>%
    left_join(df_up_selected, by = "gene") %>%
    left_join(df_down_selected, by = "gene")

gff3_df$Description <- paste(
    coalesce(gff3_df$Description, ""),
    coalesce(gff3_df$Description_down, ""),
    sep = " "
)

gff3_df$Description_down <- NULL
gff3_df <- gff3_df[!is.na(gff3_df$Description), ]

gff3_df_up <- gff3_df[gff3_df$Regulation == "Upregulated", ]
gff3_df_down <- gff3_df[gff3_df$Regulation == "Downregulated", ]

write_xlsx(gff3_df_up, "gff3_df_up.xlsx")
write_xlsx(gff3_df_down, "gff3_df_down.xlsx")

#==================================================================================

regulation_counts <- gff3_df %>%
    count(Regulation)

regulation_counts$Regulation <- factor(
  regulation_counts$Regulation,
  levels = c("Not regulated", "Downregulated", "Upregulated"),
  labels = c("Non-DE Genes", "Downregulated Genes", "Upregulated Genes")
)

colors <- c("#D9D9D9", "#E4ECC8", "#A9B291")

plot <- ggplot(regulation_counts, aes(x = Regulation, y = n, fill = Regulation)) +
    geom_col() +
    labs(
         y = "Number of Genes"
         ) +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 55, face = "bold", color = "black"),
        axis.text.y = element_text(size = 40, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 60, face = "bold", color = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)
    ) +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = n), y = 8, color = "black", size = 20, fontface = "bold")

print(plot)

ggsave("Images/regulation_counts.png", plot, width = 30, height = 30, dpi = 600, bg = "transparent")
ggsave("Images/regulation_counts.pdf", plot, width = 10, height = 8, dpi = 300)
ggsave("Images/regulation_counts.svg", plot, width = 10, height = 8, dpi = 300)

#==================================================================================