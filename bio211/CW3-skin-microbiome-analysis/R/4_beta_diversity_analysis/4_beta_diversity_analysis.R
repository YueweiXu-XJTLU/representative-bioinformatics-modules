
# Run from this script directory
################################################################################
# Preparation
################################################################################

# Load packages (Install them first if have not)
library(readr)
library(dplyr)
library(reshape2)
library(vegan)
library(ggplot2)

# Set seed for random processes
set.seed(1)

################################################################################
# 4.pre. Distance Matrices
################################################################################

# Read the previously saved data
taxa_wide <- read_tsv("../../../data/OTU_count_wide.tsv") %>% as.data.frame()
rownames(taxa_wide) <- taxa_wide$Sample
taxa_wide$Sample <- NULL

# Calculate the distance matrix
# Bray-Curtis dissimilarity
dm_bray <- vegdist(taxa_wide, method="bray") %>%  as.matrix()
# Jaccard distance
dm_jaccard <- vegdist(taxa_wide, method="jaccard", binary = TRUE) %>%  as.matrix()

################################################################################
# 4.a. PCoA
################################################################################

# Read meta information
Meta <- read_tsv("../../../data/metadata.txt")

# PCoA - bray
# Analysis
PCoA_bray <- cmdscale(dm_bray, k = 2, eig = TRUE)
DF_bray <- data.frame(Sample = row.names(PCoA_bray$points), PCoA1 = PCoA_bray$points[,1], PCoA2 = PCoA_bray$points[,2], row.names = NULL)
Eigenvalues_bray <- eigenvals(PCoA_bray)
VarianceExplained_bray <- Eigenvalues_bray /sum(Eigenvalues_bray)
VarianceExplained1_bray <- 100 * signif(VarianceExplained_bray[1], 2) # unit: %
VarianceExplained2_bray <- 100 * signif(VarianceExplained_bray[2], 2) # unit: %
PCoA_bray <- merge(DF_bray, Meta, by = "Sample", all.x = TRUE)
# Visualize
PCoA_bray$City <- factor(PCoA_bray$City, levels = c("Baoding", "Dalian" ))
PCoA_bray$Cutotype <- factor(PCoA_bray$Cutotype, levels = c("Cutotype 1", "Cutotype 2"))
(pcoa_bray_plot <- ggplot(PCoA_bray, aes(x = PCoA1, y = PCoA2, color=Cutotype)) + 
    geom_point(size=4, alpha = 0.9, aes(shape=City)) + 
    stat_ellipse(
      geom = "polygon", 
      alpha = 0.1,
      linewidth = 1,
      aes(fill = Cutotype, color = Cutotype)
    ) + 
    xlab(xlab) + 
    ylab(ylab) + 
    labs(
      title = "PCoA of Bray–Curtis Dissimilarity",
      x = paste0("PCoA1 (Variance explained: ", VarianceExplained1_bray,"%)"),
      y = paste0("PCoA2 (Variance explained: ", VarianceExplained2_bray,"%)")
    ) +
    scale_fill_manual(values = c("#71b7ed", "#f57c6e")) +
    scale_color_manual(values = c("#71b7ed", "#f57c6e")) +
    theme_bw() + 
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold"),
      plot.subtitle = element_text(face = "bold", hjust = 1)
    ))
ggsave("fig_4a_bray_pcoa.png", width = 10, height = 10, dpi = 300)

# PCoA - jaccard
# Analysis
PCoA_jaccard <- cmdscale(dm_jaccard, k = 2, eig = TRUE)
DF_jaccard <- data.frame(Sample = row.names(PCoA_jaccard$points), PCoA1 = PCoA_jaccard$points[,1], PCoA2 = PCoA_jaccard$points[,2], row.names = NULL)
Eigenvalues_jaccard <- eigenvals(PCoA_jaccard)
VarianceExplained_jaccard <- Eigenvalues_jaccard /sum(Eigenvalues_jaccard)
VarianceExplained1_jaccard <- 100 * signif(VarianceExplained_jaccard[1], 2) # unit: %
VarianceExplained2_jaccard <- 100 * signif(VarianceExplained_jaccard[2], 2) # unit: %
PCoA_jaccard <- merge(DF_jaccard, Meta, by = "Sample", all.x = TRUE)
# Visualize
PCoA_jaccard$City <- factor(PCoA_jaccard$City, levels = c("Baoding", "Dalian" ))
PCoA_jaccard$Cutotype <- factor(PCoA_jaccard$Cutotype, levels = c("Cutotype 1", "Cutotype 2"))
(pcoa_jaccard_plot <- ggplot(PCoA_jaccard, aes(x = PCoA1, y = PCoA2, color=Cutotype)) + 
    geom_point(size=4,alpha = 0.9, aes(shape=City)) + 
    stat_ellipse(
      geom = "polygon", 
      alpha = 0.1,
      linewidth = 1,
      aes(fill = Cutotype, color = Cutotype)
      ) + 
    xlab(xlab) + 
    ylab(ylab) + 
    labs(
      title = "PCoA of Binary Jaccard Distance",
      x = paste0("PCoA1 (Variance explained: ", VarianceExplained1_jaccard,"%)"),
      y = paste0("PCoA2 (Variance explained: ", VarianceExplained2_jaccard,"%)")
    ) +
    scale_fill_manual(values = c("#71b7ed", "#f57c6e")) +    
    scale_color_manual(values = c("#71b7ed", "#f57c6e")) +
    theme_bw() + 
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold"),
      plot.subtitle = element_text(face = "bold", hjust = 1)
    ))
ggsave("fig_4a_jaccard_pcoa.png", width = 10, height = 10, dpi = 300)

################################################################################
# 4.b. PERMANOVA
################################################################################

# Permanova - bray
# use `by="margin"` to prevent the influence of order
(adonis_bray <- adonis2(dm_bray ~ City + Cutotype + Age, data = Meta, by = "margin"))
# Permanova - jaccard
# use `by="margin"` to prevent the influence of order
(adonis_jaccard <- adonis2(dm_jaccard ~ City + Cutotype + Age, data = Meta, by = "margin"))

