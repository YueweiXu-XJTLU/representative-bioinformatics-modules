
# Run from this script directory
################################################################################
# Preparation
################################################################################

# Load packages (Install them first if have not)
library(tidyverse)
library(patchwork)
library(dplyr)
library(tidyr)
library(readr)

# Set seed for random processes
set.seed(1)

################################################################################
# 2.a. Spearman + FDR
################################################################################

# Read the table into R
taxa_data <- read_tsv("../../../data/OTU_adundance_wide.tsv") %>% as.data.frame()
rownames(taxa_data) <- taxa_data$Sample
taxa_data$Sample <- NULL

pollutant_data <- read_tsv("../../../data/pollutant.txt") %>% as.data.frame()
rownames(pollutant_data) <- pollutant_data$Sample
pollutant_data$Sample <- NULL

# Multiple Spearman’s rank-order correlation analysis (use a nested loop)
species_names <- colnames(taxa_data)
pah_names <- colnames(pollutant_data)
cor_list <- list()

for (pah in pah_names) {
  for (sp in species_names) {
    
    x <- as.numeric(pollutant_data[[pah]])
    y <- as.numeric(taxa_data[[sp]])
    
    ct <- cor.test(x, y, method = "spearman")
    
    cor_list[[length(cor_list) + 1]] <- data.frame(
      Species = sp,
      PAH = pah,
      rho = unname(ct$estimate),
      S_score = unname(ct$statistic),
      p_value = ct$p.value,
      n = length(x)
    )
  }
}

cor_results <- bind_rows(cor_list)

# FDR correction (Benjamini-Hochberg method)
cor_results <- cor_results %>%
  mutate(adj_q_value = p.adjust(p_value, method = "BH"))

# Write all results
write_tsv(cor_results, "spearman_all.tsv")

# find most significant species for each PAHs
cor_significant <- cor_results %>%
  group_by(PAH) %>%
  slice_min(adj_q_value, n = 1, with_ties = TRUE) %>%
  ungroup() %>%
  arrange(adj_q_value)
write_tsv(cor_significant, "spearman_significant.tsv")

################################################################################
# 2.b. Spearman + FDR
################################################################################

# Prepare data
cor_plot <- read_tsv("spearman_significant.tsv") %>% 
  arrange(adj_q_value) %>%
  head(5)

# A function for subplots
plot_scatter <- function(species, pah, rho_val, adj_q_val) {
  # Extract relative abundance & PAH
  df <- data.frame(
    RA = as.numeric(taxa_data[[species]]),
    PAH = pollutant_data[[pah]]
  )
  
  # Annotate rho and adjusted q-value
  stats_text <- paste0(
    "Spearman’s ρ = ", round(rho_val, 3),
    "\nAdjusted p-value = ", format(adj_q_val, digits = 3, scientific = TRUE)
  )
  
  # Plot (subplot)
  ggplot(df, aes(x = PAH, y = RA)) +
    geom_point(color = "#71b7ed", alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", color = "#f57c6e", fill = "#f2656f") +
    scale_y_log10() +
    labs(
      title = species,
      subtitle = stats_text,
      x = pah,
      y = "Relative Abundance (log10)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold"),
      legend.position = "none",
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(face = "bold", hjust = 1)
    )
}

# Plot subplots
plot_list <- list()
for (i in 1:nrow(cor_plot)) {
  row <- cor_plot[i, ]
  p <- plot_scatter(
    species = row$Species,
    pah = row$PAH,
    rho_val = row$rho,
    adj_q_val = row$adj_q_value
  )
  plot_list[[i]] <- p
}

# Merge and save subplots
(final_plot <- wrap_plots(plotlist = plot_list, ncol = 2))
ggsave("fig_2b_cor_scatter.png", final_plot, width = 10, height = 15, dpi = 300)
