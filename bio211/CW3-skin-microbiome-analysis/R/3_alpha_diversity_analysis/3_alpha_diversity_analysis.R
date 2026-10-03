
# Run from this script directory
################################################################################
# Preparation
################################################################################

# Load packages (Install them first if have not)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(OTUtable)
library(car)
library(viridis)
library(reshape2)
library(GUniFrac)

# Set seed for random processes
set.seed(1)

################################################################################
# 3.a. Rarefaction
################################################################################

# Generate long form data
OTU_long <- read_tsv("../../../data/OTU_count_adundance_long.tsv")
OTU_long <- OTU_long[, c(1:ncol(OTU_long)-1)]
write_tsv(OTU_long,"OTU_count_long.tsv")
OTU_long <- read_tsv("OTU_count_long.tsv")

# Generate wide form data
OTU_wide <- OTU_long %>%
  pivot_wider(
    names_from = Species,
    values_from = count,
    values_fill = 0
  )
colnames(OTU_wide) <- c("Sample", paste0("Species", 1:ncol(OTU_wide)))
write_tsv(OTU_wide,"../../../data/OTU_count_wide.tsv")
OTU_wide <- read_tsv("../../../data/OTU_count_wide.tsv") %>% as.data.frame()
rownames(OTU_wide) <- OTU_wide$Sample
OTU_wide$Sample <- NULL

# A function for calculate shannon for each depth in each sample
compute_rarefy_shannon <- function(count_sample) {
  # Depth and points in each sample
  max_depth <- sum(count_sample)
  depths <- seq(10, max_depth, length.out = 100)
  
  # Sample and calculate shannon for different depth
  shannon_values <- sapply(depths, function(d) {
    rarefied <- rrarefy(count_sample, sample = d)
    diversity(rarefied, index = "shannon")
  })
  
  # Export
  data.frame(depth = depths, shannon = shannon_values)
}

# Calculate shannon for each depth in each sample
rare_list <- lapply(rownames(OTU_wide), function(s) {
  df <- compute_rarefy_shannon(OTU_wide[s, ])
  df$Sample <- s
  df
})
rare_df <- bind_rows(rare_list)
write_tsv(rare_df, "shannon_depth_rarefaction.tsv")
rare_df <- read_tsv("shannon_depth_rarefaction.tsv")

# Choose rarefaction depth
(chosen_depth <- min(rowSums(OTU_wide)))

# Add `max_depth` col in data (for coloring)
sample_max_depth <- tapply(rare_df$depth, rare_df$Sample, max)
rare_df$max_depth <- sample_max_depth[ rare_df$Sample ]

# Plot & save
(plot_rare <- ggplot(rare_df, aes(x = depth, y = shannon, group = Sample)) +
  geom_line(aes(color = max_depth), size = 1.2, alpha = 0.9) +
  geom_vline(xintercept = chosen_depth, linetype = "dashed", size = 0.8, color = "black") +
  scale_color_viridis(option = "plasma", direction = 1) +
  theme_bw() +
  labs(
    title = "Rarefaction Curve Based on Shannon Diversity",
    x = "Sequencing Depth (Sampled Reads)",
    y = "Shannon Diversity Index",
    color = "Max Depth \nof Sample",
    subtitle = paste0("rarefaction threshold = ", chosen_depth)
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold"),
    plot.subtitle = element_text(face = "bold", hjust = 1)
  ))

ggsave("fig_3a_rarefaction_curve.png", plot_rare, width = 10, height = 10, dpi = 300)

# Generate rarefied tables
rarefied_table <- Rarefy(OTU_wide, depth = chosen_depth)$otu.tab.rff
rarefied_table_long <- as.data.frame(as.table(rarefied_table))
names(rarefied_table_long)[1:3] <- c("Sample", "Species", "Count")
rarefied_table_long <- rarefied_table_long %>%
  filter(! Count == 0)
write_tsv(rarefied_table_long,"rarefied_table_long.tsv")
rarefied_table_wide <- dcast(rarefied_table_long, Species ~ Sample, value.var = "Count", fill = 0)
write_tsv(rarefied_table_wide,"rarefied_table_wide.tsv") 

################################################################################
# 3.b. Rarefaction
################################################################################

# Calculate richness
richness <- read_tsv("rarefied_table_long.tsv") %>%
  select(Sample) %>%
  mutate(richness = 1) %>%
  group_by(Sample) %>% 
  mutate(richness= sum(richness)) %>%
  ungroup() %>%
  unique()
write_tsv(richness,"alpha_diversity_richness.tsv")

# Caculate shannon diversity、
# Transpose the data frame
alpha <- read_tsv("rarefied_table_wide.tsv")
Sample <- colnames(alpha[2:ncol(alpha)])
Species <- alpha[1:nrow(alpha), 1]
alpha[,1] <- NULL
colnames(alpha) <- NULL
alpha_t <- data.frame(t(alpha))
alpha_t[is.na(alpha_t)] <- 0
colnames(alpha_t) <- t(Species)
# Calculate the shannon diversity for global samples
Shannon <- diversity(alpha_t, "shannon")
# Merge the sample info with alpha diversity
data <- data.frame(cbind(Sample, Shannon))
# Write out
write_tsv(data, "alpha_diversity_shannon.tsv")

# Merge alpha diversity results
richness <- read_tsv("alpha_diversity_richness.tsv")
shannon <- read_tsv("alpha_diversity_shannon.tsv")
merge <- left_join(richness, shannon)
write_tsv(merge,"alpha_diversity_merge.tsv")

# Visualization

# Read the alpha diversity table
table <- read_tsv("alpha_diversity_merge.tsv")

# Reshape the data to long format and filter out zero
table <- table %>%
  pivot_longer(cols = 2:ncol(table), names_to = "alpha", values_to = "value") %>%
  filter(value > 0)

# Read metadata and merge with the alpha diversity data
meta <- read_tsv("../../../data/metadata.txt")
table <- left_join(table, meta, by = "Sample")

# Density plot
(density_plot <- ggplot(table, aes(x = value, fill = alpha)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ alpha, scales = "free") +
    scale_fill_manual(values = c("#71b7ed", "#f57c6e")) +
    labs(
     title = "Alpha Diversity Density Plots",
     x = "Value",
     y = "Density"
     ) +
    theme_minimal() +
    theme(
     plot.margin = margin(15, 15, 15, 15),
     plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
     axis.title = element_text(face = "bold", size = 12),
     axis.text = element_text(face = "bold"),
     axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
     legend.position = "none",
     strip.background = element_rect(fill = "grey95", color = NA),
     strip.text = element_text(face = "bold", size = 12),
     plot.subtitle = element_text(face = "bold", hjust = 1)
     ))
ggsave("fig_3b_alpha_diversity_density_plots.png", plot = density_plot, width = 10, height = 5, dpi = 300)

# QQ plot
(qq_plot <- ggplot(table, aes(sample = value)) +
   stat_qq() +
   stat_qq_line() +
   facet_wrap(~ alpha, scales = "free") +
   theme_minimal() +
   labs(
     title = "Alpha Diversity Q-Q Plots",
     x = "Theoretical Quantiles",
     y = "Sample Quantiles"
   ) +
   theme_minimal() +
   theme(
     plot.margin = margin(15, 15, 15, 15),
     plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
     axis.title = element_text(face = "bold", size = 12),
     axis.text = element_text(face = "bold"),
     axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
     legend.position = "none",
     strip.background = element_rect(fill = "grey95", color = NA),
     strip.text = element_text(face = "bold", size = 12),
     plot.subtitle = element_text(face = "bold", hjust = 1)
   ))
ggsave("fig_3b_alpha_diversity_qq_plots.png", plot = qq_plot, width = 10, height = 5, dpi = 300)

# Violin plot
(violin_plot <- ggplot(table, aes(Cutotype, value)) +
    geom_violin(aes(fill = Cutotype), trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_jitter(aes(color = City), width = 0.1, size = 2, alpha = 0.8) +
    theme_classic() +
    facet_wrap(~alpha, scales = "free") +
    scale_fill_manual(values = c(
      "Cutotype 1" = "#71b7ed", 
      "Cutotype 2" = "#f57c6e"
    )) +
    scale_color_manual(values = c(
      "Baoding" = "#84c3b7",
      "Dalian" = "#fae69e"
    )) +
    labs(
      title = "Alpha Diversity Violin Plots",
      x = "Cutotype",
      y = "Value"
      ) +
    theme(
      plot.margin = margin(15, 15, 15, 15),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, vjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(face = "bold", hjust = 1)
      ))
ggsave("fig_3b_alpha_diversity_violin_plots.png", violin_plot, width = 10, height = 8, dpi = 300)

################################################################################
# 3.c. ANOVA
################################################################################

# ANOVA for simultaneous effect analysis

# Data
table_richness <- subset(table, alpha == "richness")
table_shannon <- subset(table, alpha == "Shannon")

# For richness

# Assumptions
# Independence: Independence is highly related to experimental design, hence we assume the observations are independent within the single group, and data are also independent among different groups.
# Normality: need residues from the model, hence we verify this later
# Homogeneity of variances for categorical variables (Levene's Test):
leveneTest(value ~ City * Cutotype, data = table_richness, center = median)

# ANOVA
aov_richness <- aov(value ~ City * Cutotype * Age, data = table_richness)
(aov_richness_typeIII <- Anova(aov_richness, type = 3))

# Verify normality 
# Shapiro-Wilk Test
shapiro.test(residuals(aov_richness))
# QQ plot
df_resid_richness <- data.frame(residuals = residuals(aov_richness))
(qq_plot_richness <- ggplot(df_resid_richness, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line() +
    theme_minimal() +
    labs(
      title = "Richness Q-Q Plots",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal() +
    theme(
      plot.margin = margin(15, 15, 15, 15),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
      legend.position = "none",
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(face = "bold", hjust = 1)
    ))
ggsave("fig_3c_richness_qq_plots.png", plot = qq_plot_richness, width = 10, height = 10, dpi = 300)

# Interaction Plot
plot_data_richness <- table_richness %>%
  group_by(City, Cutotype) %>%
  summarise(
    mean_value = mean(value),
    se_value = sd(value) / sqrt(n()),
    .groups = "drop"
  )
(interaction_plot_richness <- ggplot(
  plot_data_richness, 
  aes(x = City, y = mean_value, group = Cutotype, color = Cutotype)
) +
  geom_point(size = 5) +
  geom_line(linewidth = 1.5) +
  geom_errorbar(
    aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
    width = 0.15,
    linewidth = 1
  ) +
  scale_color_manual(values = c("#71b7ed", "#f57c6e")) +
  labs(
    title = "Interaction Plot: City × Cutotype on Richness",
    x = "City",
    y = "Mean Richness ± SE",
    color = "Cutotype"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey95", linewidth = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(15, 15, 15, 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(face = "bold", hjust = 1)
  ))
ggsave("fig_3c_richness_interaction_plot.png", interaction_plot_richness, width = 10, height = 10, dpi = 300)

# For shannon

# Assumptions
# Independence: Independence is highly related to experimental design, hence we assume the observations are independent within the single group, and data are also independent among different groups.
# Normality: need residues from the model, hence we verify this later
# Homogeneity of variances for categorical variables (Levene's Test):
leveneTest(value ~ City * Cutotype, data = table_shannon, center = median)

# ANOVA
aov_shannon <- aov(value ~ City * Cutotype * Age, data = table_shannon)
(aov_shannon_typeIII <- Anova(aov_shannon, type = 3))

# Verify normality 
# Shapiro-Wilk Test
shapiro.test(residuals(aov_shannon))
# QQ plot
df_resid_shannon <- data.frame(residuals = residuals(aov_shannon))
(qq_plot_shannon <- ggplot(df_resid_shannon, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line() +
    theme_minimal() +
    labs(
      title = "Shannon Q-Q Plots",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal() +
    theme(
      plot.margin = margin(15, 15, 15, 15),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
      legend.position = "none",
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(face = "bold", hjust = 1)
    ))
ggsave("fig_3c_shannon_qq_plots.png", plot = qq_plot_shannon, width = 10, height = 10, dpi = 300)

# Interaction Plot
plot_data_shannon <- table_shannon %>%
  group_by(City, Cutotype) %>%
  summarise(
    mean_value = mean(value),
    se_value = sd(value) / sqrt(n()),
    .groups = "drop"
  )
(interaction_plot_shannon <- ggplot(
  plot_data_shannon, 
  aes(x = City, y = mean_value, group = Cutotype, color = Cutotype)
) +
    geom_point(size = 5) +
    geom_line(linewidth = 1.5) +
    geom_errorbar(
      aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      width = 0.15,
      linewidth = 1
    ) +
    scale_color_manual(values = c("#71b7ed", "#f57c6e")) +
    labs(
      title = "Interaction Plot: City × Cutotype on Shannon",
      x = "City",
      y = "Mean Shannon ± SE",
      color = "Cutotype"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_line(color = "grey95", linewidth = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(15, 15, 15, 15),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
      legend.position = "top",
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(face = "bold", hjust = 1)
    ))
ggsave("fig_3c_shannon_interaction_plot.png", interaction_plot_shannon, width = 10, height = 10, dpi = 300)

