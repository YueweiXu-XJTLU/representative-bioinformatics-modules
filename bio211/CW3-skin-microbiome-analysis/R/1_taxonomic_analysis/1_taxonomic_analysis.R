
# Run from this script directory
################################################################################
# Preparation
################################################################################

# Load packages (Install them first if have not)
library(tidyverse) 
library(plyr)
library(wilkoxmisc)
library(maaslin3)
library(patchwork)

# Set seed for random processes
set.seed(1)

################################################################################
# 1.a. Relative Abundance & Stacked Bar Plot
################################################################################

# Read the table into R
OTU <- read_tsv("../../../data/species_table.txt")
names(OTU)[1] <- c("Sample")

# Calculate the count of each taxon in each sample and only remove taxa not identified
OTU <- OTU %>%
  gather(key = "Species", value = "Count", -Sample) %>%
  filter(Count > 0)

# Tabulate read counts by species
OTU <- ddply(OTU, .(Sample, Species), summarise, count = sum(Count))

# Convert count to relative abundance and add column (require dplyr)
# RelativeAbundance in unit %
OTU <- ddply(OTU, .(Sample), mutate, RelativeAbundance = (count * 100) / sum(count))

# Save the processed table into a tsv file
write_tsv(OTU,"../../../data/OTU_count_adundance_long.tsv")

# Prepare the table for plotting

# Read the previously saved tsv file
OTU <- read_tsv("../../../data/OTU_count_adundance_long.tsv")

# Collapse taxa table to only 10 top phyla, genus, family, etc.
if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}
OTUTable <- collapse_taxon_table(OTU, n = 10, Rank = "Species") 

# Read the metadata table
Meta <- read_tsv("../../../data/metadata.txt")

# Merge relative abundance table and metatable together
OTUTable <- merge(OTUTable, Meta, by = "Sample", all.x = TRUE)

# Save the table in tsv file
write_tsv(OTUTable, "Top10Species.tsv")

# Plotting
OTUTable <- read_tsv("Top10Species.tsv")

# Set factor for two variables
OTUTable$Species <- factor(OTUTable$Species, levels = unique(OTUTable$Species))
OTUTable$Sample <- factor(OTUTable$Sample,levels = unique(OTUTable$Sample))

# Add sample count to Cutotype label (for polt)
OTUTable <- OTUTable %>%
  group_by(Cutotype) %>%
  mutate(Cutotype_Label = paste0(Cutotype, " (", n_distinct(Sample), " Samples)")) %>%
  ungroup()

labels <- unique(OTUTable$Cutotype_Label)
OTUTable$Cutotype_Label <- factor(OTUTable$Cutotype_Label, levels = sort(labels))

# Plot the relative abundance
(plot_stacked_bar_chart <- ggplot(OTUTable, aes(x = Sample, y = RelativeAbundance, fill = Species)) +
    geom_bar(stat="identity") +
    scale_fill_brewer(palette = "Paired") +
    facet_wrap(~Cutotype_Label, scales = "free_x") +
    scale_y_continuous(expand = c(0, 0)) + 
    ylab("Relative abundance (%)") +
    xlab("Sample") + 
    labs(title = "Microbial Community Composition by Cutotype") +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold"),
      plot.margin = margin(15, 15, 15, 15),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      legend.text = element_text(face = "bold.italic", size = 10), 
      legend.title = element_text(face = "bold"),
      legend.key.size = unit(0.5, "cm"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = 12)
    ))

# Save the plot in pdf file
ggsave("fig_1a_composition_stacked_bar.png", plot_stacked_bar_chart, width = 10, height = 6, dpi = 300)

################################################################################
# 1.b. Multivariable Association Analysis
################################################################################

# Data Preparation

# Input taxa table
OTU_abundance_long <- read_tsv("../../../data/OTU_count_adundance_long.tsv")
OTU_abundance_wide <- OTU_abundance_long %>%
  select(Sample, Species, RelativeAbundance) %>%
  pivot_wider(
    names_from = Species,
    values_from = RelativeAbundance,
    values_fill = 0
  )
input_data <- OTU_abundance_wide
input_data <- as.data.frame(input_data)
names(input_data)[1] <- c("Sample")
write_tsv(input_data,"../../../data/OTU_adundance_wide.tsv")
input_data <- read_tsv("../../../data/OTU_adundance_wide.tsv") %>% as.data.frame()
rownames(input_data) <- input_data$Sample
input_data$Sample <- NULL

# Input meta data
input_metadata <- read_tsv("../../../data/metadata.txt")
input_metadata <- as.data.frame(input_metadata)
rownames(input_metadata) <- input_metadata$Sample
input_metadata$Sample <- NULL
input_metadata <- input_metadata %>%
  mutate(
    City     = factor(City, levels = unique(input_metadata$City)),
    Cutotype = factor(Cutotype, levels = unique(input_metadata$Cutotype)),
    Age      = as.numeric(Age)
  )

# Output directory
output_dir <- "1b_maaslin3_output"

# Conduct MaAsLin3 analysis
fit_out <- maaslin3 (
  input_data = input_data,
  input_metadata = input_metadata,
  output = output_dir,
  formula = '~ City + Cutotype + Age',
  normalization = 'TSS',
  transform = 'LOG',
  augment = TRUE,
  standardize = FALSE, # remain interpretability for Age in unit of year
  max_significance = 0.01,
  median_comparison_abundance = TRUE,
  median_comparison_prevalence = TRUE,
  max_pngs = 250,
  cores = 1
  )

# Interpretation and visualization
res <- read_tsv(paste0(output_dir, "/all_results.tsv"))

# A function to search for target species
get_top_species <- function(res_table, variable_name, direction = "both", num = 2) {
  subset_df <- res_table %>% 
    filter(metadata == variable_name) %>%
    arrange(qval_individual)
  
  if(direction == "both") {
    return(subset_df %>% head(num))
  } else if (direction == "positive") {
    return(subset_df %>% filter(coef > 0) %>% head(num))
  } else if (direction == "negative") {
    return(subset_df %>% filter(coef < 0) %>% head(num))
  }
}

# City (Categorical)
top_city <- get_top_species(res, "City", "both", 2)
(top_city_species <- top_city$feature)

# Cutotype (Categorical)
top_cutotype <- get_top_species(res, "Cutotype", "both", 2)
(top_cutotype_species <- top_cutotype$feature)

# C. Age (Continuous)
top_age_pos <- get_top_species(res, "Age", "positive", 1)
(top_age_pos_species <- top_age_pos$feature)
top_age_neg <- get_top_species(res, "Age", "negative", 1)
(top_age_neg_species <- top_age_neg$feature)

# Visualization

# Box plot for categorical variables

# A function for subplots
plot_boxplot <- function(species_name, meta_col, res_row) {
  # Extract relative abundance & metadata
  plot_data <- data.frame(
    Abundance = input_data[[species_name]], 
    Group = input_metadata[[meta_col]]
  )
  
  # Annotate coefficient and q-value
  stats_text <- paste0("Coefficient = ", round(res_row$coef, 3), 
                       "\nq-value = ", format(res_row$qval_individual, scientific = TRUE, digits = 3))
  
  # Plot (subplot)
  ggplot(plot_data, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5, color = "grey50") + 
    scale_y_log10() + # for better appearance, automatically remove 0 relative abundance
    scale_fill_manual(values = c("#71b7ed", "#f57c6e")) +
    labs(
      title = species_name, 
      subtitle = stats_text, 
      x = meta_col,
      y = "Relative Abundance (Log10)", 
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
(p1 <- plot_boxplot(top_city_species[1], "City", top_city[1,]))
(p2 <- plot_boxplot(top_city_species[2], "City", top_city[2,]))
(p3 <- plot_boxplot(top_cutotype_species[1], "Cutotype", top_cutotype[1,]))
(p4 <- plot_boxplot(top_cutotype_species[2], "Cutotype", top_cutotype[2,]))

# Merge and save subplots
(p_categorical_city <- (p1 + p2))
(p_categorical_cutotype <- (p3 + p4))
ggsave("fig_1b_city_box_maaslin3.png", p_categorical_city, width = 10, height = 5, dpi = 300)
ggsave("fig_1b_cutotype_box_maaslin3.png", p_categorical_cutotype, width = 10, height = 5, dpi = 300)


# Scatter plots for continuous variable

# A function for subplots
plot_scatter <- function(species_name, res_row) {
  # Extract relative abundance & metadata
  plot_data <- data.frame(
    Abundance = input_data[[species_name]], 
    Age = input_metadata[["Age"]]
  )
  
  # Annotate coefficient and q-value
  stats_text <- paste0("Coefficient = ", round(res_row$coef, 3), 
                       "\nq-value = ", format(res_row$qval_individual, scientific = TRUE, digits = 3))
  
  # Plot (subplot)
  ggplot(plot_data, aes(x = Age, y = Abundance)) +
    geom_point(alpha = 0.6, color = "#71b7ed") +
    geom_smooth(method = "lm", color = "#f57c6e", fill = "#f2656f") +
    scale_y_log10() +
    labs(title = species_name, 
         subtitle = stats_text,
         x = "Age (Years)",
         y = "Relative Abundance (Log10)") +
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
(p5 <- plot_scatter(top_age_pos_species, top_age_pos))
(p6 <- plot_scatter(top_age_neg_species, top_age_neg))

# Merge and save subplots
(p_continuous <- p5 + p6)
ggsave("figure_1b_age_scatter_maaslin3.png", p_continuous, width = 10, height = 5, dpi = 300)
