# data preparation
observed_count <- c(43, 89, 59) # 100/100, 100/400, 400/400
expected_probability <- c(0.209869521, 0.496491324, 0.293639155) # 100/100, 100/400, 400/400

# chi-square test
chisq.test(
  x = observed_count,
  p = expected_probability,
  rescale.p = TRUE
)

# 	Chi-squared test for given probabilities
# data:  observed_count
# X-squared = 0.72187, df = 2, p-value = 0.697

library(ggplot2)
library(tidyr)
library(dplyr)

geno <- data.frame(
  Genotype = c("100/100", "100/400", "400/400"),
  Expected = c(40, 95, 56),
  Observed = c(43, 89, 59)
)

geno_long <- geno %>%
  pivot_longer(
    cols = c(Expected, Observed),
    names_to = "Type",
    values_to = "Count"
  )

p <- ggplot(geno_long, aes(x = Genotype, y = Count, fill = Type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("#71b7ed", "#f2656f")) +
  labs(
    title = "Observed and expected genotype counts at the TPA-25 locus",
    subtitle = expression(chi^2 == 0.722 * "," ~ p == 0.697),
    x = "Genotype",
    y = "Number of individuals",
    fill = "Count type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.margin = margin(15, 15, 15, 15),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.key.size = unit(0.5, "cm"),
    strip.background = element_rect(color = NA),
    strip.text = element_text(face = "bold")
  )
ggsave(
  "TPA-25.png",
  plot = p,
  width = 8,
  height = 4,
  dpi = 300
)
