library(tidyverse)

# read data
df <- read.csv("GC_Content.csv", check.names = FALSE)

# remove leading/trailing spaces from column names
colnames(df) <- trimws(colnames(df))

# convert wide to long format
plot_df <- df %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sequence",
    values_to = "GC_content"
  ) %>%
  mutate(
    Sequence = factor(
      Sequence,
      levels = c("AN8455-R", "AN8455-R-NS-S_Original", "AN8455-R-NS-S_Paired",
                 "AN8455-R-NS-L_Original", "AN8455-R-NS-L_Paired")
    )
  )

# bar plot
(p <- ggplot(plot_df, aes(x = Sequence, y = GC_content)) +
  geom_col(fill = "#71b7ed", color = "#71b7ed", width = 0.7) +
  geom_text(aes(label = scales::percent(GC_content, accuracy = 1)),
            vjust = -0.4, size = 4, fontface = "bold") +
  scale_y_continuous(
    limits = c(0, 0.7),
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = "GC content",
    title = "GC content of AN8455-R and putative non-specific annealing sites"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.margin = margin(15, 15, 15, 15),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold", size = 8),
    legend.position = "right",
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.key.size = unit(0.5, "cm"),
    strip.background = element_rect(color = NA),
    strip.text = element_text(face = "bold")
  ))

(pp <- ggplot(plot_df, aes(x = Sequence, y = GC_content)) +
    geom_col(fill = "#71b7ed", color = "#71b7ed", width = 0.7) +
    geom_text(aes(label = scales::percent(GC_content, accuracy = 1)),
              vjust = -0.4, size = 4, fontface = "bold") +
    scale_y_continuous(
      limits = c(0, 0.7),
      expand = c(0, 0),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = NULL,
      y = "GC content",
      title = "GC content of AN8455-R and putative non-specific annealing sites"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.margin = margin(15, 15, 15, 15),
      axis.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(face = "bold", size = 8),
      axis.text.y = element_text(face = "bold", size = 8),
      legend.position = "right",
      legend.text = element_text(face = "bold", size = 12),
      legend.title = element_text(face = "bold", size = 14),
      legend.key.size = unit(0.5, "cm"),
      strip.background = element_rect(color = NA),
      strip.text = element_text(face = "bold")
    ))

ggsave("GC_Content.png", plot = p, width = 8, height = 4, dpi = 300)
ggsave("GC_Content_Poster.png", plot = pp, width = 8, height = 4, dpi = 300)
