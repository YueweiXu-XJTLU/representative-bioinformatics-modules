library(tidyverse)
library(grid)

# raw data
df <- tibble(
  ring = 0:9,
  plate7_accu = c(0, 0, 4, 9, 19, 27, 28, 28, 28, 28),
  plate8_accu = c(0, 0, 2, 11, 24, 28, 33, 34, 35, 36)
)

# derive newly added colonies in each ring
df <- df %>%
  mutate(
    plate7_new = c(plate7_accu[1], diff(plate7_accu)),
    plate8_new = c(plate8_accu[1], diff(plate8_accu))
  )

# long data for cumulative plot
accu_df <- df %>%
  select(ring, plate7_accu, plate8_accu) %>%
  pivot_longer(
    cols = c(plate7_accu, plate8_accu),
    names_to = "plate",
    values_to = "colonies"
  ) %>%
  mutate(
    plate = recode(
      plate,
      plate7_accu = "1% MMS",
      plate8_accu = "2% MMS"
    )
  )

# long data for new colonies plot
new_df <- df %>%
  select(ring, plate7_new, plate8_new) %>%
  pivot_longer(
    cols = c(plate7_new, plate8_new),
    names_to = "plate",
    values_to = "colonies"
  ) %>%
  mutate(
    plate = recode(
      plate,
      plate7_new = "1% MMS",
      plate8_new = "2% MMS"
    )
  )

# colours
mms_cols <- c(
  "1% MMS" = "#71b7ed",
  "2% MMS" = "#f2656f"
)

# x-axis labels
ring_labels <- c(
  "1\n(disc)",
  "2\n(inhibition\nzone)",
  "3", "4", "5", "6", "7", "8", "9", "10"
)

# common theme
base_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.margin = margin(15, 15, 15, 15),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.key.size = unit(0.5, "cm"),
    strip.background = element_rect(color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()   # remove minor grid
  )

# dodge width for alignment
dodge_w <- 0.75

# 1. newly observed colonies in each ring
p_new <- ggplot(new_df, aes(x = factor(ring), y = colonies, fill = plate)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, alpha = 0.14) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = 0.14) +
  geom_col(
    position = position_dodge(width = dodge_w),
    width = 0.68,
    linewidth = 0.2
  ) +
  scale_fill_manual(values = mms_cols) +
  scale_x_discrete(labels = ring_labels) +
  labs(
    x = "Distance from disc centre (disc-width units)",
    y = "Newly observed revertant colonies",
    title = "Spatial distribution of newly observed revertants in the MMS spot assay"
  ) +
  base_theme

# 2. cumulative colonies
p_accu <- ggplot(accu_df, aes(x = ring, y = colonies, colour = plate, group = plate)) +
  annotate("rect", xmin = -0.5, xmax = 0.5, ymin = -Inf, ymax = Inf, alpha = 0.14) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, alpha = 0.14) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = mms_cols) +
  scale_x_continuous(
    breaks = 0:9,
    labels = ring_labels
  ) +
  labs(
    x = "Distance from disc centre (disc-width units)",
    y = "Cumulative revertant colonies",
    title = "Cumulative spatial pattern of revertants in the MMS spot assay"
  ) +
  base_theme

print(p_new)
print(p_accu)

# save
ggsave("MMS_new_revertants.png", p_new, width = 10, height = 5, dpi = 300)
ggsave("MMS_cumulative_revertants.png", p_accu, width = 10, height = 5, dpi = 300)
