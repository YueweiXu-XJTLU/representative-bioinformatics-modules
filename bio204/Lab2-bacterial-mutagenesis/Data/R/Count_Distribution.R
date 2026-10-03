library(tidyverse)
library(grid)

df_plate <- tibble(
  Plate = factor(1:12, levels = 1:12),
  Count = c(NA, 2, 8, 4, 260, 22, 28, 36, 8, 11, 225, 14)
)

p_plate <- ggplot(df_plate, aes(x = Plate, y = Count)) +
  geom_col(
    width = 0.72,
    fill = "#71b7ed",
    linewidth = 0.25,
    na.rm = TRUE
  ) +
  geom_text(
    data = df_plate %>% filter(!is.na(Count)),
    aes(label = Count),
    vjust = -0.35,
    size = 3.8,
    fontface = "bold"
  ) +
  annotate("text", x = 1, y = 8, label = "NA", fontface = "bold", size = 4) +
  labs(
    x = "Plate number",
    y = "Colony count",
    title = "Colony counts across Plates 1–12"
  ) +
  expand_limits(y = max(df_plate$Count, na.rm = TRUE) * 1.08) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
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

print(p_plate)

ggsave("Plate_colony_counts.png", p_plate, width = 8, height = 5, dpi = 300)
