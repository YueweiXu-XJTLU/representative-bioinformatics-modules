# =========================
# Gel Analysis
# =========================

library(ggplot2)
library(dplyr)

# -----------------------------
# 1. Input data
# -----------------------------
data <- read.csv("gel_data.csv", header = TRUE)

# -----------------------------
# 2. Extract ladder
# -----------------------------
ladder <- data %>%
  filter(grepl("marker", band)) %>%
  mutate(
    bp = as.numeric(gsub("marker_", "", band)),
    distance = Pixcel
  ) %>%
  select(bp, distance)

# -----------------------------
# 3. Fit standard curve
# log10(bp) ~ distance
# -----------------------------
fit <- lm(log10(bp) ~ distance, data = ladder)

# R
(r_squared <- summary(fit)$r.squared)

# -----------------------------
# 4. Extract sample bands
# -----------------------------
sample_bands <- data %>%
  filter(!grepl("marker", band)) %>%
  filter(Pixcel > 0) %>%
  mutate(
    distance = Pixcel,
    sample = band
  ) %>%
  select(sample, distance)

# -----------------------------
# 5. Predict fragment size
# -----------------------------
(sample_bands <- sample_bands %>%
  mutate(
    sample_name = c("S3", "S5_1", "S5_2", "S6", "S7_1", "S7_2", "S8"), 
    bp_est = round(10^(predict(fit, newdata = .)))
  ))
write.csv(sample_bands, "bands_data.csv")
# Desired Result:
# lane_1_1: 1124 bp
# lane_4_1: 1023 bp
# lane_6_1: 1023 bp

# -----------------------------
# 6. Plot calibration
# -----------------------------
ladder_plot_df <- ladder %>%
  mutate(log10_bp = log10(bp))

(p_calibration <- ggplot() +
  geom_point(
    data = ladder_plot_df,
    aes(x = distance, y = log10_bp),
    color = "#f2656f",
    size = 2,
    alpha = 0.5
  ) +
  geom_smooth(
    data = ladder_plot_df,
    aes(x = distance, y = log10_bp),
    method = "lm",
    se = FALSE,
    color = "#f2656f",
  ) +
  geom_point(
    data = sample_bands[sample_bands$distance < 390,],
    aes(x = distance, y = log10(bp_est)),
    color = "#71b7ed",
    size = 3,
    shape = 17
  ) +
  geom_point(
    data = sample_bands[sample_bands$distance >= 390,],
    aes(x = distance, y = log10(bp_est)),
    color = "#71b7ed",
    size = 3,
    shape = 17,
    alpha = 0.6
  ) +
  labs(
    title = "Band Size Estimation Using Calibration Curve",
    subtitle = expression(R^2 == 0.97),
    x = "Migration Distance (Pixels)",
    y = expression(log[10](bp))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(face = "bold", hjust = 0.5),
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
)

ggsave("Fig/calibration.png", p_calibration, width = 6, height = 4, dpi = 300)
