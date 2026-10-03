# =========================
# VWA fragment size analysis
# =========================

# Load packages
library(ggplot2)
library(dplyr)

# -----------------------------
# 1. Input ladder data
# -----------------------------
ladder <- data.frame(
  bp       = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517),
  distance = c(461, 380, 323, 286, 250, 233, 219, 203, 179, 170, 150, 131)
)

# -----------------------------
# 2. Fit standard curve
# log10(bp) ~ distance
# -----------------------------
fit <- lm(log10(bp) ~ distance, data = ladder)
r <- sqrt(summary(fit)$r.squared)
r

# Ladder plot data
ladder_plot_df <- ladder %>%
  mutate(
    log10_bp = log10(bp),
    type = "Ladder"
  )

# -----------------------------
# 3. Input sample band distances
# -----------------------------
sample_bands <- data.frame(
  sample = c(
    "G3-1-1", "G3-1-2", "G3-2-1", "G3-2-2", "G3-3-1", "G3-3-2", "G3-4-1", "G3-4-2",
    "G4-1-1", "G4-1-2", "G4-2-1", "G4-2-2", "G4-3-1", "G4-3-2", "G4-4-1", "G4-4-2"
  ),
  distance = c(
    404, 404, 398, 409, 401, 413, 398, 403,
    414, 421, 409, 409, 410, 410, 411, 421
  )
)

# -----------------------------
# 4. Predict fragment size
# -----------------------------
sample_bands <- sample_bands %>%
  mutate(
    log10_bp   = predict(fit, newdata = .),
    bp_est     = 10^(log10_bp),
    bp_est_int = round(bp_est)
  )

# -----------------------------
# 5. Assign VWA allele size classes
# -----------------------------
vwa_classes <- c(138, 142, 146, 150, 154, 158, 162)

nearest_class <- function(x, classes) {
  classes[which.min(abs(classes - x))]
}

sample_bands <- sample_bands %>%
  mutate(
    size_class      = sapply(bp_est, nearest_class, classes = vwa_classes),
    int_minus_class = bp_est_int - size_class
  )

# -----------------------------
# 6. Output sample_df
# only requested columns
# -----------------------------
sample_df <- sample_bands %>%
  select(sample, distance, log10_bp, bp_est, bp_est_int, size_class, int_minus_class)

print(sample_df)

# -----------------------------
# 7. Plot ladder points + fit line + sample points
# -----------------------------
sample_plot_df <- sample_bands %>%
  mutate(type = "Sample")

plot_df <- bind_rows(
  ladder_plot_df %>% select(distance, log10_bp, type),
  sample_plot_df %>% select(distance, log10_bp, type)
)

p_calibration <- ggplot() +
  # ladder points
  geom_point(
    data = ladder_plot_df,
    aes(x = distance, y = log10_bp),
    color = "#f2656f",
    size = 2,
    alpha = 0.5
  ) +
  # regression line
  geom_smooth(
    data = ladder_plot_df,
    aes(x = distance, y = log10_bp),
    method = "lm",
    se = FALSE,
    color = "#71b7ed",
    linewidth = 0.7
  ) +
  labs(
    title = "DNA ladder calibration curve",
    subtitle = expression(R^2 == 0.999),
    x = "Migration distance (pixels)",
    y = expression(log[10](bp))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 12),
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

print(p_calibration)

ggsave(
  "ladder_calibration.png",
  plot = p_calibration,
  width = 6,
  height = 4,
  dpi = 300
)

# Scaled to 100-200bp
p_zoom <- ggplot() +
  # regression line
  geom_smooth(
    data = ladder_plot_df,
    aes(x = distance, y = log10_bp),
    method = "lm",
    se = FALSE,
    color = "#71b7ed",
    linewidth = 0.7
  ) +
  # sample points
  geom_point(
    data = sample_plot_df,
    aes(x = distance, y = log10_bp),
    color = "#f2656f",
    size = 2,
    shape = 17
  ) +
  labs(
    title = "Sample positions on the zoomed calibration curve",
    subtitle = "100–200 bp",
    x = "Migration distance (pixels)",
    y = expression(log[10](bp))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.margin = margin(15, 15, 15, 15),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.key.size = unit(0.5, "cm"),
    strip.background = element_rect(color = NA),
    strip.text = element_text(face = "bold")
  ) +
  coord_cartesian(
    ylim = log10(c(100, 200)),
    xlim = c(375, 475)
  )
print(p_zoom)
ggsave(
  "zoom.png",
  plot = p_zoom,
  width = 6,
  height = 4,
  dpi = 300
)

###################
# All Groups
###################

library(dplyr)

# ladder bp size
bp <- c(1517, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100)

# ladder pixel (每组一行)
ladder_pixels <- list(
  c(125, 144, 161, 171, 195, 208, 221, 240, 273, 309, 366, 458),
  c(131, 150, 170, 179, 203, 219, 233, 250, 286, 323, 380, 461),
  c(124, 141, 160, 172, 193, 207, 222, 239, 273, 313, 370, 458),
  c(127, 144, 163, 173, 196, 210, 223, 240, 276, 310, 367, 463),
  c(128, 146, 166, 176, 198, 212, 225, 242, 278, 315, 368, 455),
  c(128, 147, 167, 176, 199, 215, 227, 246, 279, 317, 371, 453),
  c(133, 153, 172, 184, 207, 224, 238, 256, 292, 331, 392, 478),
  c(132, 150, 171, 183, 202, 220, 234, 251, 288, 326, 383, 482),
  c(120, 136, 153, 164, 184, 200, 213, 228, 260, 297, 352, 438),
  c(121, 135, 153, 163, 183, 198, 210, 227, 260, 297, 348, 435),
  c(122, 140, 158, 169, 189, 204, 217, 235, 269, 304, 360, 449),
  c(129, 146, 165, 176, 199, 214, 226, 243, 278, 313, 370, 451),
  c(128, 145, 164, 174, 195, 210, 225, 239, 276, 312, 367, 450),
  c(125, 142, 160, 169, 191, 207, 219, 235, 271, 306, 366, 449)
)

# band pixels
bands <- list(
  c(381, 381, 381, 392, 396, 390, 390, 403, 389, 396, 397, 397, 401, 408),
  c(404, 404, 398, 409, 401, 413, 398, 403, 414, 421, 409, 409, 410, 410, 411, 421),
  c(390, 390, 393, 393, 392, 400, 409, 409, 394, 394, 397, 397, 383, 383, 380, 380),
  c(397, 406, 400, 416, 399, 406, 399, 412, 390, 406, 396, 396, 379, 382, 377, 390),
  c(404, 411, 394, 408, 399, 399, 402, 402, 409, 409, 396, 396, 407, 391, 399, 403),
  c(401, 408, 410, 410, 410, 410, 399, 405, 398, 398, 408, 408, 409, 409, 406, 406),
  c(414, 429, 416, 423, 413, 420, 423, 431, 418, 418, 420, 428, 419, 427, 418, 423),
  c(393, 410, 402, 402, 409, 416, 417, 417, 424, 412, 409, 422, 409, 409, 402, 402),
  c(371, 371, 375, 383, 382, 382, 381, 391, 375, 390, 387, 387, 391, 378, 388),
  c(363, 377, 367, 367, 373, 385, 368, 384, 369, 376, 375, 387, 381, 381),
  c(375, 375, 378, 386, 379, 385, 389, 394, 389, 389, 393, 393),
  c(387, 394, 399, 406, 396, 396, 394, 407, 406, 406, 390, 390, 399, 397),
  c(415, 415, 400, 400, 397, 411, 404, 404, 394, 394, 395, 395, 391, 391),
  c(378, 387, 385, 385, 384, 391, 384, 384, 386, 402, 389, 402, 382, 394, 385, 385)
)

# VWA allele classes
classes <- c(138, 142, 146, 150, 154, 158, 162)

nearest_class <- function(x) {
  classes[which.min(abs(classes - x))]
}

# 批量计算

result <- bind_rows(lapply(seq_along(bands), function(i) {
  # build ladder dataframe
  ladder <- data.frame(
    bp = bp,
    distance = ladder_pixels[[i]]
  )

  # fit calibration model
  fit <- lm(log10(bp) ~ distance, data = ladder)

  # sample bands
  sample_bands <- data.frame(
    distance = bands[[i]]
  )

  # predict fragment size
  sample_bands <- sample_bands %>%
    mutate(
      log10_bp = predict(fit, newdata = .),
      bp_est = 10^(log10_bp),
      bp_est_int = round(bp_est),
      size_class = sapply(bp_est, nearest_class),
      error = bp_est_int - size_class,
      group = i
    )

  sample_bands %>%
    select(group, pixel = distance, estimate = bp_est_int, class = size_class, error)
}))

result_filtered <- result %>%
  filter(
    estimate <= 164,
    if_else(class == 142, abs(error) < 2, abs(error) <= 8)
  )

print(result_filtered)

# -----------------------------
# 8. Frequency table
# -----------------------------
freq_table <- result_filtered %>%
  count(class, name = "Count") %>%
  arrange(class)

print(freq_table)

# -----------------------------
# 9. Frequency bar plot
# -----------------------------
p_freq <- ggplot(freq_table, aes(x = factor(class), y = Count)) +
  geom_col(
    fill = "#71b7ed",
    width = 0.7
  ) +
  labs(
    title = "Frequency distribution of VWA fragment size classes",
    x = "Fragment size class (bp)",
    y = "Count"
  ) +
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

print(p_freq)
ggsave(
  "frequency.png",
  plot = p_freq,
  width = 6,
  height = 4,
  dpi = 300
)
# -----------------------------
# 10. Summary statistics
# -----------------------------
min_size <- min(sample_bands$size_class, na.rm = TRUE)
max_size <- max(sample_bands$size_class, na.rm = TRUE)
n_alleles <- length(unique(sample_bands$size_class))

freq_vec <- table(sample_bands$size_class)
most_common <- names(freq_vec)[which.max(freq_vec)]
rarest <- names(freq_vec)[freq_vec == min(freq_vec)]

summary_stats <- list(
  model_coefficients = coef(fit),
  model_summary = summary(fit),
  min_size = min_size,
  max_size = max_size,
  n_alleles = n_alleles,
  most_common = most_common,
  rarest = rarest
)

print(summary_stats)
