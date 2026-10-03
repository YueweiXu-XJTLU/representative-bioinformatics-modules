# import raw data
df_raw <- read.csv("../../../shared-data/kidney_disease.csv")
df <- df_raw

# data clean
# fill gaps using NA, for later data process
df[df == ""] <- NA

# clean characters
df[] <- lapply(df, function(x){
  if (is.character(x)) {
    x <- trimws(x)
    x <- gsub("\\t", "", x)
    x[x %in% c("", "?", "NA")] <- NA
  }
  return(x)
})

# clean numerics
num_cols <- c("age","bp","sg","al","su","bgr","bu","sc","sod","pot","hemo","pcv","wc","rc")
for (col in num_cols) {
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}

# change all binomial variables into one-hot (0 & 1), for better ML model performance
df$htn <- ifelse(df$htn == "no", 0L, ifelse(df$htn == "yes", 1L, NA))
df$dm <- ifelse(df$dm == "no", 0L, ifelse(df$dm == "yes", 1L, NA))
df$cad <- ifelse(df$cad == "no", 0L, ifelse(df$cad == "yes", 1L, NA))
df$appet <- ifelse(df$appet == "poor", 0L, ifelse(df$appet == "good", 1L, NA))
df$pe <- ifelse(df$pe == "no", 0L, ifelse(df$pe == "yes", 1L, NA))
df$ane <- ifelse(df$ane == "no", 0L, ifelse(df$ane == "yes", 1L, NA))
df$pcc <- ifelse(df$pcc == "notpresent", 0L, ifelse(df$pcc == "present", 1L, NA))
df$ba <- ifelse(df$ba == "notpresent", 0L, ifelse(df$ba == "present", 1L, NA))
df$rbc <- ifelse(df$rbc == "abnormal", 0L, ifelse(df$rbc == "normal", 1L, NA))
df$pc <- ifelse(df$pc == "abnormal", 0L, ifelse(df$pc == "normal", 1L, NA))

# build factors (negative class is the reference level)
df$classification <- factor(df$classification, levels = c("notckd", "ckd"))

# calculate missing rate for rows
row_missing_table <- data.frame(
  row_id = 1:nrow(df),
  missing_count = rowSums(is.na(df)),
  missing_rate = rowSums(is.na(df)) / ncol(df)
)
row_missing_table <- row_missing_table[order(-row_missing_table$missing_rate), ]
head(row_missing_table)

# check number of row with over 30% cells missing
nrow(row_missing_table[row_missing_table$missing_count > (0.3 * ncol(df)), ])
dirty_rows <- row_missing_table[row_missing_table$missing_count > (0.3 * ncol(df)), ]
df <- df[-as.vector(dirty_rows$row_id), ]

# calculate missing rate for cols
col_missing_table <- data.frame(
  col_id = 1:ncol(df),
  missing_count = colSums(is.na(df)),
  missing_rate = colSums(is.na(df)) / nrow(df),
  missing_level = rep(NA, ncol(df))
)
for (i in 1:ncol(df)) {
  if (col_missing_table$missing_rate[i] < 0.05) {
    col_missing_table$missing_level[i] <- "low"
  } else if(col_missing_table$missing_rate[i] < 0.2) {
    col_missing_table$missing_level[i] <- "moderate"
  } else {
    col_missing_table$missing_level[i] <- "high"
  }
}
col_missing_table <- col_missing_table[order(-col_missing_table$missing_rate), ]
head(col_missing_table)

# fill methods (mannually pick out)
median_filled_col_low <- c("bu", "bp", "sc", "age")
mode_filled_col_low <- c("pcc", "ba", "htn", "dm", "cad", "appet", "pe", "ane")
median_filled_col_moderate <- c("pot", "sod", "pcv", "bgr", "hemo", "su", "sg", "al")
mode_filled_col_moderate <- c("pc")
mode_filled_col_high <- c("rbc", "rc", "wc")

# missing value filling
# a function for later mode calculation
calculate_mode <- function(x) {
  unique_x <- unique(x[!is.na(x)])
  return(unique_x[which.max(tabulate(match(x, unique_x)))])
}

# for low missing_level columns
for (i in median_filled_col_low) {
  df[[i]][is.na(df[[i]])] <- median(df[[i]], na.rm = TRUE)
}
for (i in mode_filled_col_low) {
  df[[i]][is.na(df[[i]])] <- calculate_mode(df[[i]])
}

# for moderate missing_level columns
for(i in median_filled_col_moderate){
  median_ckd <- median(df[df$classification == "ckd", i], na.rm = TRUE)
  median_notckd <- median(df[df$classification == "notckd", i], na.rm = TRUE)
  df[df$classification == "ckd"    & is.na(df[, i]), i] <- median_ckd
  df[df$classification == "notckd" & is.na(df[, i]), i] <- median_notckd
}

for (i in mode_filled_col_moderate) {
  mode_ckd <- calculate_mode(df[df$classification == "ckd", i])
  mode_notckd <- calculate_mode(df[df$classification == "notckd", i])
  df[df$classification == "ckd"    & is.na(df[, i]), i] <- mode_ckd
  df[df$classification == "notckd" & is.na(df[, i]), i] <- mode_notckd
}

# for high missing_level columns
for (i in mode_filled_col_high) {
  df[[paste0(i, "_missing")]] <- ifelse(is.na(df[[i]]), 1L, 0L)
  
  mode_ckd <- calculate_mode(df[df$classification == "ckd", i])
  mode_notckd <- calculate_mode(df[df$classification == "notckd", i])
  df[df$classification == "ckd"    & is.na(df[, i]), i] <- mode_ckd
  df[df$classification == "notckd" & is.na(df[, i]), i] <- mode_notckd
}

write.csv(df, "cleaned_kidney_disease.csv")

##########################################################################################################################################################################################################################################################
# RQ1
##########################################################################################################################################################################################################################################################
library(ggplot2)
library(dplyr)
library(e1071)

# 1. Load Data
df <- read.csv("cleaned_kidney_disease.csv")

# ==========================================
# Part A: Age Analysis (Descriptive Stats & Plot)
# ==========================================

# Calculate Statistics
age_stats <- df %>%
  summarise(
    Mean = mean(age, na.rm = TRUE),
    SD = sd(age, na.rm = TRUE),
    Median = median(age, na.rm = TRUE),
    Mode = as.numeric(names(sort(table(age), decreasing = TRUE)[1])),
    Min = min(age, na.rm = TRUE),
    Max = max(age, na.rm = TRUE),
    Skewness = skewness(age, na.rm = TRUE),
    Kurtosis = kurtosis(age, na.rm = TRUE)
  )
print(age_stats)

# Plot: Histogram with Density Curve
plot_age <- ggplot(df, aes(x = age)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 5, 
                 fill = "#93C5FD", color = "#f1f3fa", alpha = 0.7) +
  geom_density(color = "#3B82F6", size = 1.2) +
  labs(title = "Distribution of Age",
       x = "Age (years)",
       y = "Density") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "#f1f3fa", color = NA),
    panel.background = element_rect(fill = "#f1f3fa", color = NA),
    axis.title = element_text(face = "bold", size = 12, color = "#374151"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "#1E3A8A"),
    axis.text = element_text(color = "#374151"),
    panel.grid.major = element_line(color = "#d1d5db", size = 0.5),
    panel.grid.minor = element_blank()
  )
print(plot_age)
ggsave("RQ1_Age_Distribution.png", plot_age, width = 8, height = 6, dpi = 300)

# ==========================================
# Part B: Classification Analysis (Freq Table & Pie Plot)
# ==========================================

# Calculate Frequency Table
class_table <- df %>%
  group_by(classification) %>%
  summarise(Count = n()) %>%
  mutate(
    Percentage = Count / sum(Count) * 100,
    Label = paste0(classification, "\n", round(Percentage, 1), "%")
  )
print(class_table)

# Plot: Pie Chart
plot_pie <- ggplot(class_table, aes(x = "", y = Percentage, fill = classification)) +
  geom_bar(stat = "identity", width = 1, color = "#f1f3fa", linewidth = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("#3B82F6", "#93C5FD")) + # Consistent Blue Palette
  labs(title = "Proportion of CKD Diagnosis",
       fill = "Group") +
  theme_void() + # Removes axes and grids
  theme(
    plot.background = element_rect(fill = "#f1f3fa", color = NA),
    panel.background = element_rect(fill = "#f1f3fa", color = NA),
    plot.title = element_text(face = "bold", size = 16, color = "#1E3A8A", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", color = "#374151")
  ) +
  # Add text labels inside pie
  geom_text(aes(label = Label), 
            position = position_stack(vjust = 0.5), 
            color = "#f1f3fa", fontface = "bold", size = 5)
print(plot_pie)
ggsave("RQ1_Classification_Pie.png", plot_pie, width = 6, height = 6, dpi = 300)


##########################################################################################################################################################################################################################################################
# RQ2
##########################################################################################################################################################################################################################################################
notckd_hemo <- subset(df, classification == "ckd")$hemo
notnotckd_hemo <- subset(df, classification == "notckd")$hemo

# samples independence
# Independence is highly related to experimental design, hence we assume the observations are independent within the single group, data are also independent among different groups.

# sample normality (Shapiro-Wilk Test): due to large N = 369, we choose alpha = 0.001
# H0: Data is normally distributed
# H1: Data is not normally distributed
(ckd_shapiro <- shapiro.test(ckd_hemo)) # p-value = 0.001053 > 0.001, we CANNOT REJECT H0, hence assume data are normally distributed
(non_ckd_shapiro <- shapiro.test(notckd_hemo)) # p-value = 0.002518 > 0.001, we CANNOT REJECT H0, hence assume data are normally distributed

# variance equality  (Levene's Test)
# H0: The variances of the groups are equal.
# H1: At least one group has a different variance.
library(car)
merged_value <- c(
  ckd_hemo, 
  notckd_hemo
  )
merged_label <- c(
  rep('ckd', length(ckd_hemo)), 
  rep('notckd', length(notckd_hemo))
)
merged_data <- data.frame(
  value = merged_value,
  label = merged_label,
  stringsAsFactors = TRUE
)
(ckd_levene <- leveneTest(value ~ label, data = merged_data, center = median)) # p-value = 0.0003288 < 0.05, we REJECT H0, hence CANNOT REJECT H1: At least one group has a different variance.

# t.test
(ckd_t_test <- t.test(notckd_hemo, ckd_hemo, alternative = "two.sided", conf.level = 0.95, var.equal = FALSE))
# t = 25.987, df = 362.72, p-value < 2.2e-16
# 95 percent confidence interval:
#   4.812369 4.135266

##########
# Hypotheses: 
#   H0: there is no significant difference in means of hemoglobin levels between CKD & non-CKD patient.
# Pre-evaluation: 
#   Independence: related to experimental design & assume true
# Normality (Shapiro-Wilk Test): 
#   large N (150 & 219)  CLT holds  use α = 0.001
# p-value = 0.0011 (CKD) & 0.0025 (non-CKD) < 0.001  H0  both data are normally distributed
# Variance Equality (Levene's Test):
# p-value = 0.0003288 < 0.05  H1  variances are un-equal
# Two Samples T-Test (non_CKD - CKD):
# t = 25.987 (df = 362.72)
# 95 CI: [4.81, 4.14] (g/dL)
# p-value < 2.2e-16  H1 
# Conclusion: There is a significant difference in means of hemoglobin levels between CKD & non-CKD patient.  Lower hemoglobin levels are a potential pathological characteristic of CKD.
##########

# violin plot
library(tidyverse)
library(ggplot2)
library(hrbrthemes)

my_colors <- c("ckd" = "#3B82F6", "notckd" = "#93C5FD")
p <- ggplot(df, aes(x = classification, y = hemo, fill = classification)) +
  geom_violin(
    trim = FALSE,            # 保留尾部数据，展示完整分布
    alpha = 0.8,             # 透明度适中，让内部线条可见
    color = NA,              # 移除小提琴最外圈的轮廓线，更干净
    quantiles = c(0.25, 0.5, 0.75),   # 指定三条线的位置
    quantile.linetype = "dashed",     # 【关键】设置为虚线
    quantile.color = "grey50",       # 线条颜色：深灰蓝 (比纯黑更柔和)
    linewidth = 0.6                   # 线条稍微细一点，更精致
  ) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "Comparative Distribution of Hemoglobin Levels (ckd VS. notckd)",
    x = NULL,
    y = "Hemoglobin (g/dL)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "#f1f3fa", color = NA),
    panel.background = element_rect(fill = "#f1f3fa", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "#374151", margin = margin(b = 10)),
    axis.title.y = element_text(face = "bold", color = "#374151", margin = margin(r = 10)),
    axis.text = element_text(color = "#374151"),
    panel.grid.major.y = element_line(color = "#d1d5db", linewidth = 0.8),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )
p
ggsave("RQ2_violin_hemo_academic.png", plot = p, width = 8, height = 6, dpi = 300)

##########################################################################################################################################################################################################################################################
# RQ3
##########################################################################################################################################################################################################################################################
library(tidyverse)
library(ggpubr)  # 用于绘制学术级交互作用图
library(car)     # 用于方差齐性检验 (Levene's Test)

# samples independence
# Independence is highly related to experimental design, hence we assume the observations are independent within the single group, data are also independent among different groups.

# sample normality (Shapiro-Wilk Test): due to large N = 369, we choose alpha = 0.001
# H0: Data is normally distributed
# H1: Data is not normally distributed
a <- df$sc[df$htn == 0 & df$dm == 0]
b <- df$sc[df$htn == 0 & df$dm == 1]
c <- df$sc[df$htn == 1 & df$dm == 0]
d <- df$sc[df$htn == 1 & df$dm == 1]
(la <- length(a)) # 213 OK
(lb <- length(b)) # 24
(lc <- length(c)) # 36
(ld <- length(d)) # 96 OK
(shapiro.test(log(a))) # p-value = 7.353e-16
(shapiro.test(log(b))) # p-value = 0.003832 OK
(shapiro.test(log(c))) # p-value = 0.1635 OK
(shapiro.test(log(d))) # p-value = 0.002434

# variance equality  (Levene's Test)
# H0: The variances of the groups are equal.
# H1: At least one group has a different variance.
library(car)
merged_value_aov <- c(
  a, b, c, d
)
merged_label_aov <- c(
  rep('00', la), 
  rep('01', lb),
  rep('10', lc),
  rep('11', ld)
)
merged_data_aov <- data.frame(
  value = merged_value_aov,
  label = merged_label_aov,
  stringsAsFactors = TRUE
)
(ckd_levene_aov <- leveneTest(value ~ label, data = merged_data_aov, center = median)) # p-value = 0.0002992 < 0.05, we REJECT H0, hence CANNOT REJECT H1: At least one group has a different variance.

# 公式 sc ~ htn * dm 等价于 sc ~ htn + dm + htn:dm
# 它可以同时检验：
# 1. 高血压是否影响 sc ?
# 2. 糖尿病是否影响 sc ?
# 3. 高血压和糖尿病对 sc 是否有交互作用 (例如: 既有高血压又有糖尿病的人, sc 升高得更剧烈)?
# 对数变换分析 (针对 sc 这种右偏数据)
ckd_aov <- aov(log(sc) ~ htn * dm, data = df)
Anova(ckd_aov, type = 3, white.adjust = TRUE)
# Df       F    Pr(>F)    
# (Intercept)   1  1.0441 0.3075386    
# htn           1 67.4402 3.800e-15 ***
#   dm            1 16.3410 6.456e-05 ***
#   htn:dm        1 11.0180 0.0009932 ***
#
# recheck
shapiro.test(residuals(ckd_aov)) # p-value = 2.753e-14

# Interactive EffectPlot
merged_data_aov$htn <- factor(ifelse(substr(merged_data_aov$label, 1, 1) == "0", "Not HTN", "HTN"), levels = c("Not HTN", "HTN"))
merged_data_aov$dm  <- factor(ifelse(substr(merged_data_aov$label, 2, 2) == "0", "Not DM", "DM"), levels = c("Not DM", "DM"))
p <- ggplot(merged_data_aov, aes(x = htn, y = log(value), color = dm, group = dm)) +
  stat_summary(fun = mean, geom = "point", size = 6) +
  stat_summary(fun = mean, geom = "line", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  scale_color_manual(values = c("#3B82F6", "#93C5FD")) +
  labs(
    title = "Interaction Plot of HTN and DM on log(SC)",
    x = "Hypertension (HTN)",
    y = "Mean log(Serum Creatinine) ± SE",
    color = "Diabetes Mellitus (DM)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "#f1f3fa", color = NA),
    panel.background = element_rect(fill = "#f1f3fa", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "#374151"),
    axis.title = element_text(face = "bold", color = "#374151"),
    axis.text = element_text(color = "#374151"),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "#d1d5db", linewidth = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
p
ggsave("RQ3_interaction_plot.png", width = 8, height = 6, dpi = 300)
# Interpretation: The concurrent presence of hypertension and diabetes has a synergistic impact on renal impairment, suggesting that comorbidity accelerates kidney dysfunction beyond additive effects .

# 效应分析
# H0: µ_no_constraint = µ_constriant
library(emmeans)
simple_dm <- emmeans(ckd_aov, pairwise ~ dm | htn)
summary(simple_dm)
# $contrasts
# htn = 0:
#   contrast  estimate    SE  df t.ratio p.value
# dm0 - dm1  -0.8423 0.157 365  -5.375  <.0001
# 
# htn = 1:
#   contrast  estimate    SE  df t.ratio p.value
# dm0 - dm1   0.0226 0.142 365   0.159  0.8741
# 
# simple_htn <- emmeans(ckd_aov, pairwise ~ htn | dm)
# summary(simple_htn)


##########################################################################################################################################################################################################################################################
# RQ4
##########################################################################################################################################################################################################################################################
library(ggpubr)

# 去除过高sc的1个离群值
df_filtered <- df[df$sc <= 50, ]

# 拟合模型
x <- "bu"
y <- "sc"
formula_str <- as.formula(paste(y, "~", x))
model <- lm(formula_str, data = df_filtered)
summary(model)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.656272   0.175532  -3.739 0.000215 ***
#   bu           0.059806   0.002367  25.267  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.236 on 366 degrees of freedom
# Multiple R-squared:  0.6356,	Adjusted R-squared:  0.6346 
# F-statistic: 638.4 on 1 and 366 DF,  p-value: < 2.2e-16
# y = 0.060 x - 0.66
# R = 0.8
# p < 2.2e-16

# 绘图
p <- ggplot(df_filtered, aes(x = .data[[x]], y = .data[[y]])) +
  # 1. 散点图
  geom_point(color = "#84a6fc", alpha = 0.5, size = 2) +
  # 2. 回归拟合线
  geom_smooth(method = "lm", color = "#3B82F6", fill = "#93C5FD", alpha = 0.3) +
  # 3. 添加方程和 R 值
  # stat_regline_equation(label.x = 150, label.y = 60, color = "#374151") +
  # stat_cor(label.x = 150, label.y = 55, color = "#374151") +
  # 4. 标签
  labs(
    title = "Linear Regression Analysis: SC VS. BU",
    x = "Blood Urea (mg/dL)",
    y = "Serum Creatinine (mg/dL)"
  ) +
  # 5. 主题
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "#f1f3fa", color = NA),
    panel.background = element_rect(fill = "#f1f3fa", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "#374151"),
    axis.title = element_text(face = "bold", color = "#374151"),
    axis.text = element_text(color = "#374151"),
    panel.grid.major.y = element_line(color = "#d1d5db", linewidth = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
p
ggsave("RQ4_linear_regression.png", width = 8, height = 6, dpi = 300)
