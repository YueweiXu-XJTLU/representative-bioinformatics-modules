# if(!require("BiocManager", quietly=TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)

# Read the Gene expression file from Step 8 and rename the columns
sampleNames <- c("C1", "C2", "C3", "T1", "T2", "T3")
ED1 <- read.table("feature.txt", header = TRUE, quote = "\t", skip = 1) # Expression Data 1

# Extract data as a matrix
names(ED1)[c(7:12)] <- sampleNames
countData <- as.matrix(ED1[7:12])
rownames(countData) <- ED1$Geneid
# Remove the genes which did not expressed
countData <- countData[rowSums(countData) != 0,]

# Generate the colData which save the condition of samples
SRR_df <- data.frame(
  id = c("C1", "C2", "C3", "T1", "T2", "T3"), 
  condition = c("control", "control", "control", "treatment", "treatment", "treatment")
  )
# Matching check
all(colnames(countData) == SRR_df$id)

####################################################################################################

# Gene differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = SRR_df,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds)
res # (treatment / control)

sum(res$padj < 0.05, na.rm = TRUE)

####################################################################################################

# Plot a MA figure (Mean-Average Plot)
# fold change ~ expression level (counts)
# grey: not DEGs & blue: DEGs (discriminated by alpha (significant level/FDR))
plotMA(res, main = "DESeq2", ylim = c(-3, 3), alpha = 0.05)

# better version (FDR<0.05 & |LFC|≥1 highlighted)

# LFC: log2FoldChange
# NS: Not Significant
# FDR: False Discovery Rate

# prepare data
ma_df <- as.data.frame(res)
alpha  <- 0.05
lfc_cut <- 1

# classify
ma_df$dir <- with(ma_df,
                  ifelse(!is.na(padj) & padj < alpha & log2FoldChange >=  lfc_cut, "Up",
                         ifelse(!is.na(padj) & padj < alpha & log2FoldChange <= -lfc_cut, "Down", "NS")))
ma_df$dir <- factor(ma_df$dir, levels = c("Up","Down","NS"))

# counts
n_up   <- sum(ma_df$dir == "Up",   na.rm = TRUE)
n_down <- sum(ma_df$dir == "Down", na.rm = TRUE)

# MA plot
library(ggplot2)
p_ma <- ggplot(ma_df, aes(x = baseMean, y = log2FoldChange, color = dir)) +
  geom_hline(yintercept = c(-1, 0, 1),
             linetype = c(2,1,2), linewidth = c(0.3,0.4,0.3),
             color = c("grey70","grey60","grey70")) +
  geom_point(alpha = 0.9, size = 0.85, na.rm = TRUE) +
  scale_x_log10() +
  scale_color_manual(
    values = c(Up = "#f57c6e", Down = "#71b7ed", NS = "grey70"),
    breaks = c("Up","Down","NS"),
    labels = c(paste0("Up (", n_up, ")"),
               paste0("Down (", n_down, ")"),
               "NS"),
    name = NULL
  ) +
  labs(x = "Mean of normalized counts",
       y = "log2 fold-change (treatment vs control)",
       title = "DESeq2 MA Plot") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.98, 0.98), # 整张图里的一个相对位置(0-1)坐标
    legend.justification = c(1,1), # 把legend的这个位置对其上述坐标
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  guides(color = guide_legend(title = NULL,
                              override.aes = list(size = 2, alpha = 1))) +
  theme(legend.key.height = grid::unit(6, "pt"),
        legend.key.width  = grid::unit(6, "pt"))
print(p_ma)
# ggsave("MAplot.png", p_ma, width = 6, height = 4.2, dpi = 300)

###################################################################################################

# Merge the result with normalized reads counts
resdata <- merge(as.data.frame(res), # DESeq2 result
                 as.data.frame(counts(dds, normalized=TRUE)), # normalized raw counts (DESeq对每个gene在各sample中normalize过的表达量)
                 by= "row.names", sort=FALSE)
colnames(resdata)[1] <- "gene_id"

# Give the tag to indicate whether genes have difference
resdata$change <- as.factor(
  ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) > 1,
         # Check whether ID are same
         ifelse(resdata$log2FoldChange > 1, "Up", "Down"),
         "NoDiff"))

# Extract up-regulated and down-regulated genes
upRegu <- resdata[which(resdata$change == "Up"),]
downRegu <- resdata[which(resdata$change == "Down"),]

print(paste0("Up-regulated genes: ", nrow(upRegu)))
print(paste0("Down-regulated genes: ", nrow(downRegu)))

# Write up- and down-regulated gene into files
# sub(pattern, replacement, x): cut prefix for gene_id
downRegu$gene_id <- sub("gene:", "", downRegu$gene_id)
upRegu$gene_id <- sub("gene:", "", upRegu$gene_id)

# write.table(upRegu$gene_id, 
            # "upRegulateGene_0.05.txt",
            # col.names = F,
            # row.names = F,
            # quote = F)
# write.table(downRegu$gene_id,
            # "downRegulateGene_0.05.txt",
            # col.names = F,
            # row.names = F,
            # quote = F)

# better version
# ## ====== sort by padj and export ======
# # 上调：保留常用列，按 padj 升序（若相同，再按 |LFC| 降序）
# up_sorted <- subset(upRegu, !is.na(padj),
#                     select = c(gene_id, log2FoldChange, padj, baseMean))
# up_sorted <- up_sorted[order(up_sorted$padj, -abs(up_sorted$log2FoldChange)), ]
# # 同步去掉 gene: 前缀（如不需要可删此行）
# up_sorted$gene_id <- sub("^gene:", "", up_sorted$gene_id)
# 
# # 下调同理
# down_sorted <- subset(downRegu, !is.na(padj),
#                       select = c(gene_id, log2FoldChange, padj, baseMean))
# down_sorted <- down_sorted[order(down_sorted$padj, -abs(down_sorted$log2FoldChange)), ]
# down_sorted$gene_id <- sub("^gene:", "", down_sorted$gene_id)
# 
# # 导出：带统计信息的 TSV（便于审阅），以及仅基因清单（便于后续富集分析）
# write.table(up_sorted,
#             "upRegulateGene_0.05.sorted_by_padj.tsv",
#             sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(down_sorted,
#             "downRegulateGene_0.05.sorted_by_padj.tsv",
#             sep = "\t", quote = FALSE, row.names = FALSE)
# 
# write.table(up_sorted$gene_id,
#             "upRegulateGene_0.05.sorted_gene_list.txt",
#             col.names = FALSE, row.names = FALSE, quote = FALSE)
# write.table(down_sorted$gene_id,
#             "downRegulateGene_0.05.sorted_gene_list.txt",
#             col.names = FALSE, row.names = FALSE, quote = FALSE)

####################################################################################################

# Volcano plot
volcano <- ggplot(data=resdata,
                  aes(x = log2FoldChange,
                      y = -log10(padj),
                      color = change)) +
  geom_point(alpha = 0.8, size = 1) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  ggtitle("DESeq2 Valcano") +
  scale_color_manual(name = "",
                     values = c("green", "red", "black"),
                     limits = c("Up", "Down", "NoDiff")) +
  geom_vline(xintercept = c(-1, 1),
             lty = 2,
             col = "gray",
             lwd = 0.5) +
  geom_hline(yintercept = -log10(0.01),
             lty = 2,
             col = "gray",
             lwd = 0.5)

volcano

# better version
library(ggplot2)
library(grid)

## 1. 处理 padj = 0，避免 -log10(0) = -Inf
min_pos_padj <- min(resdata$padj[resdata$padj > 0], na.rm = TRUE)
resdata$padj[resdata$padj == 0] <- min_pos_padj * 0.1

## 2. 保证分组顺序：Up, Down, NoDiff
resdata$change <- factor(resdata$change,
                         levels = c("Up", "Down", "NoDiff"))

## 3. 阈值线：FDR = 0.05
alpha_thr <- 0.05

## 4. Volcano plot
volcano <- ggplot(
  resdata,
  aes(x = log2FoldChange,
      y = -log10(padj),
      color = change)
) +
  geom_point(alpha = 0.9, size = 0.9, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1),
             linetype = 2, color = "grey70") +
  geom_hline(yintercept = -log10(alpha_thr),
             linetype = 2, color = "grey70") +
  scale_color_manual(
    ## 这里的名字必须和 factor 水平一致：Up / Down / NoDiff
    values = c(
      Up     = "#f57c6e",
      Down   = "#71b7ed",
      NoDiff = "grey70"
    ),
    breaks = c("Up", "Down", "NoDiff"),
    ## label 想显示成 “NS” 没问题，这只是文字
    labels = c(
      paste0("Up (",   n_up,   ")"),
      paste0("Down (", n_down, ")"),
      "NS"
    ),
    name = NULL
  ) +
  labs(
    x     = "log2 fold-change (treatment vs control)",
    y     = "-log10(adjusted p-value)",
    title = "DESeq2 Volcano Plot"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    legend.position  = c(0.98, 0.98),
    legend.justification = c(1, 1),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    axis.title.x     = element_text(face = "bold"),
    axis.title.y     = element_text(face = "bold"),
    legend.key.height = grid::unit(6, "pt"),
    legend.key.width  = grid::unit(6, "pt")
  ) +
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(size = 2, alpha = 1)
    )
  )

print(volcano)
# ggsave("Volcano.png", volcano, width = 6, height = 4.2, dpi = 300)

####################################################################################################

# ==== PCA (VST-normalized counts) – styled like your MA/Volcano ====

# 1) VST 变换（用于样本间距离 / PCA）
vsd <- vst(dds, blind = FALSE)

# 2) 取矩阵（如需仅用高变基因，可解注释下两行）
mat <- assay(vsd)
# rv <- matrixStats::rowVars(mat); mat <- mat[order(rv, decreasing = TRUE)[1:500], ]

# 3) PCA
pca <- prcomp(t(mat), scale. = FALSE)
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

# 4) 整理数据框
df <- as.data.frame(pca$x[, 1:2])
df$sample    <- rownames(df)
df$condition <- colData(vsd)$condition[match(df$sample, rownames(colData(vsd)))]
df$condition <- factor(df$condition, levels = c("control","treatment"))

# 5) 绘图（配色/标题/轴标题加粗，图例右下角）
library(ggplot2)
p <- ggplot(df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(level = 0.68, linetype = 2, linewidth = 0.5, show.legend = FALSE) +
  {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(size = 3, max.overlaps = 100, show.legend = FALSE)
    } else {
      geom_text(vjust = -1, size = 3, show.legend = FALSE)
    }
  } +
  scale_color_manual(
    values = c(control = "#f57c6e", treatment = "#71b7ed"),
    breaks = c("control","treatment"),
    labels = c("Control","Treatment"),
    name   = NULL
  ) +
  labs(
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)"),
    title = "PCA"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid        = element_blank(),
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    axis.title.x      = element_text(face = "bold"),
    axis.title.y      = element_text(face = "bold"),
    legend.position      = c(0.98, 0.98),  # 右上角
    legend.justification = c(1, 1),        # 以右上角对齐
    legend.background    = element_blank(),
    legend.key           = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  guides(color = guide_legend(title = NULL,
                            override.aes = list(size = 2, alpha = 1))) +
  theme(legend.key.height = grid::unit(6, "pt"),
        legend.key.width  = grid::unit(6, "pt"))
print(p)
# ggsave("PCA.png", p, width = 6, height = 4.2, dpi = 300)
