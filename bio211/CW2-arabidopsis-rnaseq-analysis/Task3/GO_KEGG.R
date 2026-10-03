# if(!require("BiocManager", quietly=TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)

# prepare

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

# Gene differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = SRR_df,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds)
res # (treatment / control)

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
resdata$change <- factor(resdata$change, levels = c("Up","Down","NoDiff"))

###################################################################################################

## =========================================================
## Enrichment analysis (GO BP/CC/MF + KEGG) and plotting
## =========================================================

# ---- Packages ----
# need <- c("clusterProfiler","enrichplot","org.At.tair.db",
#           "ggplot2","dplyr","stringr")
# inst <- need[!need %in% installed.packages()[,1]]
# if (length(inst)) {
#   if (!"BiocManager" %in% installed.packages()[,1]) install.packages("BiocManager")
#   BiocManager::install(inst, ask = FALSE, update = FALSE)
# }

library(clusterProfiler)
library(enrichplot)
library(org.At.tair.db)
library(ggplot2)
library(dplyr)
library(stringr)

# ---- GO enrichment ----
# ---- helpers: reorder_within/scale_y_reordered (no tidytext needed) ----
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) gsub(paste0(sep, ".+$"), "", x), ...)
}

# ---- Inputs from resdata (already built above) ----
# 去掉 "gene:" 前缀，转为 TAIR/AGI 号 (ATxGxxxxx)
all_ids <- sub("^gene:", "", resdata$gene_id)

# DEG 标记   FDR<0.05; |LFC|≥1
deg_flag <- !is.na(resdata$padj) & resdata$padj < 0.05 & abs(resdata$log2FoldChange) >= 1
up_ids   <- sub("^gene:", "", resdata$gene_id[deg_flag & resdata$log2FoldChange >  0])
down_ids <- sub("^gene:", "", resdata$gene_id[deg_flag & resdata$log2FoldChange <  0])
deg_all  <- unique(c(up_ids, down_ids))

# 背景基因：建议取所有进入差异检验/成功注释的基因
universe_ids <- unique(all_ids)

# ---- GO enrichment (BP/CC/MF; keyType=TAIR) ----
ego_bp <- enrichGO(gene          = deg_all,
                   universe      = universe_ids,
                   OrgDb         = org.At.tair.db,
                   keyType       = "TAIR",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05, qvalueCutoff = 0.05)
ego_cc <- enrichGO(gene          = deg_all,
                   universe      = universe_ids,
                   OrgDb         = org.At.tair.db,
                   keyType       = "TAIR",
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05, qvalueCutoff = 0.05)
ego_mf <- enrichGO(gene          = deg_all,
                   universe      = universe_ids,
                   OrgDb         = org.At.tair.db,
                   keyType       = "TAIR",
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05, qvalueCutoff = 0.05)

# 去冗余（语义相似度）(simplify output from enrichGO and gseGO by removing redundancy of enriched GO terms)
ego_bp <- tryCatch(simplify(ego_bp, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) ego_bp)
ego_cc <- tryCatch(simplify(ego_cc, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) ego_cc)
ego_mf <- tryCatch(simplify(ego_mf, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) ego_mf)

# ---- Plot: GO barplot (BP/CC/MF faceted, Top 10 each by GeneRatio) ----
df_bp <- transform(as.data.frame(ego_bp), ONTO = "BP")
df_cc <- transform(as.data.frame(ego_cc), ONTO = "CC")
df_mf <- transform(as.data.frame(ego_mf), ONTO = "MF")
g_all <- dplyr::bind_rows(df_bp, df_cc, df_mf)

# 计算 GeneRatio 数值
g_all$GeneRatio_num <- if (is.numeric(g_all$GeneRatio)) g_all$GeneRatio else
  vapply(strsplit(as.character(g_all$GeneRatio), "/"),
         function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))

# 各 ONTO 选 n 个（按 GeneRatio 由高到低）
g_bp_top <- g_all[g_all$ONTO == "BP", ] |>
  dplyr::arrange(dplyr::desc(GeneRatio_num)) |>
  dplyr::slice_head(n = 5)
g_cc_top <- g_all[g_all$ONTO == "CC", ] |>
  dplyr::arrange(dplyr::desc(GeneRatio_num)) |>
  dplyr::slice_head(n = 5)
g_mf_top <- g_all[g_all$ONTO == "MF", ] |>
  dplyr::arrange(dplyr::desc(GeneRatio_num)) |>
  dplyr::slice_head(n = 5)
g_top <- dplyr::bind_rows(g_bp_top, g_cc_top, g_mf_top)

# y 轴文字换行美化
g_top$Description <- stringr::str_wrap(g_top$Description, width = 28)

p_go_bar <- ggplot(g_top,
                   aes(x = Count,
                       y = reorder_within(Description, GeneRatio_num, ONTO),
                       fill = p.adjust)) +
  geom_col() +
  scale_y_reordered() +
  facet_grid(ONTO ~ ., scales = "free_y", switch = "y") +
  scale_fill_gradient(name = "p.adjust",
                      low = "#cfe8ff", high = "#084a91",
                      trans = "reverse") +
  labs(x = "Count", y = NULL, title = "GO enrichment") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey70", color = NA),
        plot.title  = element_text(hjust = .5, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.y  = element_text(face = "bold"))
p_go_bar
ggsave("GO.png", p_go_bar, width = 7, height = 8, dpi = 300)

# # ---- Plot: GO dotplot (BP/CC/MF merged, Top 10 overall) ----
# p_bp <- enrichplot::dotplot(ego_bp, showCategory = 10) +
#   labs(x = "GeneRatio", y = NULL, title = "GO BP enrichment") +
#   theme_bw(base_size = 12) + theme(panel.grid = element_blank(),
#                                    plot.title = element_text(hjust = .5, face = "bold"),
#                                    axis.title.x = element_text(face = "bold"))
# # ggsave("GO_BP_dotplot.png", p_bp, width = 8, height = 5, dpi = 300)
# 
# p_cc <- enrichplot::dotplot(ego_cc, showCategory = 10) +
#   labs(x = "GeneRatio", y = NULL, title = "GO CC enrichment") +
#   theme_bw(base_size = 12) + theme(panel.grid = element_blank(),
#                                    plot.title = element_text(hjust = .5, face = "bold"),
#                                    axis.title.x = element_text(face = "bold"))
# # ggsave("GO_CC_dotplot.png", p_cc, width = 8, height = 5, dpi = 300)
# 
# p_mf <- enrichplot::dotplot(ego_mf, showCategory = 10) +
#   labs(x = "GeneRatio", y = NULL, title = "GO MF enrichment") +
#   theme_bw(base_size = 12) + theme(panel.grid = element_blank(),
#                                    plot.title = element_text(hjust = .5, face = "bold"),
#                                    axis.title.x = element_text(face = "bold"))
# p_bp
# p_cc
# p_mf
# # ggsave("GO_dotplot.png", p_go_dot, width = 10, height = 7, dpi = 300)

# ---- KEGG enrichment (organism = 'ath'; needs ENTREZ) ----
# 若已定义 resdata：用 FDR<0.05 的 DEGs（不再额外卡 |LFC|）
deg_all <- sub("^gene:", "", resdata$gene_id[!is.na(resdata$padj) & resdata$padj < 0.05])

# 映射 TAIR → ENTREZ
map_all <- bitr(deg_all, fromType = "TAIR", toType = "ENTREZID", OrgDb = org.At.tair.db)
entrez  <- unique(na.omit(map_all$ENTREZID))

# KEGG 富集（注意 keyType 用 "ncbi-geneid" 来匹配 ENTREZID）
ekegg <- enrichKEGG(
  gene          = entrez,
  organism      = "ath",
  keyType       = "ncbi-geneid",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# 查富集条目数
n <- nrow(as.data.frame(ekegg))
n

# Dotplot
p_kegg <- enrichplot::dotplot(ekegg, showCategory = n) +
  scale_color_continuous(name = "p.adjust", guide = "legend") +
  labs(x = "GeneRatio", y = NULL, title = "KEGG pathways") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title  = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.y  = element_text(face = "bold"))
p_kegg
# ggsave("KEGG.png", p_kegg, width = 10, height = 5, dpi = 300)

# ---- Quick console summary ----
cat(sprintf("GO terms: BP=%d, CC=%d, MF=%d; KEGG=%s\n",
            nrow(as.data.frame(ego_bp)),
            nrow(as.data.frame(ego_cc)),
            nrow(as.data.frame(ego_mf)),
            if (is.null(ekegg)) "NA" else as.character(nrow(as.data.frame(ekegg)))))
