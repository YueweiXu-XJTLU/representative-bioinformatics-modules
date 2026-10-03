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
## Enrichment analysis (GO BP/CC/MF + KEGG) for Up/Down DEGs
## =========================================================

library(clusterProfiler)
library(enrichplot)
library(org.At.tair.db)
library(ggplot2)
library(dplyr)
library(stringr)

## ---------- Helper functions for GO barplots ----------

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(
    labels = function(x) gsub(paste0(sep, ".+$"), "", x),
    ...
  )
}

## ---------- Define DEG sets and background ----------

# Remove "gene:" prefix to get AGI codes (ATxGxxxxx)
all_ids <- sub("^gene:", "", resdata$gene_id)

# DEGs: FDR < 0.05 and |log2FC| >= 1
deg_flag <- !is.na(resdata$padj) &
  resdata$padj < 0.05 &
  abs(resdata$log2FoldChange) >= 1

up_ids   <- sub("^gene:", "", resdata$gene_id[deg_flag & resdata$log2FoldChange >  0])
down_ids <- sub("^gene:", "", resdata$gene_id[deg_flag & resdata$log2FoldChange <  0])

# Background: all tested / annotated genes
universe_ids <- unique(all_ids)

## =========================================================
## GO enrichment (BP / CC / MF) – Up- and Down-regulated
## =========================================================

## ----- Up-regulated DEGs -----

ego_bp_up <- enrichGO(
  gene          = up_ids,
  universe      = universe_ids,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_cc_up <- enrichGO(
  gene          = up_ids,
  universe      = universe_ids,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_mf_up <- enrichGO(
  gene          = up_ids,
  universe      = universe_ids,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Reduce redundancy
ego_bp_up <- tryCatch(simplify(ego_bp_up, cutoff = 0.7, by = "p.adjust",
                               select_fun = min), error = function(e) ego_bp_up)
ego_cc_up <- tryCatch(simplify(ego_cc_up, cutoff = 0.7, by = "p.adjust",
                               select_fun = min), error = function(e) ego_cc_up)
ego_mf_up <- tryCatch(simplify(ego_mf_up, cutoff = 0.7, by = "p.adjust",
                               select_fun = min), error = function(e) ego_mf_up)

# Combine BP/CC/MF for plotting (Up)
df_bp_up <- transform(as.data.frame(ego_bp_up), ONTO = "BP")
df_cc_up <- transform(as.data.frame(ego_cc_up), ONTO = "CC")
df_mf_up <- transform(as.data.frame(ego_mf_up), ONTO = "MF")
g_up_all <- dplyr::bind_rows(df_bp_up, df_cc_up, df_mf_up)

# Numeric GeneRatio
g_up_all$GeneRatio_num <- if (is.numeric(g_up_all$GeneRatio)) {
  g_up_all$GeneRatio
} else {
  vapply(
    strsplit(as.character(g_up_all$GeneRatio), "/"),
    function(x) as.numeric(x[1]) / as.numeric(x[2]),
    numeric(1)
  )
}

# Top 5 terms per ontology
g_bp_up_top <- g_up_all[g_up_all$ONTO == "BP", ] |>
  dplyr::arrange(dplyr::desc(GeneRatio_num)) |>
  dplyr::slice_head(n = 5)

g_cc_up_top <- g_up_all[g_up_all$ONTO == "CC", ] |>
  dplyr::arrange(dplyr::desc(GeneRatio_num)) |>
  dplyr::slice_head(n = 5)

g_mf_up_top <- g_up_all[g_up_all$ONTO == "MF", ] |>
  dplyr::arrange(dplyr::desc(GeneRatio_num)) |>
  dplyr::slice_head(n = 5)

g_up_top <- dplyr::bind_rows(g_bp_up_top, g_cc_up_top, g_mf_up_top)
g_up_top$Description <- stringr::str_wrap(g_up_top$Description, width = 28)

p_go_up <- ggplot(
  g_up_top,
  aes(x = Count,
      y = reorder_within(Description, GeneRatio_num, ONTO),
      fill = p.adjust)
) +
  geom_col() +
  scale_y_reordered() +
  facet_grid(ONTO ~ ., scales = "free_y", switch = "y") +
  scale_fill_gradient(
    name = "p.adjust",
    low  = "#cfe8ff",
    high = "#084a91",
    trans = "reverse"
  ) +
  labs(x = "Count", y = NULL, title = "GO enrichment (Up-regulated DEGs)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    strip.background = element_rect(fill = "grey70", color = NA),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.title.x    = element_text(face = "bold"),
    axis.text.y     = element_text(face = "bold")
  )

p_go_up
ggsave("GO_up.png", p_go_up, width = 7, height = 8, dpi = 300)

## ----- Down-regulated DEGs -----

ego_bp_down <- enrichGO(
  gene          = down_ids,
  universe      = universe_ids,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_cc_down <- enrichGO(
  gene          = down_ids,
  universe      = universe_ids,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_mf_down <- enrichGO(
  gene          = down_ids,
  universe      = universe_ids,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_bp_down <- tryCatch(simplify(ego_bp_down, cutoff = 0.7, by = "p.adjust",
                                 select_fun = min), error = function(e) ego_bp_down)
ego_cc_down <- tryCatch(simplify(ego_cc_down, cutoff = 0.7, by = "p.adjust",
                                 select_fun = min), error = function(e) ego_cc_down)
ego_mf_down <- tryCatch(simplify(ego_mf_down, cutoff = 0.7, by = "p.adjust",
                                 select_fun = min), error = function(e) ego_mf_down)

df_bp_down <- transform(as.data.frame(ego_bp_down), ONTO = "BP")
df_cc_down <- transform(as.data.frame(ego_cc_down), ONTO = "CC")
df_mf_down <- transform(as.data.frame(ego_mf_down), ONTO = "MF")
g_down_all <- dplyr::bind_rows(df_bp_down, df_cc_down, df_mf_down)

g_down_all$GeneRatio_num <- if (is.numeric(g_down_all$GeneRatio)) {
  g_down_all$GeneRatio
} else {
  vapply(
    strsplit(as.character(g_down_all$GeneRatio), "/"),
    function(x) as.numeric(x[1]) / as.numeric(x[2]),
    numeric(1)
  )
}

g_bp_down_top <- g_down_all[g_down_all$ONTO == "BP", ] |>
  dplyr::arrange(dplyr::desc(GeneRatio_num)) |>
  dplyr::slice_head(n = 5)

g_cc_down_top <- g_down_all[g_down_all$ONTO == "CC", ] |>
  dplyr::arrange(dplyr::desc(GeneRatio_num)) |>
  dplyr::slice_head(n = 5)

g_mf_down_top <- g_down_all[g_down_all$ONTO == "MF", ] |>
  dplyr::arrange(dplyr::desc(GeneRatio_num)) |>
  dplyr::slice_head(n = 5)

g_down_top <- dplyr::bind_rows(g_bp_down_top, g_cc_down_top, g_mf_down_top)
g_down_top$Description <- stringr::str_wrap(g_down_top$Description, width = 28)

p_go_down <- ggplot(
  g_down_top,
  aes(x = Count,
      y = reorder_within(Description, GeneRatio_num, ONTO),
      fill = p.adjust)
) +
  geom_col() +
  scale_y_reordered() +
  facet_grid(ONTO ~ ., scales = "free_y", switch = "y") +
  scale_fill_gradient(
    name = "p.adjust",
    low  = "#cfe8ff",
    high = "#084a91",
    trans = "reverse"
  ) +
  labs(x = "Count", y = NULL, title = "GO enrichment (Down-regulated DEGs)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    strip.background = element_rect(fill = "grey70", color = NA),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.title.x    = element_text(face = "bold"),
    axis.text.y     = element_text(face = "bold")
  )

p_go_down
ggsave("GO_down.png", p_go_down, width = 7, height = 8, dpi = 300)

## =========================================================
## KEGG enrichment – Up- and Down-regulated
## =========================================================

# Up-regulated DEGs → ENTREZ
map_up <- bitr(
  up_ids,
  fromType = "TAIR",
  toType   = "ENTREZID",
  OrgDb    = org.At.tair.db
)
entrez_up <- unique(na.omit(map_up$ENTREZID))

ekegg_up <- enrichKEGG(
  gene          = entrez_up,
  organism      = "ath",
  keyType       = "ncbi-geneid",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

n_up_kegg <- nrow(as.data.frame(ekegg_up))

p_kegg_up <- enrichplot::dotplot(ekegg_up, showCategory = n_up_kegg) +
  scale_color_continuous(name = "p.adjust", guide = "legend") +
  labs(x = "GeneRatio", y = NULL, title = "KEGG pathways (Up-regulated DEGs)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text.y  = element_text(face = "bold")
  )

p_kegg_up
ggsave("KEGG_up.png", p_kegg_up, width = 10, height = 5, dpi = 300)

# Down-regulated DEGs → ENTREZ
map_down <- bitr(
  down_ids,
  fromType = "TAIR",
  toType   = "ENTREZID",
  OrgDb    = org.At.tair.db
)
entrez_down <- unique(na.omit(map_down$ENTREZID))

ekegg_down <- enrichKEGG(
  gene          = entrez_down,
  organism      = "ath",
  keyType       = "ncbi-geneid",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

n_down_kegg <- nrow(as.data.frame(ekegg_down))

p_kegg_down <- enrichplot::dotplot(ekegg_down, showCategory = n_down_kegg) +
  scale_color_continuous(name = "p.adjust", guide = "legend") +
  labs(x = "GeneRatio", y = NULL, title = "KEGG pathways (Down-regulated DEGs)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text.y  = element_text(face = "bold")
  )

p_kegg_down
ggsave("KEGG_down.png", p_kegg_down, width = 10, height = 5, dpi = 300)

## ---------- Console summary ----------

cat(sprintf("GO (Up):   BP=%d, CC=%d, MF=%d\n",
            nrow(as.data.frame(ego_bp_up)),
            nrow(as.data.frame(ego_cc_up)),
            nrow(as.data.frame(ego_mf_up))))
cat(sprintf("GO (Down): BP=%d, CC=%d, MF=%d\n",
            nrow(as.data.frame(ego_bp_down)),
            nrow(as.data.frame(ego_cc_down)),
            nrow(as.data.frame(ego_mf_down))))
cat(sprintf("KEGG (Up) terms   = %d\n",   n_up_kegg))
cat(sprintf("KEGG (Down) terms = %d\n", n_down_kegg))