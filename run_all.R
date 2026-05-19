# ==============================================
# RNA-seq 复现一键脚本（修复版，解决vst报错）
# ==============================================

# 强制设置工作目录
setwd("F:/RNA-seq-Reproduction-1")

# 1. 安装并加载依赖包
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cran.r-project.org/")
}
BiocManager::install(c("DESeq2", "ggplot2", "pheatmap"), update = FALSE, ask = FALSE)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# 2. 创建目录
dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("count", recursive = TRUE, showWarnings = FALSE)

# 3. 生成足够大的模拟数据（避免vst报错）
set.seed(123)
# 生成1000个基因的计数矩阵（足够大，满足vst的nsub要求）
counts <- matrix(round(runif(4000, 1, 1000)), ncol = 4)
colnames(counts) <- c("control1", "control2", "treat1", "treat2")
rownames(counts) <- paste0("gene_", 1:nrow(counts))
write.table(counts, "count/counts.txt", sep = "\t", quote = FALSE, row.names = TRUE)
cat("✅ 模拟计数矩阵已生成到 count/counts.txt\n")

# 4. DESeq2差异分析
coldata <- data.frame(
    group = factor(c("control", "control", "treat", "treat")),
    row.names = colnames(counts)
)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "results/de_result.csv", row.names = TRUE)
cat("✅ 差异分析结果已保存到 results/de_result.csv\n")

# 5. PCA图（改用rlog，避免vst报错）
vsd <- rlog(dds) # 用rlog替代vst，更稳定
pca_plot <- plotPCA(vsd, intgroup = "group") + 
    theme_bw(base_size = 14) + 
    ggtitle("PCA of RNA-seq Samples")
ggsave("results/pca.png", pca_plot, width = 6, height = 4, dpi = 300)
cat("✅ PCA图已保存到 results/pca.png\n")

# 6. 火山图
res_df <- as.data.frame(res)
res_df$sig <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
    theme_bw(base_size = 14) +
    labs(title = "Volcano Plot of DE Genes", x = "log2(Fold Change)", y = "-log10(Adjusted p-value)")
ggsave("results/volcano.png", volcano_plot, width = 6, height = 4, dpi = 300)
cat("✅ 火山图已保存到 results/volcano.png\n")

# 7. 热图（取前50个差异最显著的基因）
top_genes <- rownames(res_df)[order(res_df$padj)[1:50]]
pheatmap(assay(vsd)[top_genes, ], 
         filename = "results/heatmap.png", 
         width = 7, height = 6, 
         show_rownames = FALSE, scale = "row")
cat("✅ 热图已保存到 results/heatmap.png\n")

cat("\n🎉 所有步骤完成！请检查 results/ 文件夹\n")