library(DESeq2)
library(ggplot2)
library(pheatmap)
dir.create("results", showWarnings = FALSE)

counts <- read.table("count/counts.txt", header=TRUE, row.names=1)
coldata <- data.frame(group=c("ctrl","ctrl","treat","treat"))
rownames(coldata) <- colnames(counts)

dd <- DESeqDataSetFromMatrix(counts, coldata, ~group)
dds <- DESeq(dd)
res <- results(dds)

# PCA
vsd <- vst(dds)
pca <- plotPCA(vsd, intgroup="group")
ggsave("results/pca.png", pca)

# 火山图
res_df <- as.data.frame(res)
res_df$sig <- ifelse(res_df$padj<0.05&abs(res_df$log2FoldChange)>1,"sig","ns")
vol <- ggplot(res_df,aes(log2FoldChange,-log10(padj),color=sig))+geom_point()
ggsave("results/volcano.png", vol)

# 热图
pheatmap(assay(vsd)[1:50,], filename="results/heatmap.png")
cat("分析+出图完成\n")