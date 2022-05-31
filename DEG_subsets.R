t3v0_up <- subset(res$t3v0,
                  res$t3v0$padj <= 0.05 & res$t3v0$log2FoldChange >= 1)
nrow(t3v0_up)

t3v0_down <- subset(res$t3v0,
                  res$t3v0$padj <= 0.05 & res$t3v0$log2FoldChange <= -1)
nrow(t3v0_down)
