# Set up a dataframe of raw gene counts and specify factors --------------------
  coldata <- read.xlsx("metadata.xlsx")
  rownames(coldata) <- coldata[, 1]
  coldata <- coldata[, c("timepoint", "animal")]
  coldata$timepoint <- factor(coldata$timepoint)
  coldata$animal <- factor(coldata$animal)
  coldata

  raw_data <- read.xlsx("1543PC_raw_gene_counts_Panubis1.xlsx")
  colnames(raw_data) <- c(
    "Chromosome", "Start", "Stop", "Strand", "Gene_Symbol", rownames(coldata)
  )
  raw_gene_counts <- as.matrix(raw_data[, 6:17])
  rownames(raw_gene_counts) <- raw_data[, 5]
  raw_gene_counts <- subset(
    raw_gene_counts, !is.na(rownames(raw_gene_counts))
  )
  head(raw_gene_counts)



# Create a DESeqDataSet --------------------------------------------------------
  cts <- round(raw_gene_counts)
  head(cts)
  dds_pca <- DESeqDataSetFromMatrix(countData = cts,
                                    colData = coldata,
                                    design = ~ animal + timepoint)

# Perform Variance stabilizing transformation ----------------------------------
  vsd <- vst(dds_pca, blind=FALSE)
  pcaData <- plotPCA(vsd,
                     intgroup=c("timepoint", "animal"),
                     returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  all_pca <- ggplot(pcaData, aes(PC1, PC2, color=timepoint, shape=animal)) +
    geom_point(size=3, stroke=1.05,
               fill=c(NA, "#F8766D", NA, "#F8766D",
                      NA, "#00BFC4", NA, "#00BFC4",
                      NA, "#7CAE00", NA, "#7CAE00")) +
    scale_shape_manual(values=c(1, 15, 2, 25)) +
    scale_color_manual(values=c('#F8766D','#7CAE00', '#00BFC4'))+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    guides(shape = guide_legend(override.aes = list(fill = "black"))) +
    theme(aspect.ratio = 1) +
    coord_fixed()

  all_pca

  ggsave("PCA plot.png", plot = all_pca,
         units="in", width=6, height=4, dpi=300)

  # The results of the PCA plot indicate that animal 32535 is an outlier at T0
  # For downstream analyses, we must remove this dataset

# Removing outlier (32535 - T0) ------------------------------------------------
  cts <- raw_gene_counts[, c(2:12)]
    # excluding outlier only at day 0
  cts <- round(cts)
  head(cts)

  coldata <- coldata[c(2:12), ]
  coldata

# PCA Plots - outlier excluded -----------------------------------------------
  dds_pca_excl <- DESeqDataSetFromMatrix(countData = cts,
                                         colData = coldata,
                                         design = ~ animal + timepoint)

# Variance stabilizing transformation --------------------------------------
  vsd_excl <- vst(dds_pca_excl, blind=FALSE)
  pcaData_excl <- plotPCA(vsd_excl,
                          intgroup=c("timepoint", "animal"),
                          returnData=TRUE)
  percentVar <- round(100 * attr(pcaData_excl, "percentVar"))
  excl_pca <- ggplot(pcaData_excl, aes(PC1, PC2, color=timepoint, shape=animal)) +
    geom_point(size=3, stroke=1.05,
               fill=c("#F8766D", NA, "#F8766D",
                      NA, "#00BFC4", NA, "#00BFC4",
                      NA, "#7CAE00", NA, "#7CAE00")) +
    scale_shape_manual(values=c(1, 15, 2, 25)) +
    scale_color_manual(values=c('#F8766D','#7CAE00', '#00BFC4'))+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    guides(shape = guide_legend(override.aes = list(fill = "black"))) +
    theme(aspect.ratio = 1) +
    coord_fixed()

  excl_pca

  ggsave("PCA plot - T0 3253 excl.png", plot = excl_pca,
         units="in", width=6, height=4, dpi=300)

  # The results of the second PCA plot indicate that animal 32535 is an outlier
  # at the two other time points as well
  # For downstream analyses, we must remove all data from this animal

# Removing outlier (32535 - all timepoints) ------------------------------------
  cts <- raw_gene_counts[, c(2:4, 6:8, 10:12)]
    # excluding outlier at all timepoints
  cts <- round(cts)
  head(cts)

  coldata <- coldata[c(1:3, 5:7, 9:11), ]
  coldata

# create a final DESeqDataSet to use for downstream analyses -------------------
  coldata$timepoint <- factor(coldata$timepoint)
  coldata$animal <- factor(coldata$animal)
     # must reassign factors before setting up DESeqDataSet

  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ animal + timepoint)
