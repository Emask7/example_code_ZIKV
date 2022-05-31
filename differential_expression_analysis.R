# Run differential expression analysis using DESeq() ---------------------------
  dds <- DESeq(dds)
  resultsNames(dds)

  # Run DESeq() again to analyze day 15 vs 0 -----------------------------------
    dds_15v3 <- DESeq(dds_15v3)
    resultsNames(dds_15v3)

# Obtain results for each comparison and organize them into a list -------------
    # Optional: Use lfcShrink() to calculate values using Bayesian shrinkage
    # estimators with apeglm

  res <- list(
    t3v0 = results(dds, name="timepoint_T3_vs_T0"),
    t15v3 = results(dds_15v3, name="timepoint_T15_vs_T3"),
    t15v0 = results(dds, name="timepoint_T15_vs_T0"),
    t3v0_shrink = lfcShrink(dds, coef="timepoint_T3_vs_T0", type="apeglm"),
    t15v3_shrink = lfcShrink(dds_15v3, coef="timepoint_T15_vs_T3", type="apeglm"),
    t15v0_shrink = lfcShrink(dds, coef="timepoint_T15_vs_T0", type="apeglm")
  )

# Create a matrix containing normalized gene count values ----------------------
  norm_counts <- counts(dds, normalized = TRUE)

# Optional: generate normalized counts matrix with Ensembl IDs for row names ---
  # (useful for downstream analyses, e.g. some gene ontology analyses)

  norm_counts_ensembl <- norm_counts
  rownames(norm_counts_ensembl) <- mapIds(org.Hs.eg.db,
                                          keys = rownames(norm_counts_ensembl),
                                          keytype = "SYMBOL",
                                          column = c("ENSEMBL"))



# Create Excel workbooks and save all results ----------------------------------
  # Normalized counts results --------------------------------------------------
    wb <- createWorkbook("ZIKV DESeq2 Normalized Counts.xlsx")

    addWorksheet(wb, "Gene_Names")
    writeData(wb, "Gene_Names", norm_counts, rowNames = TRUE)

    addWorksheet(wb, "Ensembl_IDs")
    writeData(wb, "Ensembl_IDs", norm_counts_ensembl, rowNames = TRUE)

    saveWorkbook(wb, "ZIKV DESeq2 Normalized Counts.xlsx", overwrite = TRUE)

  # Differential expression analysis results -----------------------------------
    # Recommended that you use a function for this because you have to convert
    # the results to dataframes for the row names to format correctly

    write_to_excel <- function(res, workbook, sheetname) {
      res <- as.data.frame(res)
      addWorksheet(workbook, sheetname)
      writeData(workbook, sheetname, res, rowNames = TRUE)
    }

    wb2 <- createWorkbook("ZIKV DESeq2 Results.xlsx")

    write_to_excel(res$t3v0, wb2, "T3v0")
    write_to_excel(res$t15v3, wb2, "T15v3")
    write_to_excel(res$t15v0, wb2, "T15v0")

    write_to_excel(res$t3v0_shrink, wb2, "T3v0 + lfcShrink")
    write_to_excel(res$t15v3_shrink, wb2, "T15v3 + lfcShrink")
    write_to_excel(res$t15v0_shrink, wb2, "T15v0 + lfcShrink")

    saveWorkbook(wb2, "ZIKV DESeq2 Results.xlsx", overwrite = TRUE)

