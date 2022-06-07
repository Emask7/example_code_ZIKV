t3v0_reg_of_viral_rep <- read.xlsx("t3v0_GOBP_NEGATIVE_REGULATION_OF_VIRAL_GENOME_REPLICATION.xlsx")
head(t3v0_reg_of_viral_rep)

t3v0_res_table <- as.data.frame(res$t3v0)
t3v0_res_table[, 7] <- rownames(t3v0_res_table)
colnames(t3v0_res_table) <- c(colnames(t3v0_res_table[, 1:6]), "Gene_symbol")
head(t3v0_res_table)

t3v0_reg_of_viral_rep <- inner_join(t3v0_reg_of_viral_rep, t3v0_res_table)
rownames(t3v0_reg_of_viral_rep) <- t3v0_reg_of_viral_rep$Gene_symbol
t3v0_reg_of_viral_rep <- t3v0_reg_of_viral_rep[, c(3, 7)]
head(t3v0_reg_of_viral_rep)
nrow(t3v0_reg_of_viral_rep)
