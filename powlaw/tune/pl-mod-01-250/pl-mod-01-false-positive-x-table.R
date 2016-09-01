

base_dir <- "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250"
methods <- c("spacemap", "space", "scggm")
fdr_by_xtype_files <- list.files(path = file.path(base_dir, methods), pattern = "fdr_by_xtype.rds", full.names = TRUE)
fdr_by_xtype_files <- rev(fdr_by_xtype_files)
fp_by_xtype <- lapply(fdr_by_xtype_files, readRDS)
names(fp_by_xtype) <- c("Spacemap", "SPACE", "sCGGM")
mean_fp_by_xtype <- t(sapply(fp_by_xtype, function(x) colMeans(x/rowSums(x))))
sd_fp_by_xtype <- t(sapply(fp_by_xtype, function(x) apply(x/rowSums(x),2,sd)))       
res <- matrix(paste0(signif(mean_fp_by_xtype, 4), " (", signif(sd_fp_by_xtype, 2), ")"),3,3)

colnames(res) <- c("Xhubs", "X Confounders", "X Background")
rownames(res) <- rownames(mean_fp_by_xtype)
library(xtable)
print(xtable(res), include.rownames = T)

