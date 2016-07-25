setwd("~/scratch-data/sim-spacemap/hub11-20/2016/results/n250/")
voted <- readRDS("spacemap/spacemap_top_cv_voting_performance.rds")

files <- c("spacemap/spacemap_top_cv_voting_performance.rds", 
           "space/space_top_cv_voting_performance.rds", 
           "scggm/scggm_top_cv_voting_performance.rds")

library(data.table)
library(foreach)
#Tables
res <- foreach(f = files, .combine = 'rbind') %do% { 
  voted <- readRDS(f)
  xbarv <- voted[,list(mean_mcc = mean(mcc), mean_power = mean(power), mean_fdr = mean(fdr)),by = comparison]
  sdv <- voted[,list(sd_mcc = sd(mcc), sd_power = sd(power), sd_fdr = sd(fdr)),by = comparison]
  power <- c(xbarv$mean_mcc[1], xbarv$mean_power[2:3], xbarv$mean_fdr[2:3])
  sd <- c(sdv$sd_mcc[1], sdv$sd_power[2:3], sdv$sd_fdr[2:3])
  fmt <- paste0(signif(power, 4), " (", signif(sd, 2), ")")
}
res2 <- foreach(f = files) %do% { 
  voted <- readRDS(f)
  xbarv <- voted[,list(mean_mcc = mean(mcc), mean_power = mean(power), mean_fdr = mean(fdr)),by = comparison]
  sdv <- voted[,list(sd_mcc = sd(mcc), sd_power = sd(power), sd_fdr = sd(fdr)),by = comparison]
  list(xbar = xbarv, sd = sdv)
}

colnames(res) <- c("MCC","YY Power", "XY Power", "YY FDR", "XY FDR")
rownames(res) <- c("Spacemap", "SPACE", "sCGGM")
library(xtable)
print(xtable(res), include.rownames = T)

#Comparisons
library(ggplot2)
lvoted <- lapply(files, readRDS)
methods <- c("Spacemap", "SPACE", "sCGGM")
lvoted <- lapply(1:3, function(i) { lvoted[[i]]$method <- methods[i]; lvoted[[i]] })
voted <- rbindlist(lvoted)
setkey(voted, comparison)

library(reshape2)
gvoted <- melt(as.data.frame(voted)[,c(1:4,8)])
levels(gvoted$variable)[levels(gvoted$variable)=="mcc"] <- "MCC"
levels(gvoted$variable)[levels(gvoted$variable)=="power"] <- "Power"
levels(gvoted$variable)[levels(gvoted$variable)=="fdr"] <- "FDR"
gvoted$variable <-  factor(x= gvoted$variable, levels = c("MCC", "Power", "FDR"))
fig0 <- ggplot(data = gvoted, aes(x = method, y = value)) + 
  geom_boxplot() + facet_grid(variable ~ comparison, scale = "free") + theme_bw()
ggsave(filename = "~/Dropbox/Chris_Conley/ms-spacemap/figures/hub11-20-performance-summary.png", 
       plot = fig0)

########Potential
all_files <- c("spacemap/spacemap_all_tuning_cv_vote_results.rds", 
           "space/space_all_tuning_cv_vote_results.rds", 
           "scggm/scggm_all_tuning_cv_vote_results.rds")
lall_voted <- lapply(all_files, readRDS)

#amv <- lall_voted[[1]]
#damv <- amv[1:(nrow(amv)/100),]

#visualise potential vs. selected
powerVfdr <- function(res, metaInfo = "", dotsize = 0.5, alpha = 1, 
                      powerZoom = c(0.65,1), fdrZoom = c(0, 0.3), pal = "Paired") {
  library(ggplot2)
  ggplot(aes(x=fdr, y = power, colour = method), data = res) + 
    facet_grid(. ~ comparison) +
    xlim(fdrZoom) + ylim(powerZoom) + 
    #scale_y_continuous(breaks = seq(powerZoom[1], powerZoom[2], 0.05)) + 
    #scale_x_continuous(breaks = seq(fdrZoom[1], fdrZoom[2], 0.1))  + 
    geom_point(size = dotsize, alpha = alpha) + 
    scale_color_brewer(palette = pal) + 
    #geom_smooth(size = 2, formula = y ~ s(x, bs = "cs"), method = "gam", na.rm = TRUE)  + 
    ggtitle(metaInfo) +
    theme(axis.text=element_text(size=16)) + 
    theme_bw()
}


#relabel methods in color-convenient order
mv <- as.data.frame(lvoted[[1]])
spv <- as.data.frame(lvoted[[2]])
scv <- as.data.frame(lvoted[[3]])

ordered_levels <- c("sCGGM.CV.vote", "sCGGM.All","Spacemap.All", "Spacemap.CV.vote")
mv$method <- factor(x = "Spacemap.CV.vote", levels = ordered_levels)
scv$method <- factor(x = "sCGGM.CV.vote", levels = ordered_levels)
lall_voted <- lapply(lall_voted, as.data.frame)
lall_voted[[1]]$method <- factor(x = "Spacemap.All", levels = ordered_levels)
lall_voted[[3]]$method <- factor(x = "sCGGM.All", levels = ordered_levels)
#Spacemap vs. sCGGM
cmp <- rbindlist(lall_voted[c(1,3)])
fig1 <- powerVfdr(res = cmp, powerZoom = c(0.50,1), fdrZoom = c(0.0,0.4), pal = "Paired")  +
  geom_point(data = scv, size = 2) +  
  geom_point(data = mv, size = 2) 
ggsave(filename = "~/Dropbox/Chris_Conley/ms-spacemap/figures/hub11-20-power-vs-fdr-spmapVscggm.png", 
       plot = fig1)



ordered_levels <- c("SPACE.CV.vote", "SPACE.All","Spacemap.All", "Spacemap.CV.vote")
mv$method <- factor(x = "Spacemap.CV.vote", levels = ordered_levels)
spv$method <- factor(x = "SPACE.CV.vote", levels = ordered_levels)
lall_voted <- lapply(lall_voted, as.data.frame)
lall_voted[[1]]$method <- factor(x = "Spacemap.All", levels = ordered_levels)
lall_voted[[2]]$method <- factor(x = "SPACE.All", levels = ordered_levels)
#Spacemap vs. SPACE
cmp <- as.data.frame(rbindlist(lall_voted[c(1,2)]))
fig2 <- powerVfdr(res = cmp, powerZoom = c(0.35,1), fdrZoom = c(0.0,0.4), pal = "Paired")  +
  geom_point(data = spv, size = 2) +  
  geom_point(data = mv, size = 2)  
ggsave(filename = "~/Dropbox/Chris_Conley/ms-spacemap/figures/hub11-20-power-vs-fdr-spmapVspace.png", 
       plot = fig2)
