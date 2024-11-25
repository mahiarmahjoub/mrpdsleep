## ========================================================================= ##
## =========================== Collider bias =============================== ##

dir <- c("/Users/mahiarmahjoub/OneDrive/Research/6 PD Sleep/GWAS/")
setwd(dir)

require(SlopeHunter)
library(gwasvcf)
library(gwasglue)
library(TwoSampleMR)
library(MRInstruments)
library("ieugwasr")
library(ggpubr)

## -- import result summary files ------------------------------------------ ##
# find all filenames of RDS files 
filenames <- list.files(path = ".", pattern = "res_")
filenames <- subset(filenames, !grepl("colbias",filenames))
filenames <- subset(filenames, grepl("rs",filenames))
filenames <- subset(filenames, !grepl("_PD",filenames))

# import RDS
res_all <- list()
for (i in 1:length(filenames)){
  res_all[[i]] <- readRDS(filenames[i])
}
#res_all <- sapply(filenames, function(x){readRDS(x)})

# change names
names(res_all) <- filenames
names(res_all) <- sub(".RDS","",names(res_all))
names(res_all) <- sub("res_","",names(res_all))
res_all <- lapply(res_all, function(x){
  if(x$mr$nsnp[1]==2){
    x$mr <- x$mr[nrow(x$mr):1,]; x
  } else {x}
}
)


## -- table of PD and progression markers ---------------------------------- ##
i<-5
# initialise dataframe 
colbias_table <- res_all[[i]]$iv[,
                                 c("SNP", "beta.outcome","se.outcome","pval.outcome")]
colbias_table$analysis <- names(res_all)[i]
colbias_table <- colbias_table[!duplicated(colbias_table$SNP),]

# loop to include beta_exp & se_exp for all IVs used in dataframe
# for (i in 2:length(res_all)){
#   colbias_temp <- res_all[[i]]$iv[,c("SNP","beta.outcome","se.outcome","pval.outcome")]
#   colbias_temp$analysis <- names(res_all)[i]
#   colbias_temp <- colbias_temp[!duplicated(colbias_temp$SNP),]
#   colbias_table <- rbind(colbias_table, colbias_temp)
# }
# 
# # add columns for PD estimates 
# colbias_table <- cbind(colbias_table, 
#                        data.frame(matrix(NA, nrow(colbias_table), ncol = 3, 
#                                          dimnames = list(NULL,c("beta.incid",
#                                                                 "se.incid",
#                                                                 "pval.incid"))
#                        )))



# find PD / incidence association of ivs 
# n<-1
# for (i in unique(colbias_table$analysis)){
#   temp <- colbias_table[colbias_table$analysis==i,]
   incid <- ieugwasr::associations(variants = colbias_table$SNP, id = "ieu-b-7")
   incid <- incid[!duplicated(incid$rsid),]
   colbias_table$beta.incid <- incid$beta
   colbias_table$se.incid <- incid$se
   colbias_table$pval.incid <- incid$p
#   
#   n <- n+nrow(incid)-1
# }


# SlopeHunter
# colbias_res <- list()
# for (i in unique(colbias_table$analysis)){

  colbias_res <- hunt(
      dat = colbias_table, #colbias_table[i==colbias_table$analysis,],
      snp_col = "SNP",
      xbeta_col = "beta.incid",
      xse_col = "se.incid",
      xp_col = "pval.incid",
      ybeta_col = "beta.outcome",
      yse_col = "se.outcome",
      yp_col = "pval.outcome",
      xp_thresh = 1, #0.001,
      init_pi = 0.6,
      init_sigmaIP = 1e-05,
      Bootstrapping = TRUE,
      M = 50,
      seed = 777,
      Plot = TRUE,
      show_adjustments = TRUE
    )
#}

# save list and plot
ggsave(paste0("plot_colbias_",colbias_table$analysis[1],".pdf"))
saveRDS(colbias_res, 
        file=paste0("res_colbias_",colbias_table$analysis[1],".RDS"))
