## ========================================================================= ##
## ======================== Mendelian Randomisation ======================== ##
## ===================== collate all sleep vs PD results =================== ##

library(phenoscanner)
library(MendelianRandomization)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(magrittr)
library(tidyverse)
library(TwoSampleMR)
library(ggpubr)
library(progress)
library(topGO)
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
#library(gprofiler2)
library(forester)
.libPaths(c(.libPaths(), "~/Documents/myRlib/"))

dir <- c("D://OneDrive/Research/6 PD Sleep/GWAS/")
#dir <- c("/Users/macbook/OneDrive - The University of Sydney (Students)/Research/6 PD Sleep/GWAS/")
setwd(dir)
source("rscript_functions_manual_mranalysis.R")

## -- import result summary files ------------------------------------------ ##
# find all filenames of RDS files 
filenames <- list.files(path = ".", pattern = "res_")
#filenames <- subset(filenames, !grepl("ACC",filenames))
#filenames <- subset(filenames, !grepl("rs",filenames))

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
                      if(x$mr$nsnp[1]==2){x$mr <- x$mr[nrow(x$mr):1,]; x} else {x}
                    }
              )

## -- collate all instruments ---------------------------------------------- ##
#extract IVs for each run
# snps_all <- lapply(res_all, function(x){x$iv$SNP})
# # max IV per run & fill rest with NA
# nIVmax <- max(unlist(lapply(snps_all, length)))
# snps_all <- lapply(snps_all, function(x){c(x, rep(NA,nIVmax-length(x)))})
# # place IV SNPS into df
# snps_all_df <- as.data.frame(snps_all)
# # colnames(snps_all_df) <- gsub(".RDS","",colnames(snps_all_df))
# # colnames(snps_all_df) <- gsub("res_","",colnames(snps_all_df))
# # export csv
# write_csv(snps_all_df, "snps_instruments_all.csv")

methods <- c("IVW","Egger","WM","MRPc","CM","lasso")

## -- BETA (new) ----------------------------------------------------------- ##
# prepare dataframe
mr_beta_colnames <- c("exp","out","rsid_rm","beta","beta_lower",
                      "beta_upper","pval","method","nsnp","nsnp_final")
mr_beta <- as.data.frame(matrix(NA, nrow = length(filenames)*length(methods), 
                                ncol = length(mr_beta_colnames), 
                                dimnames = list(NULL,mr_beta_colnames)))
mr_beta$exp <- rep(sapply(filenames, 
                          function(x){strsplit(x, split = "_")[[1]][2]}, 
                          USE.NAMES = F), 
                   (rep(length(methods), length(filenames)))
)
mr_beta$out <- rep(sapply(filenames, 
                          function(x){sub(".RDS","",
                                          strsplit(x, split = "_")[[1]][3])}, 
                          USE.NAMES = F), 
                   (rep(length(methods), length(filenames)))
)
mr_beta$rsid_rm <- rep(sapply(filenames, 
                          function(x){sub(".RDS","",
                                          strsplit(x, split = "_")[[1]][4])}, 
                          USE.NAMES = F), 
                   (rep(length(methods), length(filenames)))
)
mr_beta$method <- rep(methods,length(filenames))
# create a combined out_method variable (for viz purposes)
mr_beta$out_method <- apply(mr_beta[,c("out","method")], 1, 
                            function(x){paste(x,collapse = ".")})


# loop through res_all
index <- 1; index_het <- 1
pb <- progress::progress_bar$new(total = length(res_all))

for (i in 1:length(res_all)){
  pb$tick()
  
  #IVW
  mr_beta$beta[index] <- res_all[[i]]$mr$b[1]
  mr_beta$beta_lower[index] <- res_all[[i]]$mr$b[1] - 1.96*res_all[[i]]$mr$se[1]
  mr_beta$beta_upper[index] <- res_all[[i]]$mr$b[1] + 1.96*res_all[[i]]$mr$se[1]
  mr_beta$pval[index] <- res_all[[i]]$mr$pval[1] < 0.05
  mr_beta$nsnp[index] <- res_all[[i]]$mr$nsnp[1]
  mr_beta$nsnp_final[index] <- 0
  
  #Egger
  mr_beta$beta[index+1] <- res_all[[i]]$mr$b[2]
  mr_beta$beta_lower[index+1] <- res_all[[i]]$mr$b[2] - 1.96*res_all[[i]]$mr$se[2]
  mr_beta$beta_upper[index+1] <- res_all[[i]]$mr$b[2] + 1.96*res_all[[i]]$mr$se[2]
  mr_beta$pval[index+1] <- res_all[[i]]$mr$pval[2] < 0.05
  mr_beta$nsnp[index+1] <- res_all[[i]]$mr$nsnp[2]
  mr_beta$nsnp_final[index+1] <- 0
  
  #WM
  mr_beta$beta[index+2] <- res_all[[i]]$mr$b[3]
  mr_beta$beta_lower[index+2] <- res_all[[i]]$mr$b[3] - 1.96*res_all[[i]]$mr$se[3]
  mr_beta$beta_upper[index+2] <- res_all[[i]]$mr$b[3] + 1.96*res_all[[i]]$mr$se[3]
  mr_beta$pval[index+2] <- res_all[[i]]$mr$pval[3] < 0.05
  mr_beta$nsnp[index+2] <- res_all[[i]]$mr$nsnp[3]
  mr_beta$nsnp_final[index+2] <- 0
  
  #MRPc
  if (class(res_all[[i]]$mrp)=="list"){
    mr_beta$beta[index+3] <- res_all[[i]]$mrp$`Main MR results`$`Causal Estimate`[2]
    mr_beta$beta_lower[index+3] <- mr_beta$beta[index+3] - 1.96*res_all[[i]]$mrp$`Main MR results`$Sd[2]
    mr_beta$beta_upper[index+3] <- mr_beta$beta[index+3] + 1.96*res_all[[i]]$mrp$`Main MR results`$Sd[2]
    mr_beta$pval[index+3] <- res_all[[i]]$mrp$`Main MR results`$`P-value`[2] < 0.05
    mr_beta$nsnp[index+3] <- res_all[[i]]$mr$nsnp[3]
    mr_beta$nsnp_final[index+3] <- res_all[[i]]$mr$nsnp[3]-length(res_all[[i]]$mrp$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  } else {
    mr_beta$beta[index+3] <- NA
    mr_beta$beta_lower[index+3] <- NA
    mr_beta$beta_upper[index+3] <- NA
    mr_beta$pval[index+3] <- NA
    mr_beta$nsnp[index+3] <- res_all[[i]]$mr$nsnp[1]
    mr_beta$nsnp_final[index+3] <- NA
  }
  
  #CM
  if (class(res_all[[i]]$mrcm)=="MRConMix"){
    mr_beta$beta[index+4] <- res_all[[i]]$mrcm@Estimate
    mr_beta$beta_lower[index+4] <- if(length(res_all[[i]]$mrcm@CILower)<2){
                                            res_all[[i]]$mrcm@CILower
                                      } else if (between(res_all[[i]]$mrcm@Estimate, 
                                                         res_all[[i]]$mrcm@CILower[1],
                                                         res_all[[i]]$mrcm@CIUpper[1])) {
                                            res_all[[i]]$mrcm@CILower[1]
                                        } else {
                                            res_all[[i]]$mrcm@CILower[2]
                                          }
    mr_beta$beta_upper[index+4] <- if(length(res_all[[i]]$mrcm@CIUpper)<2){
                                        res_all[[i]]$mrcm@CIUpper
                                      } else if (between(res_all[[i]]$mrcm@Estimate, 
                                                         res_all[[i]]$mrcm@CILower[1],
                                                         res_all[[i]]$mrcm@CIUpper[1])) {
                                        res_all[[i]]$mrcm@CIUpper[1]
                                      } else {
                                        res_all[[i]]$mrcm@CIUpper[2]
                                      }
    mr_beta$pval[index+4] <- res_all[[i]]$mrcm@Pvalue < 0.05
    mr_beta$nsnp[index+4] <- res_all[[i]]$mr$nsnp[3]
    mr_beta$nsnp_final[index+4] <- length(res_all[[i]]$mrcm@Valid)
  } else {
    mr_beta$beta[index+4] <- NA
    mr_beta$beta_lower[index+4] <- NA
    mr_beta$beta_upper[index+4] <- NA
    mr_beta$pval[index+4] <- NA
    mr_beta$nsnp[index+4] <- res_all[[i]]$mr$nsnp[1]
    mr_beta$nsnp_rm[index+4] <- NA
  }
  
  #lasso
  if (class(res_all[[i]]$mrlasso)=="MRLasso"){
    mr_beta$beta[index+5] <- res_all[[i]]$mrlasso@Estimate
    mr_beta$beta_lower[index+5] <- if(length(res_all[[i]]$mrlasso@CILower)<2){
                          res_all[[i]]$mrlasso@CILower
                        } else if (between(res_all[[i]]$mrlasso@Estimate, 
                                           res_all[[i]]$mrlasso@CILower[1],
                                           res_all[[i]]$mrlasso@CIUpper[1])) {
                          res_all[[i]]$mrlasso@CILower[1]
                        } else {
                          res_all[[i]]$mrlasso@CILower[2]
                        }
    mr_beta$beta_upper[index+5] <- if(length(res_all[[i]]$mrlasso@CIUpper)<2){
                          res_all[[i]]$mrlasso@CIUpper
                        } else if (between(res_all[[i]]$mrlasso@Estimate, 
                                           res_all[[i]]$mrlasso@CILower[1],
                                           res_all[[i]]$mrlasso@CIUpper[1])) {
                          res_all[[i]]$mrlasso@CIUpper[1]
                        } else {
                          res_all[[i]]$mrlasso@CIUpper[2]
                        }
    mr_beta$pval[index+5] <- res_all[[i]]$mrlasso@Pvalue < 0.05
    mr_beta$nsnp[index+5] <- res_all[[i]]$mr$nsnp[3]
    mr_beta$nsnp_final[index+5] <- length(res_all[[i]]$mrlasso@Valid)
  } else {
    mr_beta$beta[index+5] <- NA
    mr_beta$beta_lower[index+5] <- NA
    mr_beta$beta_upper[index+5] <- NA
    mr_beta$pval[index+5] <- NA
    mr_beta$nsnp[index+5] <- res_all[[i]]$mr$nsnp[1]
    mr_beta$nsnp_final[index+5] <- NA
  }
    
  index <- index + length(methods)
  index_het <- index_het + 2
  
}


write.table(mr_beta, "mr_beta_all.txt", quote = F,row.names = F, sep = "\t")


## -- BETA (old) ----------------------------------------------------------- ##
# initiate data.frame - starting with IVW
# mr_beta_colnames <- c("exp","out","beta","beta_lower",
#                       "beta_upper","pval","method","nsnp","nsnp_rm")
# mr_beta <- as.data.frame(matrix(NA, nrow = length(filenames), 
#                                 ncol = length(mr_beta_colnames), 
#                                 dimnames = list(NULL,mr_beta_colnames)))
# mr_beta$exp <- sapply(filenames, 
#                       function(x){strsplit(x, split = "_")[[1]][2]})
# mr_beta$out <- sapply(filenames, 
#                       function(x){sub(".RDS","",
#                                       strsplit(x, split = "_")[[1]][3])}, 
#                       USE.NAMES = F)
# mr_beta$method <- rep(methods[1],length(filenames))
# mr_beta$beta <- unlist(lapply(res_all, function(x){x$mr$b[1]}))
# mr_beta$beta_lower <- unlist(lapply(res_all, function(x){x$mr$se[1]}))
# mr_beta$beta_upper <- unlist(lapply(res_all, function(x){x$mr$se[1]}))
# mr_beta$pval <- unlist(lapply(res_all, function(x){x$mr$pval[1]}))
# mr_beta$nsnp <- unlist(lapply(res_all, function(x){
#   if(class(x$mr)=="data.frame"){
#     x$mr$nsnp[3]
#   } else {
#       NA
#     }
#   })
#   )
# mr_beta$nsnp_rm <- 0
# 
# 
# 
# # loop through each method 
# for (i in 2:length(methods)){
#   
#   # Egger
#   if(methods[i]=="Egger"){
#     temp_beta <- as.data.frame(matrix(NA, nrow = length(filenames), 
#                                       ncol = length(mr_beta_colnames), 
#                                       dimnames = list(NULL,mr_beta_colnames)))
#     temp_beta$exp <- sapply(filenames, 
#                             function(x){strsplit(x, split = "_")[[1]][2]})
#     temp_beta$out <- sapply(filenames, 
#                             function(x){sub(".RDS","",
#                                             strsplit(x, split = "_")[[1]][3])}, 
#                             USE.NAMES = F)
#     temp_beta$method <- rep(methods[i],length(filenames))
#     temp_beta$beta <- unlist(lapply(res_all, function(x){x$mr$b[2]}))
#     temp_beta$beta_lower <- unlist(lapply(res_all, function(x){x$mr$se[2]}))
#     temp_beta$beta_upper <- unlist(lapply(res_all, function(x){x$mr$se[2]}))
#     temp_beta$pval <- unlist(lapply(res_all, function(x){x$mr$pval[2]}))
#     temp_beta$nsnp <- unlist(lapply(res_all, function(x){
#       if(class(x$mr)=="data.frame"){
#         x$mr$nsnp[1]
#       } else {
#         NA
#       }
#     })
#     )
#     temp_beta$nsnp_rm <- 0
#   }
#   
#   ## WM
#   if(methods[i]=="WM"){
#   temp_beta <- as.data.frame(matrix(NA, nrow = length(filenames), 
#                                   ncol = length(mr_beta_colnames), 
#                                   dimnames = list(NULL,mr_beta_colnames)))
#   temp_beta$exp <- sapply(filenames, 
#                         function(x){strsplit(x, split = "_")[[1]][2]})
#   temp_beta$out <- sapply(filenames, 
#                         function(x){sub(".RDS","",
#                                         strsplit(x, split = "_")[[1]][3])}, 
#                         USE.NAMES = F)
#   temp_beta$method <- rep(methods[i],length(filenames))
#   temp_beta$beta <- unlist(lapply(res_all, function(x){x$mr$b[3]}))
#   temp_beta$beta_lower <- unlist(lapply(res_all, function(x){x$mr$se[3]}))
#   temp_beta$beta_upper <- unlist(lapply(res_all, function(x){x$mr$se[3]}))
#   temp_beta$pval <- unlist(lapply(res_all, function(x){x$mr$pval[3]}))
#   temp_beta$nsnp <- unlist(lapply(res_all, function(x){
#     if(class(x$mr)=="data.frame"){
#       x$mr$nsnp[1]
#     } else {
#       NA
#     }
#   })
#   )
#   temp_beta$nsnp_rm <- 0
#   }
#   
#   
#   # MR-PRESSO-corrected
#   if(methods[i]=="MRPc"){
#     temp_beta <- as.data.frame(matrix(NA, nrow = length(filenames), 
#                                       ncol = length(mr_beta_colnames), 
#                                       dimnames = list(NULL,mr_beta_colnames)))
#     temp_beta$exp <- sapply(filenames, 
#                             function(x){strsplit(x, split = "_")[[1]][2]})
#     temp_beta$out <- sapply(filenames, 
#                             function(x){sub(".RDS","",
#                                             strsplit(x, split = "_")[[1]][3])}, 
#                             USE.NAMES = F)
#     temp_beta$method <- rep(methods[i],length(filenames))
#     temp_beta$beta <- unlist(lapply(res_all, 
#                                     function(x){
#                                       if(class(x$mrp)=="logical"){
#                                         NA 
#                                       } else {
#                                         x$mrp$`Main MR results`$`Causal Estimate`[2]
#                                         }
#                                         } 
#                                         ))
#     temp_beta$beta_lower <- unlist(lapply(res_all, 
#                                           function(x){
#                                             if(class(x$mrp)=="logical"){
#                                               NA 
#                                             } else {
#                                             x$mrp$`Main MR results`$Sd[2]
#                                             }}))
#     temp_beta$beta_upper <- unlist(lapply(res_all, 
#                                           function(x){
#                                             if(class(x$mrp)=="logical"){
#                                               NA 
#                                             } else {
#                                               x$mrp$`Main MR results`$Sd[2]}}))
#     temp_beta$pval <- unlist(lapply(res_all, function(x){
#                                             if(class(x$mrp)=="logical"){
#                                               NA 
#                                             } else {
#                                               x$mrp$`Main MR results`$`P-value`[2]
#                                               }}))
#     temp_beta$nsnp_rm <- unlist(lapply(res_all, 
#                                        function(x){
#                                          if(class(x$mrp)=="logical"){
#                                            NA 
#                                          } else {
#                                            length(x$mrp$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)}}))
#     temp_beta$nsnp <- unlist(lapply(res_all, function(x){
#                                                 if(class(x$mr)=="data.frame"){
#                                                   x$mr$nsnp[1]
#                                                 } else {
#                                                   NA
#                                                 }
#                                               })
#                                               )
#                                             }
#   
#   # Contamination mixture 
#   if(methods[i]=="CM"){
#     temp_beta <- as.data.frame(matrix(NA, nrow = length(filenames), 
#                                       ncol = length(mr_beta_colnames), 
#                                       dimnames = list(NULL,mr_beta_colnames)))
#     temp_beta$exp <- sapply(filenames, 
#                             function(x){strsplit(x, split = "_")[[1]][2]})
#     temp_beta$out <- sapply(filenames, 
#                             function(x){sub(".RDS","",
#                                             strsplit(x, split = "_")[[1]][3])}, 
#                             USE.NAMES = F)
#     temp_beta$method <- rep(methods[i],length(filenames))
#     
#     temp_beta$beta <- unlist(lapply(res_all, function(x){
#                                                   if(class(x$mrcm)=="MRConMix"){
#                                                     x$mrcm@Estimate
#                                                     } else {
#                                                       NA
#                                                     }
#                                                 }))
# 
#     temp_beta$beta_lower <- unlist(lapply(res_all, 
#                                           function(x){
#                                             if(class(x$mrcm)=="MRConMix"){
#                                               if(length(x$mrcm@CILower)<2){
#                                                 x$mrcm@CILower
#                                               } else if (between(x$mrcm@Estimate, 
#                                                                  x$mrcm@CILower[1], 
#                                                                  x$mrcm@CIUpper[1])) {
#                                                 x$mrcm@CILower[1]
#                                               } else {
#                                                 x$mrcm@CILower[2]
#                                               }
#                                             } else {
#                                               NA
#                                             }
#                                             }))
#     temp_beta$beta_upper <- unlist(lapply(res_all, 
#                                           function(x){
#                                             if(class(x$mrcm)=="MRConMix"){
#                                               if(length(x$mrcm@CIUpper)<2){
#                                                 x$mrcm@CIUpper
#                                               } else if (between(x$mrcm@Estimate, 
#                                                                  x$mrcm@CILower[1], 
#                                                                  x$mrcm@CIUpper[1])) {
#                                                 x$mrcm@CIUpper[1]
#                                               } else {
#                                                 x$mrcm@CIUpper[2]
#                                               }
#                                             } else {
#                                               NA
#                                             }
#                                           }))
#     
#     temp_beta$pval <- unlist(lapply(res_all, function(x){
#                                                   if(class(x$mrcm)=="MRConMix"){
#                                                     x$mrcm@Pvalue
#                                                     } else {
#                                                       NA
#                                                     }}))
#     temp_beta$nsnp_rm <- unlist(lapply(res_all, 
#                                        function(x){
#                                          if(class(x$mrcm)=="MRConMix"){
#                                            (x$mrcm@SNPs-length(x$mrcm@Valid))
#                                            } else {
#                                              NA
#                                            }}))
#     temp_beta$nsnp <- unlist(lapply(res_all, function(x){
#       if(class(x$mr)=="data.frame"){
#         x$mr$nsnp[1]
#       } else {
#         NA
#       }
#     })
#     )
#   }
#   
#   # MR-Lasso
#   if(methods[i]=="lasso"){
#     temp_beta <- as.data.frame(matrix(NA, nrow = length(filenames), 
#                                       ncol = length(mr_beta_colnames), 
#                                       dimnames = list(NULL,mr_beta_colnames)))
#     temp_beta$exp <- sapply(filenames, 
#                             function(x){strsplit(x, split = "_")[[1]][2]})
#     temp_beta$out <- sapply(filenames, 
#                             function(x){sub(".RDS","",
#                                             strsplit(x, split = "_")[[1]][3])}, 
#                             USE.NAMES = F)
#     temp_beta$method <- rep(methods[i],length(filenames))
#     
#     temp_beta$beta <- unlist(lapply(res_all, function(x){
#       if(class(x$mrlasso)!="MRLasso"){NA} else {x$mrlasso@Estimate}}))
#     temp_beta$beta_lower <- unlist(lapply(res_all, function(x){
#       if(class(x$mrlasso)!="MRLasso"){NA} else {x$mrlasso@CILower}}))
#     temp_beta$beta_upper <- unlist(lapply(res_all, function(x){
#       if(class(x$mrlasso)!="MRLasso"){NA} else {x$mrlasso@CIUpper}}))
#     temp_beta$pval <- unlist(lapply(res_all, function(x){
#       if(class(x$mrlasso)!="MRLasso"){NA} else {x$mrlasso@Pvalue}}))
#     temp_beta$nsnp_rm <- unlist(lapply(res_all, 
#                                        function(x){
#                                          if(class(x$mrlasso)!="MRLasso"){
#                                            NA
#                                            } else {
#                                            (x$mrlasso@SNPs-x$mrlasso@Valid)
#                                              }}))
#     #temp_beta$nsnp <- unlist(lapply(res_all, function(x){x$mrcm@SNPs}))
#     temp_beta$nsnp <- unlist(lapply(res_all, function(x){
#       if(class(x$mr)=="data.frame"){
#         x$mr$nsnp[1]
#       } else {
#         NA
#       }
#     })
#     )
#   }
#   
#   # collate
#   mr_beta <- rbind(mr_beta, temp_beta)
#   
# }
# 
# 
# # turn se into CI 
# mr_beta$beta_lower[1:(4*length(filenames))] <- mr_beta$beta[1:(4*length(filenames))] - 
#   1.96*mr_beta$beta_lower[1:(4*length(filenames))]
#                                                     
# mr_beta$beta_upper[1:(4*length(filenames))] <- mr_beta$beta[1:(4*length(filenames))] + 
#   1.96*mr_beta$beta_upper[1:(4*length(filenames))]
# 
# # turn pval into boolean 
# mr_beta$pval <- mr_beta$pval < 0.05
# 
# # create a combined out_method variable (for viz purposes)
# mr_beta$out_method <- apply(mr_beta[,c("out","method")], 1, 
#                             function(x){paste(x,collapse = ".")})



## -- PLEIOTROPY ----------------------------------------------------------- ##
plt_colnames <- c("exp","out","yint","lower","upper","pval")
plt <- as.data.frame(matrix(NA, nrow = length(filenames), 
                                ncol = length(plt_colnames), 
                                dimnames = list(NULL,plt_colnames)))
plt$exp <- sapply(filenames, 
                      function(x){strsplit(x, split = "_")[[1]][2]})
plt$out <- sapply(filenames, 
                      function(x){sub(".RDS","",
                                      strsplit(x, split = "_")[[1]][3])}, 
                      USE.NAMES = F)
plt$yint <- unlist(lapply(res_all, function(x){x$pleiotropy$egger_intercept}))
plt$lower <- unlist(lapply(res_all, function(x){x$pleiotropy$se}))
plt$upper <- unlist(lapply(res_all, function(x){x$pleiotropy$se}))
plt$pval <- unlist(lapply(res_all, function(x){x$pleiotropy$pval}))

# convert se to CI 
plt$lower <- plt$yint - abs(1.96*plt$lower)
plt$upper <- plt$yint + abs(1.96*plt$upper)

# turn pval into boolean 
plt$pval <- plt$pval < 0.05



## -- HETEROGENEITY -------------------------------------------------------- ##
# IVW
het_colnames <- c("exp","out","method","Q",
                      "Qdf","pval","Isq")
het <- as.data.frame(matrix(NA, nrow = length(filenames), 
                                ncol = length(het_colnames), 
                                dimnames = list(NULL,het_colnames)))
het$exp <- sapply(filenames, 
                      function(x){strsplit(x, split = "_")[[1]][2]})
het$out <- sapply(filenames, 
                      function(x){sub(".RDS","",
                                      strsplit(x, split = "_")[[1]][3])}, 
                      USE.NAMES = F)
het$method <- rep(methods[1],length(filenames))
het$Q <- unlist(lapply(res_all, function(x){x$het$Q[2]}))
het$Qdf <- unlist(lapply(res_all, function(x){x$het$Q_df[2]}))
het$Isq <- unlist(lapply(res_all, function(x){x$isq[2]}))
het$pval <- unlist(lapply(res_all, function(x){x$het$Q_pval[2]}))

# Egger
temp_het <- as.data.frame(matrix(NA, nrow = length(filenames), 
                            ncol = length(het_colnames), 
                            dimnames = list(NULL,het_colnames)))
temp_het$exp <- sapply(filenames, 
                  function(x){strsplit(x, split = "_")[[1]][2]})
temp_het$out <- sapply(filenames, 
                  function(x){sub(".RDS","",
                                  strsplit(x, split = "_")[[1]][3])}, 
                  USE.NAMES = F)
temp_het$method <- rep(methods[2],length(filenames))
temp_het$Q <- unlist(lapply(res_all, function(x){x$het$Q[1]}))
temp_het$Qdf <- unlist(lapply(res_all, function(x){x$het$Q_df[1]}))
temp_het$Isq <- unlist(lapply(res_all, function(x){x$isq[1]}))
temp_het$pval <- unlist(lapply(res_all, function(x){x$het$Q_pval[1]}))

# combine 
het <- rbind(het, temp_het)

# turn pval into boolean 
het$pval <- het$pval < 0.05



## -- forest plot ---------------------------------------------------------- ##
for (i in unique(mr_beta$exp)){
  mr_beta %>%
    filter(exp==i) %>%
    ggplot(aes(y=out_method, x=beta, xmin=beta_lower, xmax=beta_upper,
               col=method, lty=pval,
               label=ifelse(pval,"*",""))) +
    geom_point() + 
    geom_text() + 
    geom_errorbarh(height=.1) + 
    geom_text(vjust=0) + 
    facet_wrap(~exp, ncol=1, scales = "free_x") + 
    geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) + 
    coord_cartesian(xlim =c(min(mr_beta$beta_lower), max(mr_beta$beta_upper)))
  
  ggsave(paste0("plot_summary_beta_",i,".pdf"), width = 10, height = 15)
}


# UPDRS, MOCA, MMSE only 
# for (i in unique(mr_beta$exp)){
#   mr_beta %>%
#     filter(exp==i) %>%
#     filter(grepl("UPDRS",out) | grepl("MOCA",out) | 
#              grepl("MMSE",out) | grepl("PD",out)) %>%
#     ggplot(aes(y=out_method, x=beta, xmin=beta_lower, xmax=beta_upper,
#                col=method, lty=pval,
#                label=ifelse(pval,"*",""))) +
#     geom_point() + 
#     geom_text() + 
#     geom_errorbarh(height=.1) + 
#     geom_text(vjust=0) + 
#     facet_wrap(~exp, ncol=1, scales = "free_x") + 
#     geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)  
#   coord_cartesian(xlim =c(min(mr_beta$beta_lower, na.rm = T), max(mr_beta$beta_upper, na.rm = T)))
#   
#   ggsave(paste0("plot_summary_beta_",i,"_selected.pdf"), 
#          width = 10, height = 10)
# }

# forest plot - MMSE, MOCA and UPDRS3 only 
mr_beta %>%
  # filter(out=="contMOCA" | out=="contMMSE" | out=="contUPDRS3" | out=="contUPDRS") %>%
   filter(exp=="sleepduration" | exp=="insomnia") %>%
  # ggplot(aes(y=out_method, x=beta, xmin=beta_lower, xmax=beta_upper,
  #            col=method, lty=pval,
  #            label=ifelse(pval,"*",""))) +
  ggplot(aes(y=out_method, x=beta, xmin=beta_lower, xmax=beta_upper,
             col=method, lty=pval,
             label=ifelse(pval,"*",""))) +
  geom_point() + 
  geom_text() + 
  # facet_grid(rows = vars(out)) + 
  geom_errorbarh(height=.1) + 
  geom_text(vjust=0) + 
  #facet_grid(rows = var(out)) + 
  facet_wrap(~exp, ncol=2, scales = "free_x") + 
  #facet_wrap(~out, ncol=1, scales = "free_x") + 
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)  
  #scale_x_continuous(limits = c(-75,75))
  #coord_cartesian(xlim =c(min(mr_beta$beta_lower), max(mr_beta$beta_upper)))

ggsave(paste0("plot_summary_beta_ALL.pdf"), 
       width = 10, height = 8)


## -- pleiotropy 
plt %>%
  filter(!grepl("ACC", exp)) %>%
  ggplot(aes(y=out, x=yint, xmin=lower, xmax=upper, col=pval)) +
  geom_point() + 
  geom_errorbarh(height=.1) + 
  facet_wrap(~exp, ncol=1, scales = "free") + 
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)  
  #coord_cartesian(xlim =c(-4, 4))

ggsave(paste0("plot_plt_summary.pdf")) #, width = 20, height = 10)

## -- pleiotropy - ACC only
# plt %>%
#   filter(grepl("ACC", exp)) %>%
#   ggplot(aes(y=out, x=yint, xmin=lower, xmax=upper, col=pval)) +
#   geom_point() + 
#   geom_errorbarh(height=.1) + 
#   facet_wrap(~exp, ncol=2, scales = "free_x") + 
#   geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) + 
#   coord_cartesian(xlim =c(-4, 4))
# 
# ggsave(paste0("plot_plt_summary_accel.pdf"), width = 10, height = 20)

## -- heterogeneity - Isq
het$Isq_mod <- ifelse(het$Isq<0, 0,het$Isq )
het %>%
  #filter(!grepl("ACC", exp)) %>%
  ggplot(aes(y=out, x=Isq, fill=method, label=ifelse(pval==T,"*",""))) +
  geom_col(width=.5, position = "dodge") +
  geom_text(hjust=0) + 
  facet_wrap(~exp, ncol=1, scales = "free") +
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +  
  coord_cartesian(xlim =c(0, 100))

ggsave(paste0("plot_het_Isq_summary.pdf"))#, width = 10, height = 20)

## -- heterogeneity - Q-score
het %>%
  #filter(!grepl("ACC", exp)) %>%
  ggplot(aes(y=out, x=Q, fill=method, label=ifelse(pval==T,"*",""))) +
  geom_col(width=.5, position = "dodge") +
  geom_text(hjust=0) + 
  facet_wrap(~exp, ncol=1, scales = "free") +
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)  
  #coord_cartesian(xlim =c(0, 100))

ggsave(paste0("plot_het_Q_summary.pdf"))#, width = 10, height = 20)

## -- heterogeneity - ACC only 
# het %>%
#   filter(grepl("ACC", exp)) %>%
#   ggplot(aes(y=out, x=Isq, fill=method, label=ifelse(pval==T,"*",""))) +
#   geom_col(width=.5, position = "dodge") +
#   geom_text(hjust=0) + 
#   facet_wrap(~exp, ncol=2, scales = "free_x") +
#   geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) + 
#   coord_cartesian(xlim =c(-100, 100))
# 
# ggsave(paste0("plot_het_summary_accel.pdf"), width = 10, height = 20)


## -- beta - MR - version 2 - separate files based on 'METHOD' 
# for (i in methods){
#   mr_beta %>%
#     filter(method==i) %>%
#     ggplot(aes(y=out, x=beta, xmin=beta_lower, xmax=beta_upper,col=pval,
#                label=ifelse(pval,"*",""))) +
#     geom_point() + 
#     geom_text() + 
#     geom_errorbarh(height=.1) + 
#     facet_wrap(~exp, ncol=2, scales = "free_x") + 
#     geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) + 
#     coord_cartesian(xlim =c(-100, 100))
#     
#   ggsave(paste0("plot_summary_beta_",i,".pdf"), width = 10, height = 13)
# }


## -- beta - MR - version 1 - separate files based on 'EXPOSURE' 
# no MR-PRESSO
# for (i in unique(mr_beta$exp)){
#   mr_beta %>%
#     filter(method == "IVW" | method == "Egger" | method == "WM") %>%
#     filter(exp==i) %>%
#     ggplot(aes(y=out_method, x=beta, xmin=beta_lower, xmax=beta_upper,
#                col=method, lty=pval,
#                label=ifelse(pval,"*",""))) +
#     geom_point() + 
#     geom_text() + 
#     geom_errorbarh(height=.1) + 
#     geom_text(vjust=0) + 
#     facet_wrap(~exp, ncol=2, scales = "free_x") + 
#     geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) + 
#     coord_cartesian(xlim =c(-100, 100))
#   
#   #ggsave(paste0("plot_summary_beta_nomrp_",i,".pdf"), width = 10, height = 15)
# }


