## ========================================================================= ##
## ======================== Mendelian Randomisation ======================== ##
## ===================== collate all sleep vs PD results =================== ##

library(phenoscanner)
library(MRPRESSO)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(magrittr)
library(tidyverse)
#.libPaths(c(.libPaths(), "~/Documents/myRlib/"))

dir <- c("D://OneDrive - The University of Sydney (Students)/Research/6 PD Sleep/GWAS/")
#dir <- c("/Users/macbook/OneDrive - The University of Sydney (Students)/Research/6 PD Sleep/GWAS/")
setwd(dir)

## -- import result summary files ------------------------------------------ ##
# find all filenames of RDS files 
filenames <- list.files(path = ".", pattern = "res_")
# import RDS
res_all <- sapply(filenames, function(x){readRDS(x)})


## -- BETA ----------------------------------------------------------------- ##
methods <- c("IVW","Egger","WM","MRPc","CM","lasso")

# initiate data.frame - starting with IVW
mr_beta_colnames <- c("exp","out","beta","beta_lower",
                      "beta_upper","pval","method")
mr_beta <- as.data.frame(matrix(NA, nrow = length(filenames), 
                                ncol = length(mr_beta_colnames), 
                                dimnames = list(NULL,mr_beta_colnames)))
mr_beta$exp <- sapply(filenames, 
                      function(x){strsplit(x, split = "_")[[1]][2]})
mr_beta$out <- sapply(filenames, 
                      function(x){sub(".RDS","",
                                      strsplit(x, split = "_")[[1]][3])}, 
                      USE.NAMES = F)
mr_beta$method <- rep(methods[1],length(filenames))
mr_beta$beta <- unlist(lapply(res_all, function(x){x$mr$b[1]}))
mr_beta$beta_lower <- unlist(lapply(res_all, function(x){x$mr$se[1]}))
mr_beta$beta_upper <- unlist(lapply(res_all, function(x){x$mr$se[1]}))
mr_beta$pval <- unlist(lapply(res_all, function(x){x$mr$pval[1]}))


# loop through each method 
for (i in 2:length(methods)){
  temp_beta <- as.data.frame(matrix(NA, nrow = length(filenames), 
                                  ncol = length(mr_beta_colnames), 
                                  dimnames = list(NULL,mr_beta_colnames)))
  temp_beta$exp <- sapply(filenames, 
                        function(x){strsplit(x, split = "_")[[1]][2]})
  temp_beta$out <- sapply(filenames, 
                        function(x){sub(".RDS","",
                                        strsplit(x, split = "_")[[1]][3])}, 
                        USE.NAMES = F)
  temp_beta$method <- rep(methods[i],length(filenames))
  temp_beta$beta <- unlist(lapply(res_all, function(x){x$mr$b[i]}))
  temp_beta$beta_lower <- unlist(lapply(res_all, function(x){x$mr$se[i]}))
  temp_beta$beta_upper <- unlist(lapply(res_all, function(x){x$mr$se[i]}))
  temp_beta$pval <- unlist(lapply(res_all, function(x){x$mr$pval[i]}))
  
  # MR-PRESSO
  if(methods[i]=="MRP"){
    temp_beta <- as.data.frame(matrix(NA, nrow = length(filenames), 
                                      ncol = length(mr_beta_colnames), 
                                      dimnames = list(NULL,mr_beta_colnames)))
    temp_beta$exp <- sapply(filenames, 
                            function(x){strsplit(x, split = "_")[[1]][2]})
    temp_beta$out <- sapply(filenames, 
                            function(x){sub(".RDS","",
                                            strsplit(x, split = "_")[[1]][3])}, 
                            USE.NAMES = F)
    temp_beta$method <- rep(methods[i],length(filenames))
    temp_beta$beta <- unlist(lapply(res_all, function(x){x$mrp$`Causal Estimate`[1]}))
    temp_beta$beta_lower <- unlist(lapply(res_all, function(x){x$mrp$Sd[1]}))
    temp_beta$beta_upper <- unlist(lapply(res_all, function(x){x$mrp$Sd[1]}))
    temp_beta$pval <- unlist(lapply(res_all, function(x){x$mrp$`P-value`[1]}))
  }
  
  # MR-PRESSO-corrected
  if(methods[i]=="MRPc"){
    temp_beta <- as.data.frame(matrix(NA, nrow = length(filenames), 
                                      ncol = length(mr_beta_colnames), 
                                      dimnames = list(NULL,mr_beta_colnames)))
    temp_beta$exp <- sapply(filenames, 
                            function(x){strsplit(x, split = "_")[[1]][2]})
    temp_beta$out <- sapply(filenames, 
                            function(x){sub(".RDS","",
                                            strsplit(x, split = "_")[[1]][3])}, 
                            USE.NAMES = F)
    temp_beta$method <- rep(methods[i],length(filenames))
    temp_beta$beta <- unlist(lapply(res_all, function(x){x$mrp$`Causal Estimate`[2]}))
    temp_beta$beta_lower <- unlist(lapply(res_all, function(x){x$mrp$Sd[2]}))
    temp_beta$beta_upper <- unlist(lapply(res_all, function(x){x$mrp$Sd[2]}))
    temp_beta$pval <- unlist(lapply(res_all, function(x){x$mrp$`P-value`[2]}))
  }
  
  # collate
  mr_beta <- rbind(mr_beta, temp_beta)
  
}


# turn se into CI 
mr_beta$beta_lower <- mr_beta$beta - 1.96*mr_beta$beta_lower
mr_beta$beta_upper <- mr_beta$beta + 1.96*mr_beta$beta_upper

# turn pval into boolean 
mr_beta$pval <- mr_beta$pval < 0.05

# create a combined out_method variable (for viz purposes)
mr_beta$out_method <- apply(mr_beta[,c("out","method")], 1, 
                            function(x){paste(x,collapse = ".")})



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
## -- beta - MR - version 1 - separate files based on 'EXPOSURE' 
# no MR-PRESSO
for (i in unique(mr_beta$exp)){
  mr_beta %>%
    filter(method == "IVW" | method == "Egger" | method == "WM") %>%
    filter(exp==i) %>%
    ggplot(aes(y=out_method, x=beta, xmin=beta_lower, xmax=beta_upper,
               col=method, lty=pval,
               label=ifelse(pval,"*",""))) +
    geom_point() + 
    geom_text() + 
    geom_errorbarh(height=.1) + 
    geom_text(vjust=0) + 
    facet_wrap(~exp, ncol=2, scales = "free_x") + 
    geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) + 
    coord_cartesian(xlim =c(-100, 100))
  
  ggsave(paste0("plot_summary_beta_nomrp_",i,".pdf"), width = 10, height = 15)
}

# with MR-PRESSO
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
    facet_wrap(~exp, ncol=2, scales = "free_x") + 
    geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) + 
    coord_cartesian(xlim =c(-100, 100))
  
  ggsave(paste0("plot_summary_beta_",i,".pdf"), width = 10, height = 15)
}


## -- beta - MR - version 2 - separate files based on 'METHOD' 
for (i in methods){
  mr_beta %>%
    filter(method==i) %>%
    ggplot(aes(y=out, x=beta, xmin=beta_lower, xmax=beta_upper,col=pval)) +
    geom_point() + 
    geom_errorbarh(height=.1) + 
    facet_wrap(~exp, ncol=2, scales = "free_x") + 
    geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) + 
    coord_cartesian(xlim =c(-100, 100))
    
  ggsave(paste0("plot_summary_beta_",i,".pdf"), width = 10, height = 13)
}


## -- pleiotropy 
plt %>%
  ggplot(aes(y=out, x=yint, xmin=lower, xmax=upper, col=pval)) +
  geom_point() + 
  geom_errorbarh(height=.1) + 
  facet_wrap(~exp, ncol=2, scales = "free_x") + 
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) + 
  coord_cartesian(xlim =c(-4, 4))

ggsave(paste0("plot_plt_summary.pdf"), width = 10, height = 13)


## -- heterogeneity 
het %>%
  ggplot(aes(y=out, x=Isq, fill=method, label=ifelse(pval==T,"*",""))) +
  geom_col(width=.5, position = "dodge") +
  geom_text(hjust=0) + 
  facet_wrap(~exp, ncol=2, scales = "free_x") +
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) + 
  coord_cartesian(xlim =c(-100, 100))

ggsave(paste0("plot_het_summary.pdf"), width = 10, height = 13)



