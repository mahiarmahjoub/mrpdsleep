#!/usr/bin/env Rscript
## ========================================================================= ##
## ======================== Mendelian Randomisation ======================== ##
## =================================TwoSampleMR============================= ##
## ============================ sleep IEU vs PD IEU  ======================= ##

## all exposure and outcome datasets from IEU 

.libPaths(c(.libPaths(), "~/Documents/myRlib/"))

library(gwasvcf)
library(gwasglue)
library(biomaRt)
library(TwoSampleMR)
library(MendelianRandomization)
library(magrittr)
library(phenoscanner)
library(MRPRESSO)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(progress)
library(ggpubr)

gwasvcf::set_bcftools()
gwasvcf::set_plink()

#dir <- c("D://OneDrive - The University of Sydney (Students)/Research/6 PD Sleep/GWAS/")
dir <- c("/Users/macbook/OneDrive - The University of Sydney (Students)/Research/6 PD Sleep/GWAS/")
setwd(dir)

# filenames
name_exposure <- c("insomnia","sleepduration","chronotype","naps")
name_outcome <- c("PD")
filename_exposure <- c("ukb-b-3957","ukb-b-4424","ukb-b-4956","ukb-b-4616")
filename_outcome <- c("ieu-b-7")


## -- setup parameters ----------------------------------------------------- ##
# LD proxy search 
r2_proxy <- 0.8
kb_proxy <- 10000
# clumping - post-harmonisation
clump_kb <- 10000
clump_r2 <- 0.001
# Phenoscanner SNP search - pvalue cutoff
pheno_pval <- 5e-8
# collate params 
param <- c(r2_proxy=r2_proxy, kb_proxy=kb_proxy, 
           clump_kb=clump_kb, clump_r2=clump_r2, pheno_pval=pheno_pval)


## == LOOP - MR analysis pipeline ========================================== ##
# initialise grid noting errors in loop 
error_grid <- expand.grid(name_outcome, name_exposure) 
error_grid$error <- rep(NA, nrow(error_grid))
index_error <- 1
# initialise progressbar 
pbi <- progress_bar$new(total = length(name_exposure))

for (i in 1:length(name_exposure)){
  pbi$tick()
  skip_to_next <- FALSE
  tryCatch({
## -- exposure ------------------------------------------------------------- ##
  print("Importing exposure ...")
  exposure_dat <- extract_instruments(outcomes = filename_exposure[i], 
                                    p1 = 5e-8, clump = F)
  IVs <- exposure_dat$SNP

pbj <- progress_bar$new(total = length(name_outcome))
for (j in 1:length(name_outcome)){
  pbj$tick()
  print(paste0("* ",name_exposure[i]," vs ",name_outcome[j]," *"))
  skip_to_next <- FALSE
  tryCatch({
## -- outcome -------------------------------------------------------------- ##
  print("Importing outcome & proxy search ...")
  outcome_dat <- extract_outcome_data(
    snps = IVs,
    outcomes = filename_outcome[j], 
    proxies = T, rsq = r2_proxy
  )


## -- harmonise data - exposure and outcome -------------------------------- ##
  print("Harmonisation ...")
    dat <- harmonise_data(
    exposure_dat = exposure_dat, 
    outcome_dat = outcome_dat,
    action = 3
  )

  # remove palindromic snps prior to clumping
  dat <- subset(dat, mr_keep == TRUE)

## -- clump final harmonised data ------------------------------------------ ##
  print("Clumping ...")
  dat_clumped <- clump_data(dat = dat, clump_kb = clump_kb, clump_r2 = clump_r2)
  dat_clumped <- subset(dat_clumped, !duplicated(dat_clumped)) #remove duplicates

  
## -- MR analysis ---------------------------------------------------------- ##
  print("MR analysis ...")

  ## -- IVW, WM, Egger
  res <- mr(dat_clumped, 
          method_list = c("mr_ivw","mr_egger_regression",
                          "mr_weighted_median" ))
  
  ## -- cont mix 
  res_conmix <- mr_conmix(mr_input(bx = dat_clumped$beta.exposure, 
                                   bxse =  dat_clumped$se.exposure, 
                                   by = dat_clumped$beta.outcome, 
                                   byse = dat_clumped$se.outcome))
  
  ## -- MR lasso 
  res_lasso <- mr_lasso(mr_input(bx = dat_clumped$beta.exposure, 
                                 bxse =  dat_clumped$se.exposure, 
                                 by = dat_clumped$beta.outcome, 
                                 byse = dat_clumped$se.outcome))
    
  ## -- MR-PRESSO
  print("MR-PRESSO ...")
  mr_pres <- mr_presso(BetaOutcome = "beta.outcome",
                     BetaExposure = "beta.exposure",
                     SdOutcome = "se.outcome", 
                     SdExposure = "se.exposure",
                     OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                     data = dat_clumped, NbDistribution = 500,
                     SignifThreshold = 0.05, seed = 123)
  

## -- sensitivity analysis ------------------------------------------------- ##
  print("Sensitivity analysis ...")
  het <- mr_heterogeneity(dat_clumped)
  isq <- (het$Q - het$Q_df) *100 / het$Q; names(isq) <- het$method
  pleio <- mr_pleiotropy_test(dat_clumped)

  # -- directionality test
  steiger <- directionality_test(dat_clumped)

  # -- single snp & LOO analysis 
  res_single <- mr_singlesnp(dat_clumped) # single SNP analysis
  # leave-one-out analysis
  res_loo <- mr_leaveoneout(dat_clumped) # leave-one-out analysis


  # -- plots 
  p1 <- mr_scatter_plot(res, dat_clumped);  # scatter plot 
  p2 <- mr_forest_plot(res_single) # forest plot 
  p3 <- mr_leaveoneout_plot(res_loo) # LOO plot 
  p4 <- mr_funnel_plot(res_single)# funnel plot 
  
  ggsave(paste0("plots_",name_exposure[i],"_",name_outcome[j],".pdf"), 
         ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],ncol = 2, nrow = 2), 
         width = 15, height = 10)


## -- combine and save all MR results -------------------------------------- ##
  saveRDS(list(param=param, iv=dat_clumped, mr=res, mr_cm=res_conmix,
               mrlasso=res_lasso, mrp=mr_pres$`Main MR results`,
               het=het, isq=isq, pleiotropy=pleio, steiger=steiger, 
               res_single=res_single, loo=res_loo),
        file = paste0("res_",name_exposure[i],"_",name_outcome[j],".RDS"))
  
}, error = function(e) { 
  skip_to_next <- TRUE
  error_grid$error[index_error] <- skip_to_next})
  
  if(skip_to_next) { 
    index_error <- index_error + 1
    next } else {
      index_error <- index_error + 1
    } 

}
  }, error = function(e) { 
    skip_to_next <- TRUE
    error_grid$error[index_error:(index_error+length(name_exposure))] <- skip_to_next
    })
  
  if(skip_to_next) { 
    index_error <- index_error + length(name_exposure)
    next } 
}

# export error summary 
saveRDS(error_grid, "error_report_batch1.RDS")
