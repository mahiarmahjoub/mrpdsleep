#!/usr/bin/env Rscript
## ========================================================================= ##
## ======================== Mendelian Randomisation ======================== ##
## =================================TwoSampleMR============================= ##
## ======================= sleep IEU vs PD progression ===================== ##

.libPaths(c(.libPaths(), "~/Documents/myRlib/"))

library(gwasvcf)
library(gwasglue)
#library(biomaRt)
library(TwoSampleMR)
library(MendelianRandomization)
library(magrittr)
#library(phenoscanner)
library(MRPRESSO)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(ggpubr)
library(progress)

gwasvcf::set_bcftools()
gwasvcf::set_plink()

#dir <- c("D://OneDrive - The University of Sydney (Students)/Research/6 PD Sleep/GWAS/")
dir <- c("/Users/macbook/OneDrive - The University of Sydney (Students)/Research/6 PD Sleep/GWAS/")
setwd(dir)
source("rscript_functions_manual_mranalysis.R")

all_files <- list.files(".","_mod.txt")

# filenames
name_exposure <- c("sleepduration","insomnia") #,"chronotype","naps")
name_outcome <- sapply(all_files, 
                       function(x){paste(strsplit(x,"_")[[1]][1:2], 
                                         collapse = "")}, USE.NAMES = F)
name_outcome <- gsub("cont","",name_outcome)

filename_exposure <- c("ukb-b-4424","ukb-b-3957")#,"ukb-b-4956","ukb-b-4616")
filename_outcome <- all_files


## -- setup parameters ----------------------------------------------------- ##
r2_proxy <- 0.8
kb_proxy <- 10000
# clumping - post-harmonisation
clump_kb <- 10000
clump_r2 <- 0.001
# pval cutoff - Phenoscanner SNP search
pheno_pval <- 5e-8
# pval cutoff - IV selection 
iv_pval <- 5e-8
# collate params 
param <- c(r2_proxy=r2_proxy, kb_proxy=kb_proxy, 
           clump_kb=clump_kb, clump_r2=clump_r2, 
           pheno_pval=pheno_pval, iv_pval=iv_pval
)

## == LOOP - MR analysis pipeline ========================================== ##
# initialise progressbar 
pbi <- progress_bar$new(total = length(name_exposure))

for (i in 1:length(name_exposure)){
  pbi$tick()
  print("Importing exposure ...")
  ## -- exposure ----------------------------------------------------------- ##
  # note that clumping and pval threshold is enforced - default atm
  exposure_dat <- extract_instruments(outcomes = filename_exposure[i], 
                                      p1 = iv_pval, clump = F)
  IVs <- exposure_dat$SNP
  
  
  pbj <- progress_bar$new(total = length(name_outcome))
  for (j in 1:length(filename_outcome)){#1:length(name_outcome)){
    pbj$tick()

    ## -- outcome ---------------------------------------------------------- ##
    print(paste0("* ",name_exposure[i]," vs ",name_outcome[j]," *"))
    print("Importing outcome & proxy search ...")
    
    outcome_dat <- read_outcome_data(
      #snps = IVs,
      filename = filename_outcome[j],
      sep = "\t",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "allele1",
      other_allele_col = "allele0",
      eaf_col = "MAF",
      pval_col = "P", 
      chr_col = "chr", 
      pos_col = "pos", 
      samplesize_col = "N"
    )
    
    # remove absent SNPs
    outcome_dat <- subset(outcome_dat, outcome_dat$SNP!="")

    ## -- convert outcome df to vcf 
    out_vcf <- outcome_dat %$%
      create_vcf(chrom=chr.outcome, pos=pos.outcome,
                 nea=other_allele.outcome, ea=effect_allele.outcome,
                 snp=SNP, ea_af=eaf.outcome,
                 effect=beta.outcome, se=se.outcome,
                 #pval=-log10(pval.outcome),
                 pval=pval.outcome,
                 n=samplesize.outcome,
                 name="outcome")
    VariantAnnotation::writeVcf(out_vcf, file=paste0("vcf_",name_exposure[i],"_",
                                  name_outcome[j],".vcf"))
    # note in VCF: LP=-log10(pval)


    ## -- LD proxies - outcome data ---------------------------------------- ##
    out_vcf_proxy <- gwasvcf::query_gwas(paste0("vcf_",name_exposure[i],"_",
                                                     name_outcome[j],".vcf"),
                                         rsid = IVs, 
                                         proxies="yes", 
                                         bfile="EUR",
                                         #bfile="data_maf0.01_rs_ref", 
                                         tag_kb = kb_proxy, 
                                         tag_r2 = r2_proxy)
    outcome_dat_prox <- gwasglue::gwasvcf_to_TwoSampleMR(out_vcf_proxy,
                                                         "outcome")


    ## -- harmonise data - exposure and outcome ---------------------------- ##
    print("Harmonisation ...")
    dat <- harmonise_data(
      exposure_dat = exposure_dat, 
      outcome_dat = outcome_dat_prox,
      action = 3
    )

    # remove palindromic snps prior to clumping
    dat <- subset(dat, mr_keep == TRUE)
    

    ## -- clump final harmonised data -------------------------------------- ##
    print("Clumping ...")
    dat_clumped <- clump_data(dat = dat, clump_kb = clump_kb, 
                              clump_r2 = clump_r2)
    dat_clumped$outcome <- name_outcome[j]
    dat_clumped <- subset(dat_clumped, 
                          !duplicated(dat_clumped)) #remove duplicates

    ## -- manual harmonisation  -------------------------------------------- ##
    dat_clumped <- manual_harmonise_data(dat_clumped)


    ## -- removal of outlier RSIDs - manual -------------------------------- ##
    #dat_clumped <- dat_clumped[-4,]
    
    ## -- MR analysis ------------------------------------------------------ ##
    print("MR & sensitivity analysis ...")
    
    if (nrow(dat_clumped)==1 | nrow(dat_clumped)==0){
      res <- c(); res$b <- rep(NA, 3)
      res$se <- res$pval <- res$b
      
      het <- c(); het$Q <- het$Q_df <- het$Q_pval <- c(NA,NA) 
      isq <- (het$Q - het$Q_df) *100 / het$Q
      
      plt <- c(); plt$egger_intercept <- plt$se <- plt$pval <- NA
      
    } else if (nrow(dat_clumped)==2) { 
      res <-  mr(dat_clumped, 
                 method_list = c("mr_ivw","mr_egger_regression",
                                 "mr_weighted_median"))
      res <- rbind(res, res, res)
      res$b[2:length(methods)] <- NA
      res$se[2:length(methods)] <- NA
      res$pval[2:length(methods)] <- NA
      res$method <- c("WM","Egger","IVW")
      
      plt <- mr_pleiotropy_test(dat_clumped)
      
      het <- mr_heterogeneity(dat_clumped)
      het <- rbind(rep(NA, ncol(het)), het)
      isq <- (het$Q - het$Q_df) *100 / het$Q
      
    } else {
      res <-  mr(dat_clumped, 
                 method_list = c("mr_ivw","mr_egger_regression",
                                 "mr_weighted_median")) 
      
      ## -- heterogeneity & pleiotropy
      plt <- mr_pleiotropy_test(dat_clumped)
      het <- mr_heterogeneity(dat_clumped)
      isq <- (het$Q - het$Q_df) *100 / het$Q; names(isq) <- het$method
    }
    
    res$exposure <- name_exposure[i]
    res$outcome <- name_outcome[j]
    
    ## -- cont mix 
    if(nrow(dat_clumped)>1){
      res_conmix <- mr_conmix(mr_input(bx = dat_clumped$beta.exposure, 
                                       bxse =  dat_clumped$se.exposure, 
                                       by = dat_clumped$beta.outcome, 
                                       byse = dat_clumped$se.outcome))
    } else {
      res_conmix <- NA
    }
    
    ## -- MR lasso 
    if(nrow(dat_clumped)>2){
      res_lasso <- mr_lasso(mr_input(bx = dat_clumped$beta.exposure, 
                                     bxse =  dat_clumped$se.exposure, 
                                     by = dat_clumped$beta.outcome, 
                                     byse = dat_clumped$se.outcome))
    } else {
      res_lasso <- NA 
    }
    
    ## -- MR-PRESSO
    if(nrow(dat_clumped)>3){
      mr_pres <- mr_presso(BetaOutcome = "beta.outcome", 
                           BetaExposure = "beta.exposure",
                           SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                           OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                           data = dat_clumped, NbDistribution = 1000, 
                           SignifThreshold = 0.05, seed = 123)
    } else {
      mr_pres <- NA 
    }
    
    # -- directionality test 
    steiger <- directionality_test(dat_clumped)
    
    # -- single snp & LOO analysis 
    res_single <- mr_singlesnp(dat_clumped) # single SNP analysis
    # leave-one-out analysis
    res_loo <- mr_leaveoneout(dat_clumped) # leave-one-out analysis
    # assign outcome phenotype 
    res_single$outcome <- name_outcome[j]; res_loo$outcome <- name_outcome[j]
    
    # -- plots 
    p1 <- mr_scatter_plot(res, dat_clumped) # scatter plot
    p2 <- mr_forest_plot(res_single) # forest plot
    p3 <- mr_leaveoneout_plot(res_loo) # LOO plot
    p4 <- mr_funnel_plot(res_single) # funnel plot

    ggsave(paste0("plots_",name_exposure[i],"_",name_outcome[j],".pdf"),
           ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],ncol = 2, nrow = 2),
           width = 15, height = 10)
    #with removal rsid
    # ggsave(paste0("plots_",name_exposure[i],"_",name_outcome[j],"_rs16959955rs2236295.pdf"),
    #        ggarrange(p1[[1]],p2[[1]],p3[[1]],p4[[1]],ncol = 2, nrow = 2),
    #        width = 15, height = 10)



## -- combine all MR results ----------------------------------------------- ##
    saveRDS(list(param=param, iv=dat_clumped, mr=res,
                 mrcm=res_conmix,  mrlasso=res_lasso,
                 het=het, isq=isq,
                 pleiotropy=plt, mrp=mr_pres,
                 steiger=steiger,
                 res_single=res_single, loo=res_loo),
            file = paste0("res_",name_exposure[i],"_",name_outcome[j],".RDS"))
    #with removal rsid
    # saveRDS(list(param=param, iv=dat_clumped, mr=res,
    #              mrcm=res_conmix,  mrlasso=res_lasso,
    #              het=het, isq=isq,
    #              pleiotropy=plt, mrp=mr_pres,
    #              steiger=steiger,
    #              res_single=res_single, loo=res_loo),
    #         file = paste0("res_",name_exposure[i],"_",name_outcome[j],"_rs16959955rs2236295.RDS"))


  }
}






