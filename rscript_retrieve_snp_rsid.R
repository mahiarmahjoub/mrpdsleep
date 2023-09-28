#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## ================ relabelleing SNPs to rsID - GWAS ss ==================== ##

# need to have two files in the same directory: ENSEMBL + ref.txt

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least two arguments needed: sumstats + ref + outcome/exposure", 
       call.=FALSE)}
# } else if (length(args)==1) {
#   # default output file
#   args[2] = "out.txt"
# }

#.libPaths(c(.libPaths(), "~/Documents/myRlib/"))
#library(biomaRt)
#library(progress)
#n <- 100
dir <- c("/Users/macbook/OneDrive - The University of Sydney (Students)/Research/6 PD Sleep/GWAS/")
#dir <- c("D://OneDrive - The University of Sydney (Students)/Research/6 PD Sleep/GWAS/")


## -- setup biomaRt -------------------------------------------------------- ##
#ensembl <- readRDS("ensemble.RDS")
#ensembl <- useEnsembl(biomart = "snps", dataset = "hsapiens_snp")

## -- import summary stats dataset ----------------------------------------- ##
print("Importing ...")
sumstat <- read.table(args[1], header = 1, stringsAsFactors = F)
#refdat <- read.table(args[2], header = 1, stringsAsFactors = F, sep = ",")
refdat <- readRDS(args[2])

# only include refdat present in sumstats 
refdat_mod <- subset(refdat, SNP %in% sumstat$SNP)


# separate ref to those with and without RSID 
#refdat_mod_rsid <- subset(refdat_mod, RSID!="")
#refdat_mod_norsid <- subset(refdat_mod, RSID=="")


## -- obtain RSID ---------------------------------------------------------- ##
# print("Obtaining RSID ...")
# coords <- apply(refdat_mod_norsid[,c("CHR","START")], 
#                 1, function(x){paste(x,collapse = ":")})
# coords <- paste0(coords, ":", refdat_mod_norsid$START)
# 
# # obtain rsid from database
# index1 <- 1; index2 <- n
# new_rsid <- rep(NA, nrow(refdat_mod_norsid))
# 
# pb <- progress_bar$new(total = ((nrow(refdat_mod_norsid) %/% n)+1))
# for (i in 1:((nrow(refdat_mod_norsid) %/% n)+1)){
#   
#   if(i==((nrow(refdat_mod_norsid) %/% n)+1)){index2<-nrow(refdat_mod_norsid)}
#   
#   temp_rsid <- sapply(coords[index1:index2], function(x){
#     getBM(attributes = c('refsnp_id'), 
#           filters = 'chromosomal_region',
#           values = x, mart = ensembl)})
#   
#   # paste rsid if multiple present per snp_location
#   new_rsid[index1:index2] <- unlist(lapply(temp_rsid, 
#                                            function(x){paste(x,
#                                                              collapse = ":")}),
#                                     use.names = F)
# 
#   index1 <- index1 + n; index2 <- index2 + n; pb$tick()
# 
# }
# 
# 
# # new_rsid <- lapply(new_rsid, function(x){
# #   if(length(x)!=1){NA}else{x}}
# # )
# 
# refdat_mod_norsid$RSID<- new_rsid

## -- putting REF dataset back together ------------------------------------ ##
# refdat_mod <- rbind(refdat_mod_rsid, refdat_mod_norsid)


## -- reformat sumstats ---------------------------------------------------- ##
print("Reformating ...")
sumstat <- cbind(sumstat, 
                 chrpos = sumstat$SNP, 
                 allele1 = rep(NA,nrow(sumstat)), 
                 allele0 = rep(NA,nrow(sumstat)),
                 MAF = rep(NA,nrow(sumstat)), 
                 chr = rep(NA,nrow(sumstat)), 
                 pos = rep(NA,nrow(sumstat))
                 )
                 #type = rep(args[3],nrow(sumstat)))

# change rowname to SNP_position 
rownames(sumstat) <- sumstat$SNP
# order sumstat similar to refdata (asc SNP_position)
sumstat <- sumstat[refdat_mod$SNP,]

## -- obtain rsid, allele0/1, MAF ------------------------------------------ ##
sumstat$SNP <- refdat_mod$RSID
sumstat$allele1 <- refdat_mod$ALT
sumstat$allele0 <- refdat_mod$REF
sumstat$MAF <- refdat_mod$MAF
sumstat$chr <- refdat_mod$CHR
sumstat$pos <- refdat_mod$START

## -- extract chr and pos into individual columns -------------------------- ##
# sumstat$chr <- sapply(sumstat$chrpos, function(x){
#   as.numeric(strsplit(x,":")[[1]][1])
# })
# sumstat$pos <- sapply(sumstat$chrpos, function(x){
#   as.numeric(strsplit(x,":")[[1]][2])
# })

## -- find missing rsid ---------------------------------------------------- ##
# row index of snps with missing rsid 
#missing_rsid_index <- which(sumstat$SNP=="")
# place SNP_pos instead of missing rsid 
#sumstat$SNP[missing_rsid_index] <- sumstat$pos[missing_rsid_index]

## -- export sumstat file -------------------------------------------------- ##
print("Writing ...")
#"Iwaki_PDprogression_MOCA_sumstat.txt"
new_name <- gsub("\\.","_mod.", args[1])
write.table(sumstat, file=new_name, quote = F, row.names = F, 
            col.names = T, sep = "\t")

print("C'est fini.")




