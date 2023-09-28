## -- in-house functions for MR analysis 

## -- manual harmonisation - after harmonise_data() function 
manual_harmonise_data <- function(x){
  
  for (i in 1:nrow(x)){
    beta.exp <- x$beta.exposure[i]
    e.allele <- x$effect_allele.exposure[i] 
    o.allele <- x$other_allele.exposure[i]
    if (beta.exp<0){
      x$beta.exposure[i] <- abs(x$beta.exposure[i])
      x$beta.outcome[i] <- -x$beta.outcome[i]
      x$effect_allele.exposure[i] <- o.allele
      x$other_allele.exposure[i] <- e.allele
      x$effect_allele.outcome[i] <- o.allele
      x$other_allele.outcome[i] <- e.allele
    }
  }
  
  return(x)
  
}

## -- manual IVW - MR analysis 
# input x = harmonised and clumped SNP dataset

manual_mr_analysis_ivw <- function(x){
  res <- lm(x$beta.outcome ~ x$beta.exposure -1, 
             weights=x$se.outcome^-2)
  
  return(res)
}

