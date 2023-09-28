
library(forester)
library(extrafont)
library(tidyverse)
#https://github.com/rdboyes/forester 


dir <- c("/Users/macbook/OneDrive - The University of Sydney (Students)/Research/6 PD Sleep/GWAS/")
setwd(dir)

table <- read_delim("mr_beta_sleepduration.txt")

# indent the subgroup if there is a number in the placebo column
table$Outcome <- ifelse(is.na(table$`# IVs`), 
                         paste0("   ", table$Outcome),
                         table$Outcome
                         )

# change pval TRUE/FALSE to Yes/No
table$`Stat Sig.` <- ifelse((table$`Stat Sig.`)==TRUE, 
                         "Yes*", "No")

# use forester to create the table with forest plot
forester(left_side_data = table[,c(1,6, 5)], 
         estimate = table$beta,
         ci_low = table$beta_lower,
         ci_high = table$beta_upper,
         display = TRUE, 
         render_as = "pdf",
         xlim = c(-50, 50), 
         font_family = "sans",
         estimate_precision = 1,
         estimate_col_name = "Estimate (95% CI)",
         arrows = TRUE, 
         arrow_labels = c("Protective", "Worsening"), 
         file_path = "plots_forestplot_sleepduration.pdf")
