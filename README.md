# mrpdsleep
Mendelian Randomisation - sleep vs PD progression

rscript_MR_pipeline_2_IEU_noIEU.R
- *** main pipeline used for MR study ***
- IV selection
- proxy search using gwasvcf package and European 1000 Genomes Linkage Disequilibrium (LD) reference panel 
- clumping post-harmonisation
- MR analysis incl IVW, Egger, WM, MRPresso, Cont Mixt, MR-Lasso
- Heterogeneity (Q-cochrane score) + Egger y-int
- Scatter plots, LOO analysis

rscript_MR_results_summary.R
- pipeline for forest plot generation
- calculation of F-statistics

rscript_retrieve_snp_rsid.R
- pipeline to change chr:pos formatting to RSID for each variant
- for Iwaki et al PD progression meta-GWAS study

rscript_collider.R
- perform collider bias analysis using SlopeHunter R package by Mahmoud et al

res_XXX.RDS files
- contains collation of MR and sensitivity analysis for each exp vs outcome combination
- RDS files containing collider bias analysis from SlopeHunter are in res_colbias_XXX.RDS files

Further data can be obtained on request.
