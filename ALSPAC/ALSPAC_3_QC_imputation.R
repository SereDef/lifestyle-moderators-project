library(mice)
library(miceadds)

# Data location 
if(!exists('datapath')) { datapath <- dirname(file.choose()) } # volumes = X:/Psychology/ResearchProjects/EWalton/ (EarlyCause/WP3/SD/Data)
# Output path 
impqcpath <- file.path(dirname(datapath),'QC')

# Date (for marking output file)
date <- format(Sys.Date(), '%d%m%y')

# Moderation imputation file 
MODset <- readRDS(file.path(datapath, 'imp_qp_mod_only_210223.rds'))

# Identify the variables which were actually imputed at this round (the moderators)
imputed <- names(MODset$data[, colSums(is.na(MODset$data)) != 0])
# Select only imputed variables (this uses miceadds to subset and trnasform back and forth)
MODset <- datlist2mids( subset_datlist( mids2datlist(MODset), select=c('IDC',imputed) ))

# Density plot does not work very well with for loops so little trick:
# for (var in imputed) { cat(paste('densityplot(MODset, ~', var, ')', '\n')) }
pdf(file.path(impqcpath, paste0('ALSPAC_imp-vs-obs-mods_qp_',date,'.pdf')))
densityplot(MODset, ~ exercise ) 
densityplot(MODset, ~ exercise_bin ) 
densityplot(MODset, ~ exercise_age ) 
densityplot(MODset, ~ exercise_08.2y ) 
densityplot(MODset, ~ exercise_09.7y ) 
densityplot(MODset, ~ exercise_11.7y ) 
densityplot(MODset, ~ exercise_13.2y ) 
densityplot(MODset, ~ exercise_14.7y ) 
densityplot(MODset, ~ exercise_15.4y ) 
densityplot(MODset, ~ exercise_16.1y ) 
densityplot(MODset, ~ exercise_17.0y ) 
densityplot(MODset, ~ sleep_bedtime ) 
densityplot(MODset, ~ sleep_waketime ) 
densityplot(MODset, ~ sleep_hr ) 
densityplot(MODset, ~ sleep_hr_bin ) 
densityplot(MODset, ~ sleep_hr_age ) 
densityplot(MODset, ~ sleep_bedt_04.8y ) 
densityplot(MODset, ~ sleep_waket_04.8y ) 
densityplot(MODset, ~ sleep_04.8y ) 
densityplot(MODset, ~ sleep_bedt_05.8y ) 
densityplot(MODset, ~ sleep_waket_05.8y ) 
densityplot(MODset, ~ sleep_05.8y ) 
densityplot(MODset, ~ sleep_bedt_06.8y ) 
densityplot(MODset, ~ sleep_waket_06.8y ) 
densityplot(MODset, ~ sleep_06.8y ) 
densityplot(MODset, ~ sleep_bedt_09.7y ) 
densityplot(MODset, ~ sleep_waket_09.7y ) 
densityplot(MODset, ~ sleep_09.7y ) 
densityplot(MODset, ~ med_diet ) 
densityplot(MODset, ~ med_diet_bin ) 
densityplot(MODset, ~ med_diet_age ) 
densityplot(MODset, ~ med_diet_veg ) 
densityplot(MODset, ~ med_diet_veg_bin ) 
densityplot(MODset, ~ med_diet_leg ) 
densityplot(MODset, ~ med_diet_leg_bin ) 
densityplot(MODset, ~ med_diet_fru ) 
densityplot(MODset, ~ med_diet_fru_bin ) 
densityplot(MODset, ~ med_diet_cer ) 
densityplot(MODset, ~ med_diet_cer_bin ) 
densityplot(MODset, ~ med_diet_fis ) 
densityplot(MODset, ~ med_diet_fis_bin ) 
densityplot(MODset, ~ med_diet_mea ) 
densityplot(MODset, ~ med_diet_mea_bin ) 
densityplot(MODset, ~ med_diet_dai ) 
densityplot(MODset, ~ med_diet_dai_bin ) 
densityplot(MODset, ~ med_diet_7.5y ) 
densityplot(MODset, ~ med_diet_7.5y_veg ) 
densityplot(MODset, ~ med_diet_7.5y_veg_bin ) 
densityplot(MODset, ~ med_diet_7.5y_leg ) 
densityplot(MODset, ~ med_diet_7.5y_leg_bin ) 
densityplot(MODset, ~ med_diet_7.5y_fru ) 
densityplot(MODset, ~ med_diet_7.5y_fru_bin ) 
densityplot(MODset, ~ med_diet_7.5y_cer ) 
densityplot(MODset, ~ med_diet_7.5y_cer_bin ) 
densityplot(MODset, ~ med_diet_7.5y_fis ) 
densityplot(MODset, ~ med_diet_7.5y_fis_bin ) 
densityplot(MODset, ~ med_diet_7.5y_mea ) 
densityplot(MODset, ~ med_diet_7.5y_mea_bin ) 
densityplot(MODset, ~ med_diet_7.5y_dai ) 
densityplot(MODset, ~ med_diet_7.5y_dai_bin ) 
dev.off()

# ==============================================================================
# ELS project imputed file (with exposure, covariates and outcomes)
ELSset <- readRDS(file.path(datapath, 'imputation_list_sample.rds'))

# Check dimentions
dim(ELSset$data); dim(MODset$data)
# Check overlap in IDCs
identical(ELSset$data$IDC, MODset$data$IDC)
setdiff(ELSset$data$IDC, MODset$data$IDC)

imp = mice:::cbind.mids(ELSset, MODset)

# ==========================================================
# Save 
saveRDS(imp, file.path(datapath, paste0('imp_sample_merged_',date,'.rds')))
