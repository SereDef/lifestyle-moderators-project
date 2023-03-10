library(mice)
library(miceadds)

# Data location 
if(!exists('datapath')) { datapath <- dirname(file.choose()) } # volumes = X:/Psychology/ResearchProjects/EWalton/ (EarlyCause/WP3/SD/Data)
# Output path 
impqcpath <- file.path(datapath,'QC')

# Date (for marking output file)
date <- format(Sys.Date(), '%d%m%y')

# Moderation imputation file 
MODset <- readRDS(file.path(datapath, 'imp_qp_mod_only_220223.rds'))

# Identify the variables which were actually imputed at this round (the moderators)
imputed <- names(MODset$data[, colSums(is.na(MODset$data)) != 0])
# Select only imputed variables (this uses miceadds to subset and trnasform back and forth)
MODset <- datlist2mids( subset_datlist( mids2datlist(MODset), select=c('IDC',imputed) ))

# Density plot does not work very well with for loops so little trick:
# for (var in imputed) { cat(paste('densityplot(MODset, ~', var, ')', '\n')) }
pdf(file.path(impqcpath, paste0('imp-vs-obs-mods_qp_',date,'.pdf')))
densityplot(MODset, ~ exercise ) 
densityplot(MODset, ~ exercise_bin ) 
densityplot(MODset, ~ exercise_age ) 
densityplot(MODset, ~ exercise_sports_6y ) 
densityplot(MODset, ~ exercise_sports_10y ) 
densityplot(MODset, ~ exercise_walk.cycle_6y ) 
densityplot(MODset, ~ exercise_walk.cycle_10y ) 
densityplot(MODset, ~ exercise_outdoor_4y ) 
densityplot(MODset, ~ exercise_outdoor_10y ) 
densityplot(MODset, ~ sleep_hr ) 
densityplot(MODset, ~ sleep_SR ) 
densityplot(MODset, ~ sleep_hr_age ) 
densityplot(MODset, ~ sleep_hr_bin ) 
densityplot(MODset, ~ sleep_SR_bin ) 
densityplot(MODset, ~ sleep_probs_1.5y ) 
densityplot(MODset, ~ sleep_probs_3y ) 
densityplot(MODset, ~ sleep_probs_6y ) 
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
dev.off()

# ==============================================================================
# ELS project imputed file (with exposure, covariates and outcomes)
ELSset <- readRDS(file.path(datapath, 'imputation_list_sample.rds'))

# Check dimentions
dim(ELSset$data); dim(MODset$data)
# Check overlap in IDCs
identical(ELSset$data$IDC, MODset$data$IDC)
setdiff(ELSset$data$IDC, MODset$data$IDC)
# check <- data.frame(ELSset$data$IDC, MODset$data$IDC)

# Ok these do not correspond so I cannot cbind
ELStmp <- complete(ELSset, 'long', include=T)
# mod <- complete(MODset, 'long', include=T)
order <- MODset$data$IDC

for (m in 0:30) {
  d <- ELStmp[ELStmp$.imp==m, ]
  d <- d %>% dplyr::slice(match(order, d$IDC))
  ELStmp[ELStmp$.imp==m, ] <- d
}
ELSset <- as.mids(ELStmp); rm(ELStmp)

identical(ELSset$data$IDC, MODset$data$IDC)

# Now i can merge in peace
imp = mice:::cbind.mids(ELSset, MODset)

# ==========================================================
# Save 
saveRDS(imp, file.path(datapath, paste0('GenR_imp_sample_merged_',date,'.rds')))
