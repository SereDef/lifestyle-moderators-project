library(dplyr)
library(mice)

# Data location 
# if(!exists('datapath')) { datapath <- dirname(file.choose()) }
datapath <- '~/Desktop/Bath_visit/Data/-GenR/'
# Date (for marking output file)
date <- format(Sys.Date(), '%d%m%y')

# ------------------------------------------------------------------------------
# Take last imputed dataset
impset <- complete(readRDS(file.path(datapath, 'imputation_list_sample.rds')), 30)
# Moderator file
moderat_aux <- readRDS(file.path(datapath, 'GenR_moderator_vars.rds'))

# Merge everything together
data <- merge(impset, moderat_aux, by='IDC', all.x=T)

# Add total ELS variable 
data$ELS <- rowSums(data[,c('prenatal_stress','postnatal_stress')])

# ==============================================================================
# IMPUTATION ===================================================================
# ==============================================================================

imp0 <- mice(data, maxit = 0)
imp0$loggedEvents

# We use default method (pmm for continuous and logisic / polyreg for factors)
meth <- imp0$method # make.method(data, defaultMethod = rep("pmm", 4))

# Setting up passive imputation where appropriate:
# EXERCISE
meth['exercise_bin'] <- '~as.factor(ifelse(exercise %in% c("4 week",">=5 week"), 1, 0))'

# SLEEP
meth['sleep_hr_bin'] <- '~as.factor(ifelse(sleep_hr < 9 | sleep_hr > 12, 0, 1))'
meth['sleep_SR_bin'] <- '~as.factor(ifelse(sleep_SR < 9 | sleep_SR > 12, 0, 1))'

# DIET
diet_formula <- function(var) {
  if (var %in% c('med_diet_mea','med_diet_dai')) { cutoff <- '<' } else { cutoff <- '>=' }
  f = paste("~ifelse(", var, cutoff, "median(",var,"), 1, 0)")
  return(f)
}

meth['med_diet_veg_bin'] <- diet_formula('med_diet_veg')
meth['med_diet_leg_bin'] <- diet_formula('med_diet_leg')
meth['med_diet_fru_bin'] <- diet_formula('med_diet_fru')
meth['med_diet_cer_bin'] <- diet_formula('med_diet_cer')
meth['med_diet_fis_bin'] <- diet_formula('med_diet_fis')
meth['med_diet_mea_bin'] <- diet_formula('med_diet_mea')
meth['med_diet_dai_bin'] <- diet_formula('med_diet_dai')

meth['med_diet'] <- '~I(med_diet_veg_bin + med_diet_leg_bin + med_diet_fru_bin + med_diet_cer_bin + med_diet_fis_bin + med_diet_mea_bin + med_diet_dai_bin)'
meth['med_diet_bin'] <- diet_formula('med_diet')

pm <- mice::quickpred(data, mincor=0.1, exclude =c('IDC', 'twin', 
                      'prenatal_stress_z', 'postnatal_stress_z','intern_score_13_z', 
                      'total_fat_13_z', 'andr_fat_mass_13_z', 'tot_fat_percent_13_z'))
write.csv(pm, file.path(datapath,'QC', paste0('QuickPred_',date,'.csv')))

# Imputation 
imp <- mice(data, method = meth, predictorMatrix = pm, m = 30, maxit = 60) 

# ==========================================================
# Save 
saveRDS(imp, file.path(datapath, paste0('imp_qp_mod_only_',date,'.rds')))
