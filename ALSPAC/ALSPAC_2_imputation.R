
library(dplyr)
library(mice)

# Data location 
if(!exists('datapath')) { datapath <- dirname(file.choose()) } # volumes = X:/Psychology/ResearchProjects/EWalton/ (EarlyCause/WP3/SD/Data)

# Date (for marking output file)
date <- format(Sys.Date(), '%d%m%y')

# ------------------------------------------------------------------------------
# Take last (30th) imputed dataset
impset <- complete(readRDS(file.path(datapath, 'imputation_list_sample.rds')), 30)
# Moderator file
moderat_aux <- readRDS(file.path(datapath, 'moderator_vars.rds'))

# Merge everything together
data <- merge(impset, moderat_aux, by='IDC', all.x=T)

# Reorder data so that merging after imputation is easy 
# data <- data %>% dplyr::slice(match(impset$IDC, data$IDC))

# Add total ELS variable 
data$ELS <- rowSums(data[,c('prenatal_stress','postnatal_stress')])

# ==============================================================================
# IMPUTATION ===================================================================
# ==============================================================================
# Initiate
imp0 <- mice(data, maxit = 0)
imp0$loggedEvents

# We use default method (pmm for continuous and logisic / polyreg for factors)
meth <- imp0$method

# Setting up passive imputation where appropriate:
# EXERCISE
meth['exercise_bin'] <- '~as.factor(ifelse(exercise %in% c("4-6 week","daily"), 1, 0))'

# SLEEP
meth['sleep_hr'] <- '~I(24 - sleep_bedtime + sleep_waketime)'
meth['sleep_hr_bin'] <- '~as.factor(ifelse(sleep_hr < 9 | sleep_hr > 12, 0, 1))'

# DIET
diet_formula <- function(var) {
  if (var %in% c('med_diet_mea','med_diet_dai')) { cutoff <- '<' } else if (var=='med_diet_fis') { cutoff <- '>' 
  } else { cutoff <- '>=' }
  f = paste("~ifelse(", var, cutoff, "median(",var,"), 1, 0)")
  return(f)
}

meth['med_diet_7.5y_veg_bin'] <- diet_formula('med_diet_7.5y_veg')
meth['med_diet_7.5y_leg_bin'] <- diet_formula('med_diet_7.5y_leg')
meth['med_diet_7.5y_fru_bin'] <- diet_formula('med_diet_7.5y_fru')
meth['med_diet_7.5y_cer_bin'] <- diet_formula('med_diet_7.5y_cer')
meth['med_diet_7.5y_fis_bin'] <- diet_formula('med_diet_7.5y_fis')
meth['med_diet_7.5y_mea_bin'] <- diet_formula('med_diet_7.5y_mea')
meth['med_diet_7.5y_dai_bin'] <- diet_formula('med_diet_7.5y_dai')

meth['med_diet_7.5y'] <- '~I(med_diet_7.5y_veg_bin + med_diet_7.5y_leg_bin + med_diet_7.5y_fru_bin + med_diet_7.5y_cer_bin + med_diet_7.5y_fis_bin + med_diet_7.5y_mea_bin + med_diet_7.5y_dai_bin)'

meth['med_diet_veg_bin'] <- diet_formula('med_diet_veg')
meth['med_diet_leg_bin'] <- diet_formula('med_diet_leg')
meth['med_diet_fru_bin'] <- diet_formula('med_diet_fru')
meth['med_diet_cer_bin'] <- diet_formula('med_diet_cer')
meth['med_diet_fis_bin'] <- diet_formula('med_diet_fis')
meth['med_diet_mea_bin'] <- diet_formula('med_diet_mea')
meth['med_diet_dai_bin'] <- diet_formula('med_diet_dai')

meth['med_diet'] <- '~I(med_diet_veg_bin + med_diet_leg_bin + med_diet_fru_bin + med_diet_cer_bin + med_diet_fis_bin + med_diet_mea_bin + med_diet_dai_bin)'
meth['med_diet_bin'] <- diet_formula('med_diet')

# Rely on QuickPred to determine the regression models
pm <- mice::quickpred(data, mincor=0.1, exclude =c('IDC', 'twin', 'sibling', 
                      'prenatal_stress_z', 'postnatal_stress_z', 'intern_score_13_z', 
                      'total_fat_13_z', 'andr_fat_mass_13_z', 'tot_fat_percent_13_z'))
write.csv(pm, file.path(dirname(datapath),'QC',paste0('QuickPred_',date,'.csv')))

# Imputation 
imp <- mice(data, method = meth, qppredictorMatrix = pm, m = 30, maxit = 60) 

# ==========================================================
# Save 
saveRDS(imp, file.path(datapath, paste0('imp_qp_mod_only_',date,'.rds')))
