# Load needed packages
pack <- c('mice','miceadds','nnet', 'openxlsx')
invisible(lapply(pack, require, character.only = T));

# Data location 
if(!exists('datapath')) { datapath <- dirname(file.choose()) } # volumes = X:/Psychology/ResearchProjects/EWalton/ (EarlyCause/WP3/SD/Data)
# Output path 
respath <- file.path(dirname(datapath),'Results')
# Date (for marking output file)
date <- format(Sys.Date(), "%d%m%y")

sample_imp_set <- readRDS(file.path(datapath,'imp_sample_merged_220223.rds'))

# ------------------------------------------------------------------------------

# Add missing variables (reverse coded binary moderators)
impdat <- complete(sample_imp_set, action="long", include = T)
# reverse binary moderator variables
impdat$exercise_binR <- ifelse(impdat$exercise_bin == 1, 0, 1) 
impdat$sleep_hr_binR <- ifelse(impdat$sleep_hr_bin == 1, 0, 1)
impdat$med_diet_binR <- ifelse(impdat$med_diet_bin == 1, 0, 1)
# Add total ELS score in case it's missing
impdat$ELS <- rowSums(impdat[,c('prenatal_stress','postnatal_stress')])
# Include a numeric version of the exercise factor, to transform into z score 
impdat$exercise_fac <- impdat$exercise
impdat$exercise <- as.numeric(impdat$exercise)

sample_imp_set <- as.mids(impdat); rm(impdat)

# Add standardized scores for main variables
domains <- c('pre_life_events', 'pre_contextual_risk', 'pre_parental_risk', 'pre_interpersonal_risk', 
             'post_life_events', 'post_contextual_risk', 'post_parental_risk', 'post_interpersonal_risk', 'post_direct_victimization')
foodgroups <- c('med_diet_veg','med_diet_leg','med_diet_fru','med_diet_cer','med_diet_fis','med_diet_mea','med_diet_dai')

trans <- c('ELS','exercise','sleep_hr','med_diet', domains, foodgroups)
sample_imp_set <- datlist2mids( scale_datlist( mids2datlist(sample_imp_set), 
                                 orig_var = trans, trafo_var = paste0(trans, '_z')))

# Save each imputed setfor portenial plottining
for (m in 0:30) { write.csv(complete(sample_imp_set,m), file.path(datapath,'byimp',paste0('imp',m,'.csv'))) }

# ------------------------------------------------------------------------------

pool_fit <- function(modr, outc='risk_groups_perc', exp='ELS_z') {
  # Define covariates
  covs <- '+ age_child + sex + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking'
  # Define outcome recoding status for model name
  # rec <- ifelse(outc=='risk_groups_perc', 'H', 'M') # assign(paste('base',name,rec,'ref',sep='_'), 
  # IF risk outcome fit a multimod model else lm
  if (startsWith(outc, 'risk')) {
    fit <- with(sample_imp_set, nnet::multinom(as.formula(paste(outc,'~',exp,modr,covs)), model = T));
  } else { 
    fit <- with(sample_imp_set, lm(as.formula(paste(outc,'~',exp,modr,covs)))); 
  }
  p_fit <- mice::pool(fit) # pool results 
  mod <- summary(p_fit) # extract relevant information
  mod[,-c(1,2)] <- round(mod[,-c(1,2)],4)
  mod$sign <- ifelse(mod$p.value < 0.05, '*', '') # add a column to highlight significant terms
  
  if (startsWith(outc, 'risk')) {
    # make group comparisons easier to read
    if (endsWith(outc, 'REC')) { 
      levels(mod$y.level) <- c("M:healthy","M:intern", "M:fatmas")
    } else { levels(mod$y.level) <- c("H:intern", "H:fatmas", "H:multim") }
    mod$OR  <- round(exp(mod$estimate), 4)
    mod$lci <- round(exp((mod$estimate) - 1.96*mod$std.error), 4)
    mod$uci <- round(exp((mod$estimate) + 1.96*mod$std.error), 4)
    mod$AIC <- c(mean(p_fit$glanced$AIC), rep(NA, nrow(mod)-1)) # add a column for AIC
  } else {
    mod$lci <- round((mod$estimate - 1.96*mod$std.error), 4)
    mod$uci <- round((mod$estimate + 1.96*mod$std.error), 4)
    mod$rsq <- c(pool.r.squared(fit)[1], rep(NA, nrow(mod)-1)) # add a column for R2
    mod$rsq_adj <- c(pool.r.squared(fit, adjusted = TRUE)[1], rep(NA, nrow(mod)-1)) # adjusted R2
  }
  # print(mod)
  return(mod)
}
# ------------------------------------------------------------------------------

# baseline analyses with healthy as the reference group
exerc_main <- pool_fit('* exercise_z')
sleep_main <- pool_fit('* sleep_hr_z')
mdiet_main <- pool_fit('* med_diet_z')

# binary moderator analyses with healthy as the reference group
exerc_bin1 <- pool_fit('* exercise_bin')
sleep_bin1 <- pool_fit('* sleep_hr_bin')
mdiet_bin1 <- pool_fit('* med_diet_bin')
# reversed binary moderator analyses with healthy as the reference group
exerc_bin2 <- pool_fit('* exercise_binR')
sleep_bin2 <- pool_fit('* sleep_hr_binR')
mdiet_bin2 <- pool_fit('* med_diet_binR')

# baseline analyses with comorbid as the reference group
exerc_comR <- pool_fit('* exercise_z', outc='risk_groups_perc_REC')
sleep_comR <- pool_fit('* sleep_hr_z', outc='risk_groups_perc_REC')
mdiet_comR <- pool_fit('* med_diet_z', outc='risk_groups_perc_REC')

# prenatal exposure assessment
exerc_pren <- pool_fit('* exercise_z', exp='prenatal_stress_z')
sleep_pren <- pool_fit('* sleep_hr_z', exp='prenatal_stress_z')
mdiet_pren <- pool_fit('* med_diet_z', exp='prenatal_stress_z')
# postnatal exposure assessment
exerc_post <- pool_fit('* exercise_z', exp='postnatal_stress_z')
sleep_post <- pool_fit('* sleep_hr_z', exp='postnatal_stress_z')
mdiet_post <- pool_fit('* med_diet_z', exp='postnatal_stress_z')

# Internalizing outcome
exerc_intr <- pool_fit('* exercise_z', outc='intern_score_13_z')
sleep_intr <- pool_fit('* sleep_hr_z', outc='intern_score_13_z')
mdiet_intr <- pool_fit('* med_diet_z', outc='intern_score_13_z')
# Fat mass outcome 
exerc_fatm <- pool_fit('* exercise_z', outc='tot_fat_percent_13_z')
sleep_fatm <- pool_fit('* sleep_hr_z', outc='tot_fat_percent_13_z')
mdiet_fatm <- pool_fit('* med_diet_z', outc='tot_fat_percent_13_z')

# Internalizing outcome binary moderator 
exerc_intr_bin1 <- pool_fit('* exercise_bin', outc='intern_score_13_z')
sleep_intr_bin1 <- pool_fit('* sleep_hr_bin', outc='intern_score_13_z')
mdiet_intr_bin1 <- pool_fit('* med_diet_bin', outc='intern_score_13_z')
# Internalizing outcome reversed binary moderator
exerc_intr_bin2 <- pool_fit('* exercise_binR', outc='intern_score_13_z')
sleep_intr_bin2 <- pool_fit('* sleep_hr_binR', outc='intern_score_13_z')
mdiet_intr_bin2 <- pool_fit('* med_diet_binR', outc='intern_score_13_z')

# Fat mass outcome binary moderator 
exerc_fatm_bin1 <- pool_fit('* exercise_bin', outc='tot_fat_percent_13_z')
sleep_fatm_bin1 <- pool_fit('* sleep_hr_bin', outc='tot_fat_percent_13_z')
mdiet_fatm_bin1 <- pool_fit('* med_diet_bin', outc='tot_fat_percent_13_z')
# Fat mass outcome reversed binary moderator
exerc_fatm_bin2 <- pool_fit('* exercise_binR', outc='tot_fat_percent_13_z')
sleep_fatm_bin2 <- pool_fit('* sleep_hr_binR', outc='tot_fat_percent_13_z')
mdiet_fatm_bin2 <- pool_fit('* med_diet_binR', outc='tot_fat_percent_13_z')

# Internalizing outcome - prenatal stress
exerc_intr_pren <- pool_fit('* exercise_z', outc='intern_score_13_z',exp='prenatal_stress_z')
sleep_intr_pren <- pool_fit('* sleep_hr_z', outc='intern_score_13_z',exp='prenatal_stress_z')
mdiet_intr_pren <- pool_fit('* med_diet_z', outc='intern_score_13_z',exp='prenatal_stress_z')
# Fat mass outcome - prenatal stress
exerc_fatm_pren <- pool_fit('* exercise_z', outc='tot_fat_percent_13_z',exp='prenatal_stress_z')
sleep_fatm_pren <- pool_fit('* sleep_hr_z', outc='tot_fat_percent_13_z',exp='prenatal_stress_z')
mdiet_fatm_pren <- pool_fit('* med_diet_z', outc='tot_fat_percent_13_z',exp='prenatal_stress_z')

# Internalizing outcome - postatal stress
exerc_intr_post <- pool_fit('* exercise_z', outc='intern_score_13_z',exp='postnatal_stress_z')
sleep_intr_post <- pool_fit('* sleep_hr_z', outc='intern_score_13_z',exp='postnatal_stress_z')
mdiet_intr_post <- pool_fit('* med_diet_z', outc='intern_score_13_z',exp='postnatal_stress_z')
# Fat mass outcome - postatal stress
exerc_fatm_post <- pool_fit('* exercise_z', outc='tot_fat_percent_13_z',exp='postnatal_stress_z')
sleep_fatm_post <- pool_fit('* sleep_hr_z', outc='tot_fat_percent_13_z',exp='postnatal_stress_z')
mdiet_fatm_post <- pool_fit('* med_diet_z', outc='tot_fat_percent_13_z',exp='postnatal_stress_z')

# ----- exploratory followup ---------------------------------------------------
for (d in domains) {
  assign(paste0('exerc_',d), pool_fit('* exercise_z',exp=paste0(d,'_z')))
  assign(paste0('exerc_intr',d), pool_fit('* exercise_z',outc='intern_score_13_z',exp=paste0(d,'_z')))
  assign(paste0('exerc_fatm',d), pool_fit('* exercise_z', outc='tot_fat_percent_13_z',exp=paste0(d,'_z')))
  
  assign(paste0('sleep_',d), pool_fit('* sleep_hr_z',exp=paste0(d,'_z')))
  assign(paste0('sleep_intr',d), pool_fit('* sleep_hr_z',outc='intern_score_13_z',exp=paste0(d,'_z')))
  assign(paste0('sleep_fatm',d), pool_fit('* sleep_hr_z', outc='tot_fat_percent_13_z',exp=paste0(d,'_z')))
  
  assign(paste0('mdiet_',d), pool_fit('* med_diet_z',exp=paste0(d,'_z')))
  assign(paste0('mdiet_intr',d), pool_fit('* med_diet_z',outc='intern_score_13_z',exp=paste0(d,'_z')))
  assign(paste0('mdiet_fatm',d), pool_fit('* med_diet_z', outc='tot_fat_percent_13_z',exp=paste0(d,'_z')))
}

for (f in foodgroups) {
  assign(paste0(f), pool_fit(paste0('* ',f,'_z')))
  assign(paste0(f,'_intr'), pool_fit(paste0('* ',f,'_z'),outc='intern_score_13_z'))
  assign(paste0(f,'_fatm'), pool_fit(paste0('* ',f,'_z'),outc='tot_fat_percent_13_z'))
}

# No interaction models --------------------------------------------------------
exerc_noin <- pool_fit('+ exercise_z')
sleep_noin <- pool_fit('+ sleep_hr_z')
mdiet_noin <- pool_fit('+ med_diet_z')

# ------------------------------------------------------------------------------
# ls()[grepl('exerc|sleep|mdiet', ls())]
# for (m in c('main','bin1','bin2','intr','fatm','pren','post','comR','noin')) {
#   names = paste(c('exerc','sleep','mdiet'), m, sep="_")
#   for (n in names) { cat("'",n,"'= ",n,", ", sep='')}
# }

modls <- list(
  'exerc_main'= exerc_main, 'sleep_main'= sleep_main, 'mdiet_main'= mdiet_main, 
  'exerc_bin1'= exerc_bin1, 'sleep_bin1'= sleep_bin1, 'mdiet_bin1'= mdiet_bin1, 
  'exerc_bin2'= exerc_bin2, 'sleep_bin2'= sleep_bin2, 'mdiet_bin2'= mdiet_bin2, 
  'exerc_intr'= exerc_intr, 'sleep_intr'= sleep_intr, 'mdiet_intr'= mdiet_intr, 
  'exerc_fatm'= exerc_fatm, 'sleep_fatm'= sleep_fatm, 'mdiet_fatm'= mdiet_fatm, 
  'exerc_pren'= exerc_pren, 'sleep_pren'= sleep_pren, 'mdiet_pren'= mdiet_pren, 
  'exerc_post'= exerc_post, 'sleep_post'= sleep_post, 'mdiet_post'= mdiet_post, 
  'exerc_comR'= exerc_comR, 'sleep_comR'= sleep_comR, 'mdiet_comR'= mdiet_comR, 
  'exerc_noin'= exerc_noin, 'sleep_noin'= sleep_noin, 'mdiet_noin'= mdiet_noin,
  'exerc_intr_bin1'= exerc_intr_bin1, 'sleep_intr_bin1'= sleep_intr_bin1, 'mdiet_intr_bin1'= mdiet_intr_bin1, 
  'exerc_intr_bin2'= exerc_intr_bin2, 'sleep_intr_bin2'= sleep_intr_bin2, 'mdiet_intr_bin2'= mdiet_intr_bin2, 
  'exerc_fatm_bin1'= exerc_fatm_bin1, 'sleep_fatm_bin1'= sleep_fatm_bin1, 'mdiet_fatm_bin1'= mdiet_fatm_bin1, 
  'exerc_fatm_bin2'= exerc_fatm_bin2, 'sleep_fatm_bin2'= sleep_fatm_bin2, 'mdiet_fatm_bin2'= mdiet_fatm_bin2, 
  'exerc_intr_pren'= exerc_intr_pren, 'sleep_intr_pren'= sleep_intr_pren, 'mdiet_intr_pren'= mdiet_intr_pren, 
  'exerc_intr_post'= exerc_intr_post, 'sleep_intr_post'= sleep_intr_post, 'mdiet_intr_post'= mdiet_intr_post,
  'exerc_fatm_pren'= exerc_fatm_pren, 'sleep_fatm_pren'= sleep_fatm_pren, 'mdiet_fatm_pren'= mdiet_fatm_pren, 
  'exerc_fatm_post'= exerc_fatm_post, 'sleep_fatm_post'= sleep_fatm_post, 'mdiet_fatm_post'= mdiet_fatm_post)

openxlsx::write.xlsx(modls, file = file.path(respath, paste0('ALSPAC_Results_',date,'.xlsx')), 
                     overwrite=T)
# ------------------------------------------------------------------------------
# names = ls()[grepl('exerc|sleep|diet', ls())]
# for (n in names) { cat("'",n,"'= ",n,",\n", sep='')}
suppl <- list(
  'med_diet_cer'= med_diet_cer,'med_diet_cer_fatm'= med_diet_cer_fatm,'med_diet_cer_intr'= med_diet_cer_intr,
  'med_diet_dai'= med_diet_dai,'med_diet_dai_fatm'= med_diet_dai_fatm,'med_diet_dai_intr'= med_diet_dai_intr,
  'med_diet_fis'= med_diet_fis,'med_diet_fis_fatm'= med_diet_fis_fatm,'med_diet_fis_intr'= med_diet_fis_intr,
  'med_diet_fru'= med_diet_fru,'med_diet_fru_fatm'= med_diet_fru_fatm,'med_diet_fru_intr'= med_diet_fru_intr,
  'med_diet_leg'= med_diet_leg,'med_diet_leg_fatm'= med_diet_leg_fatm,'med_diet_leg_intr'= med_diet_leg_intr,
  'med_diet_mea'= med_diet_mea,'med_diet_mea_fatm'= med_diet_mea_fatm,'med_diet_mea_intr'= med_diet_mea_intr,
  'med_diet_veg'= med_diet_veg,'med_diet_veg_fatm'= med_diet_veg_fatm,'med_diet_veg_intr'= med_diet_veg_intr,
  'exerc_fatmpost_contextual_risk'= exerc_fatmpost_contextual_risk,
  'exerc_fatmpost_direct_v'= exerc_fatmpost_direct_victimization,
  'exerc_fatmpost_interp_risk'= exerc_fatmpost_interpersonal_risk,
  'exerc_fatmpost_life_events'= exerc_fatmpost_life_events,
  'exerc_fatmpost_parental_risk'= exerc_fatmpost_parental_risk,
  'exerc_fatmpre_contextual_risk'= exerc_fatmpre_contextual_risk,
  'exerc_fatmpre_interp_risk'= exerc_fatmpre_interpersonal_risk,
  'exerc_fatmpre_life_events'= exerc_fatmpre_life_events,
  'exerc_fatmpre_parental_risk'= exerc_fatmpre_parental_risk,
  'exerc_intrpost_contextual_risk'= exerc_intrpost_contextual_risk,
  'exerc_intrpost_direct_v'= exerc_intrpost_direct_victimization,
  'exerc_intrpost_interp_risk'= exerc_intrpost_interpersonal_risk,
  'exerc_intrpost_life_events'= exerc_intrpost_life_events,
  'exerc_intrpost_parental_risk'= exerc_intrpost_parental_risk,
  'exerc_intrpre_contextual_risk'= exerc_intrpre_contextual_risk,
  'exerc_intrpre_interp_risk'= exerc_intrpre_interpersonal_risk,
  'exerc_intrpre_life_events'= exerc_intrpre_life_events,
  'exerc_intrpre_parental_risk'= exerc_intrpre_parental_risk,
  'exerc_post_contextual_risk'= exerc_post_contextual_risk,
  'exerc_post_direct_v'= exerc_post_direct_victimization,
  'exerc_post_interp_risk'= exerc_post_interpersonal_risk,
  'exerc_post_life_events'= exerc_post_life_events,
  'exerc_post_parental_risk'= exerc_post_parental_risk,
  'exerc_pre_contextual_risk'= exerc_pre_contextual_risk,
  'exerc_pre_interp_risk'= exerc_pre_interpersonal_risk,
  'exerc_pre_life_events'= exerc_pre_life_events,
  'exerc_pre_parental_risk'= exerc_pre_parental_risk,
  'mdiet_fatmpost_contextual_risk'= mdiet_fatmpost_contextual_risk,
  'mdiet_fatmpost_direct_v'= mdiet_fatmpost_direct_victimization,
  'mdiet_fatmpost_interp_risk'= mdiet_fatmpost_interpersonal_risk,
  'mdiet_fatmpost_life_events'= mdiet_fatmpost_life_events,
  'mdiet_fatmpost_parental_risk'= mdiet_fatmpost_parental_risk,
  'mdiet_fatmpre_contextual_risk'= mdiet_fatmpre_contextual_risk,
  'mdiet_fatmpre_interp_risk'= mdiet_fatmpre_interpersonal_risk,
  'mdiet_fatmpre_life_events'= mdiet_fatmpre_life_events,
  'mdiet_fatmpre_parental_risk'= mdiet_fatmpre_parental_risk,
  'mdiet_intrpost_contextual_risk'= mdiet_intrpost_contextual_risk,
  'mdiet_intrpost_direct_v'= mdiet_intrpost_direct_victimization,
  'mdiet_intrpost_interp_risk'= mdiet_intrpost_interpersonal_risk,
  'mdiet_intrpost_life_events'= mdiet_intrpost_life_events,
  'mdiet_intrpost_parental_risk'= mdiet_intrpost_parental_risk,
  'mdiet_intrpre_contextual_risk'= mdiet_intrpre_contextual_risk,
  'mdiet_intrpre_interp_risk'= mdiet_intrpre_interpersonal_risk,
  'mdiet_intrpre_life_events'= mdiet_intrpre_life_events,
  'mdiet_intrpre_parental_risk'= mdiet_intrpre_parental_risk,
  'mdiet_post_contextual_risk'= mdiet_post_contextual_risk,
  'mdiet_post_direct_v'= mdiet_post_direct_victimization,
  'mdiet_post_interp_risk'= mdiet_post_interpersonal_risk,
  'mdiet_post_life_events'= mdiet_post_life_events,
  'mdiet_post_parental_risk'= mdiet_post_parental_risk,
  'mdiet_pre_contextual_risk'= mdiet_pre_contextual_risk,
  'mdiet_pre_interp_risk'= mdiet_pre_interpersonal_risk,
  'mdiet_pre_life_events'= mdiet_pre_life_events,
  'mdiet_pre_parental_risk'= mdiet_pre_parental_risk,
  'sleep_fatmpost_contextual_risk'= sleep_fatmpost_contextual_risk,
  'sleep_fatmpost_direct_v'= sleep_fatmpost_direct_victimization,
  'sleep_fatmpost_interp_risk'= sleep_fatmpost_interpersonal_risk,
  'sleep_fatmpost_life_events'= sleep_fatmpost_life_events,
  'sleep_fatmpost_parental_risk'= sleep_fatmpost_parental_risk,
  'sleep_fatmpre_contextual_risk'= sleep_fatmpre_contextual_risk,
  'sleep_fatmpre_interp_risk'= sleep_fatmpre_interpersonal_risk,
  'sleep_fatmpre_life_events'= sleep_fatmpre_life_events,
  'sleep_fatmpre_parental_risk'= sleep_fatmpre_parental_risk,
  'sleep_intrpost_contextual_risk'= sleep_intrpost_contextual_risk,
  'sleep_intrpost_direct_v'= sleep_intrpost_direct_victimization,
  'sleep_intrpost_interp_risk'= sleep_intrpost_interpersonal_risk,
  'sleep_intrpost_life_events'= sleep_intrpost_life_events,
  'sleep_intrpost_parental_risk'= sleep_intrpost_parental_risk,
  'sleep_intrpre_contextual_risk'= sleep_intrpre_contextual_risk,
  'sleep_intrpre_interp_risk'= sleep_intrpre_interpersonal_risk,
  'sleep_intrpre_life_events'= sleep_intrpre_life_events,
  'sleep_intrpre_parental_risk'= sleep_intrpre_parental_risk,
  'sleep_post_contextual_risk'= sleep_post_contextual_risk,
  'sleep_post_direct_v'= sleep_post_direct_victimization,
  'sleep_post_interp_risk'= sleep_post_interpersonal_risk,
  'sleep_post_life_events'= sleep_post_life_events,
  'sleep_post_parental_risk'= sleep_post_parental_risk,
  'sleep_pre_contextual_risk'= sleep_pre_contextual_risk,
  'sleep_pre_interp_risk'= sleep_pre_interpersonal_risk,
  'sleep_pre_life_events'= sleep_pre_life_events,
  'sleep_pre_parental_risk'= sleep_pre_parental_risk)

openxlsx::write.xlsx(suppl, file = file.path(respath, paste0('ALSPAC_SupplementaryResults_',date,'.xlsx')), 
                     overwrite=T)
