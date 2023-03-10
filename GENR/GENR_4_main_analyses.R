# Load needed packages
pack <- c('mice','miceadds','nnet', 'openxlsx')
invisible(lapply(pack, require, character.only = T));

#Data location 
if(!exists('datapath')) { datapath <- dirname(file.choose()) } # volumes = X:/Psychology/ResearchProjects/EWalton/ (EarlyCause/WP3/SD/Data)

respath <- file.path(dirname(datapath),'Results')
date <- format(Sys.Date(), "%d%m%y")

sample_imp_set <- readRDS(file.path(datapath,'GenR_imp_sample_merged_220223.rds'))

# ------------------------------------------------------------------------------

# Add missing variables
impdat <- complete(sample_imp_set, action="long", include = T)
# reverse binary moderator variables
impdat$exercise_binR <- ifelse(impdat$exercise_bin == 1, 0, 1) 
impdat$sleep_hr_binR <- ifelse(impdat$sleep_hr_bin == 1, 0, 1)
impdat$sleep_SR_binR <- ifelse(impdat$sleep_SR_bin == 1, 0, 1)
impdat$med_diet_binR <- ifelse(impdat$med_diet_bin == 1, 0, 1)
# Add total ELS score in case it's missing
impdat$ELS <- rowSums(impdat[,c('prenatal_stress','postnatal_stress')])
# Include a numeric version of the exercise factor, to transform into z score 
impdat$exercise_fac <- impdat$exercise
impdat$exercise <- as.numeric(impdat$exercise)

sample_imp_set <- as.mids(impdat); rm(impdat)

# Add standardized scores for main variables
trans <- c('ELS','exercise','sleep_hr','sleep_SR','med_diet')
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
sleep_main <- pool_fit('* sleep_hr_z'); sleepSR_main <- pool_fit('* sleep_SR_z');  
mdiet_main <- pool_fit('* med_diet_z')

# binary moderator analyses with healthy as the reference group
exerc_bin1 <- pool_fit('* exercise_bin')
sleep_bin1 <- pool_fit('* sleep_hr_bin'); sleepSR_bin1 <- pool_fit('* sleep_SR_bin')
mdiet_bin1 <- pool_fit('* med_diet_bin')
# reversed binary moderator analyses with healthy as the reference group
exerc_bin2 <- pool_fit('* exercise_binR')
sleep_bin2 <- pool_fit('* sleep_hr_binR'); sleepSR_bin2 <- pool_fit('* sleep_SR_binR')
mdiet_bin2 <- pool_fit('* med_diet_binR')

# baseline analyses with comorbid as the reference group
exerc_comR <- pool_fit('* exercise_z', outc='risk_groups_perc_REC')
sleep_comR <- pool_fit('* sleep_hr_z', outc='risk_groups_perc_REC'); sleepSR_comR <- pool_fit('* sleep_SR_z', outc='risk_groups_perc_REC')
mdiet_comR <- pool_fit('* med_diet_z', outc='risk_groups_perc_REC')

# prenatal exposure assessment
exerc_pren <- pool_fit('* exercise_z', exp='prenatal_stress_z')
sleep_pren <- pool_fit('* sleep_hr_z', exp='prenatal_stress_z'); sleepSR_pren <- pool_fit('* sleep_SR_z', exp='prenatal_stress_z')
mdiet_pren <- pool_fit('* med_diet_z', exp='prenatal_stress_z')
# postnatal exposure assessment
exerc_post <- pool_fit('* exercise_z', exp='postnatal_stress_z')
sleep_post <- pool_fit('* sleep_hr_z', exp='postnatal_stress_z'); sleepSR_post <- pool_fit('* sleep_SR_z', exp='postnatal_stress_z')
mdiet_post <- pool_fit('* med_diet_z', exp='postnatal_stress_z')

# Internalizing outcome
exerc_intr <- pool_fit('* exercise_z', outc='intern_score_13_z')
sleep_intr <- pool_fit('* sleep_hr_z', outc='intern_score_13_z'); sleepSR_intr <- pool_fit('* sleep_SR_z', outc='intern_score_13_z')
mdiet_intr <- pool_fit('* med_diet_z', outc='intern_score_13_z')
# Fat mass outcome 
exerc_fatm <- pool_fit('* exercise_z', outc='tot_fat_percent_13_z')
sleep_fatm <- pool_fit('* sleep_hr_z', outc='tot_fat_percent_13_z'); sleepSR_fatm <- pool_fit('* sleep_SR_z', outc='tot_fat_percent_13_z')
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

# No interaction models 
exerc_noin <- pool_fit('+ exercise_z')
sleep_noin <- pool_fit('+ sleep_hr_z'); sleepSR_noin <- pool_fit('+ sleep_SR_z')
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
'exerc_fatm_post'= exerc_fatm_post, 'sleep_fatm_post'= sleep_fatm_post, 'mdiet_fatm_post'= mdiet_fatm_post,
'sleepSR_main'= sleepSR_main, 
'sleepSR_bin1'= sleepSR_bin1, 
'sleepSR_bin2'= sleepSR_bin2, 
'sleepSR_intr'= sleepSR_intr, 
'sleepSR_fatm'= sleepSR_fatm, 
'sleepSR_pren'= sleepSR_pren, 
'sleepSR_post'= sleepSR_post, 
'sleepSR_comR'= sleepSR_comR, 
'sleepSR_noin'= sleepSR_noin)

respath <- '~/Desktop/Bath_visit/Results'
openxlsx::write.xlsx(modls, file = file.path(respath, paste0('GENR_Results_',date,'.xlsx')), 
                     overwrite=T)

