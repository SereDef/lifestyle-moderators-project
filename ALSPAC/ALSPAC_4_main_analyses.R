# Load needed packages
pack <- c('mice','stringr','miceadds','nnet','MNLpred','openxlsx')
invisible(lapply(pack, require, character.only = T));

# Data location 
if(!exists('datapath')) { datapath <- dirname(file.choose()) } # volumes = X:/Psychology/ResearchProjects/EWalton/ (EarlyCause/WP3/SD/Data)
# Output path 
respath <- file.path(dirname(datapath),'Results')
# Date (for marking output file)
date <- format(Sys.Date(), "%d%m%y")

sample_mids <- readRDS(file.path(datapath,'imp_sample_merged_220223.rds'))

# ------------------------------------------------------------------------------

# Add missing variables (reverse coded binary moderators)
longdat <- complete(sample_mids, action="long", include = T)

# Simplify names 
rep_str = c('exercise'='exerc','sleep_hr'='sleep','med_diet'='mdiet', 
            'intern_score_13'='intern', 'tot_fat_percent_13'='adipos', 'risk_groups_perc'='comorb', 
            'pre_life_events'='pre_le', 'pre_contextual_risk'='pre_cr', 'pre_parental_risk'='pre_pr', 'pre_interpersonal_risk'='pre_ir', 
           'post_life_events'='pos_le','post_contextual_risk'='pos_cr','post_parental_risk'='pos_pr','post_interpersonal_risk'='pos_ir', 
           'post_direct_victimization'='pos_dv')
names(longdat) <- stringr::str_replace_all(names(longdat), rep_str)

# reverse binary moderator variables
longdat$exerc_binR <- ifelse(longdat$exerc_bin == 1, 0, 1) 
longdat$sleep_binR <- ifelse(longdat$sleep_bin == 1, 0, 1)
longdat$mdiet_binR <- ifelse(longdat$mdiet_bin == 1, 0, 1)
# Include a numeric version of the exercise factor, to transform into z score 
longdat$exerc_fac <- longdat$exerc
longdat$exerc <- as.numeric(longdat$exerc)

# Add total ELS score in case it's missing
longdat$ELS <- rowSums(longdat[,c('prenatal_stress','postnatal_stress')])

sample_mids <- as.mids(longdat)

# Add standardized scores for main variables
mods = c('exerc','sleep','mdiet')
domains <- c('pre_le','pre_cr','pre_pr','pre_ir','pos_le','pos_cr','pos_pr','pos_ir','pos_dv')
foodgroups <- paste0('mdiet_', c('veg','leg','fru','cer','fis','mea','dai'))

to_trans <- c('ELS', mods, domains, foodgroups)
sample_mids <- datlist2mids( scale_datlist( mids2datlist(sample_mids), 
                             orig_var = to_trans, trafo_var = paste0(to_trans, '_z')))
# Create total lifestyle score
longdat <- complete(sample_mids,'long', include=T)
longdat$lifes <- rowSums(longdat[, paste0(mods,'_z')])
sample_mids <- as.mids(longdat)

# Save finished mids 
saveRDS(sample_mids, file.path(datapath,'imp_sample.rds'))
# Save each imputed setfor portenial plottining
for (m in 0:30) { write.csv(complete(sample_mids,m), file.path(datapath,'byimp',paste0('imp',m,'.csv'))) }

# ------------------------------------------------------------------------------

pool_fit <- function(mods=c('exerc','sleep','mdiet'), z_or_bin='_z', interact='*',
                     outc='comorb', exp='ELS_z') {
  cat(outc, '\n')
  # Define covariates
  covs <- '+ age_child + sex + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking'
  # Initiate stack
  models = data.frame()
  
  for (modr in paste0(mods,z_or_bin)) {
    # IF outcome is comorbidity group fit a multimod model else lm
    if (startsWith(outc, 'comorb')) { 
      fit <- with(sample_mids, nnet::multinom(as.formula(paste(outc,'~',exp, interact, modr, covs)), 
                                              model=T, trace=F));
    } else { 
      fit <- with(sample_mids, lm(as.formula(paste(outc,'~',exp, interact, modr, covs)))); }
    
    p_fit <- mice::pool(fit) # pool results 
    mod <- summary(p_fit) # extract relevant information
    mod[,-c(1,2)] <- round(mod[,-c(1,2)],4)
    mod$sign <- ifelse(mod$p.value < 0.05, '*', '') # add a column to highlight significant terms
    
    cat('--------', as.character(mod[nrow(mod),'term']), '--> p =', round(mod$p.value[nrow(mod)], 3), '\n')
    
    if (startsWith(outc, 'comorb')) {
      # make group comparisons easier to read
      if (endsWith(outc, 'REC')) { levels(mod$y.level) <- c("C:healthy","C:intern", "C:adipos")
      } else { levels(mod$y.level) <- c("H:intern", "H:adipos", "H:comorb") }
      mod$OR  <- round(exp(mod$estimate), 4)
      mod$lci <- round(exp((mod$estimate) - 1.96*mod$std.error), 4)
      mod$uci <- round(exp((mod$estimate) + 1.96*mod$std.error), 4)
      mod$AIC <- c(mean(p_fit$glanced$AIC), rep(NA, nrow(mod)-1)) # add a column for AIC
    } else {
      mod$lci <- round((mod$estimate - 1.96*mod$std.error), 4)
      mod$uci <- round((mod$estimate + 1.96*mod$std.error), 4)
      mod$rsq <- round(c(pool.r.squared(fit)[1], rep(NA, nrow(mod)-1)),4) # add a column for R2
      mod$rsq_adj <- round(c(pool.r.squared(fit, adjusted = TRUE)[1], rep(NA, nrow(mod)-1)),4) # adjusted R2
    }
    mod <- cbind(rep(modr, nrow(mod)),mod)
    
    models <- rbind(models, mod, rep(NA, ncol(mod)))
  }
  names(models)[1] = 'model' # print(mod)
  return(models)
}
# ------------------------------------------------------------------------------
sink(file.path(respath,'summary_int_pvals.txt'))
# Baseline analyses (comorbidity [vs. healthy], internalizing and adiposity)
els_com <- pool_fit()
els_int <- pool_fit(outc='intern_z')
els_fat <- pool_fit(outc='adipos_z')

# prenatal and postnatal exposure assessment
pre_com <- pool_fit(exp='prenatal_stress_z')
pre_int <- pool_fit(exp='prenatal_stress_z', outc='intern_z')
pre_fat <- pool_fit(exp='prenatal_stress_z', outc='adipos_z')
pos_com <- pool_fit(exp='postnatal_stress_z')
pos_int <- pool_fit(exp='postnatal_stress_z', outc='intern_z')
pos_fat <- pool_fit(exp='postnatal_stress_z', outc='adipos_z')

# EXTRA ------------------------------------------------------------------------

# Total lifestyle score 
ls_els_com <- pool_fit(mods=c('lifes'), z_or_bin='')
ls_els_int <- pool_fit(mods=c('lifes'), z_or_bin='', outc='intern_z')
ls_els_fat <- pool_fit(mods=c('lifes'), z_or_bin='', outc='adipos_z')

ls_pre_com <- pool_fit(mods=c('lifes'), z_or_bin='', exp='prenatal_stress_z')
ls_pre_int <- pool_fit(mods=c('lifes'), z_or_bin='', exp='prenatal_stress_z', outc='intern_z')
ls_pre_fat <- pool_fit(mods=c('lifes'), z_or_bin='', exp='prenatal_stress_z', outc='adipos_z')
ls_pos_com <- pool_fit(mods=c('lifes'), z_or_bin='', exp='postnatal_stress_z')
ls_pos_int <- pool_fit(mods=c('lifes'), z_or_bin='', exp='postnatal_stress_z', outc='intern_z')
ls_pos_fat <- pool_fit(mods=c('lifes'), z_or_bin='', exp='postnatal_stress_z', outc='adipos_z')

# Individual stress domains 
for (d in domains) {
  assign(paste0(d,'_com'), pool_fit(exp=paste0(d,'_z')))
  assign(paste0(d,'_int'), pool_fit(exp=paste0(d,'_z'), outc='intern_z'))
  assign(paste0(d,'_fat'), pool_fit(exp=paste0(d,'_z'), outc='adipos_z'))
  assign(paste0('ls_',d,'_com'), pool_fit(exp=paste0(d,'_z'), mods=c('lifes'), z_or_bin=''))
  assign(paste0('ls_',d,'_int'), pool_fit(exp=paste0(d,'_z'), mods=c('lifes'), z_or_bin='', outc='intern_z'))
  assign(paste0('ls_',d,'_fat'), pool_fit(exp=paste0(d,'_z'), mods=c('lifes'), z_or_bin='', outc='adipos_z'))
}

# Food groups
fg_els_com <- pool_fit(mods=foodgroups)
fg_els_int <- pool_fit(mods=foodgroups, outc='intern_z')
fg_els_fat <- pool_fit(mods=foodgroups, outc='adipos_z')

# comorbid as the reference group
els_comR <- pool_fit(outc='comorb_REC')

# Binary moderator (comorbidity [vs. healthy], internalizing and adiposity)
b1_els_com <- pool_fit(z_or_bin='_bin')
b1_els_int <- pool_fit(z_or_bin='_bin', outc='intern_z')
b1_els_fat <- pool_fit(z_or_bin='_bin', outc='adipos_z')
b2_els_com <- pool_fit(z_or_bin='_binR') # reversed binary moderator
b2_els_int <- pool_fit(z_or_bin='_binR',outc='intern_z')
b2_els_fat <- pool_fit(z_or_bin='_binR',outc='adipos_z')

sink()

# No interaction models 
ad_els_com <- pool_fit(interact='+')
ad_els_int <- pool_fit(interact='+', outc='intern_z')
ad_els_fat <- pool_fit(interact='+', outc='adipos_z')


# ------------------------------------------------------------------------------
names = ls()[grepl('int|fat|com', ls())]
for (n in names) { cat("'",n,"'= ",n,", ", sep='')}

modls <- list('els_com'= els_com, 'els_fat'= els_fat, 'els_int'= els_int, 
              'pre_com'= pre_com, 'pre_fat'= pre_fat, 'pre_int'= pre_int, 
              'pos_com'= pos_com, 'pos_fat'= pos_fat, 'pos_int'= pos_int, 
   
   'b1_els_com'= b1_els_com, 'b1_els_fat'= b1_els_fat, 'b1_els_int'= b1_els_int, 
   'b2_els_com'= b2_els_com, 'b2_els_fat'= b2_els_fat, 'b2_els_int'= b2_els_int, 
              
   'els_comR'= els_comR,
              
   'ad_els_com'= ad_els_com, 'ad_els_fat'= ad_els_fat, 'ad_els_int'= ad_els_int)

openxlsx::write.xlsx(modls, file = file.path(respath, paste0('ALSPAC_Results_',date,'.xlsx')), 
                     overwrite=T)
# ------------------------------------------------------------------------------
# names = ls()[grepl('_com|_int|_fat', ls())]
# for (n in names) { cat("'",n,"'= ",n,",\n", sep='')}
suppl <- list(
  'ls_els_com'= ls_els_com, 'ls_els_fat'= ls_els_fat, 'ls_els_int'= ls_els_int,
  'ls_pre_com'= ls_pre_com, 'ls_pre_fat'= ls_pre_fat, 'ls_pre_int'= ls_pre_int, 
  'ls_pos_com'= ls_pos_com, 'ls_pos_fat'= ls_pos_fat, 'ls_pos_int'= ls_pos_int,
  
  'pre_le_com'= pre_le_com, 'pre_le_fat'= pre_le_fat, 'pre_le_int'= pre_le_int, 
  'pre_cr_com'= pre_cr_com, 'pre_cr_fat'= pre_cr_fat, 'pre_cr_int'= pre_cr_int,
  'pre_pr_com'= pre_pr_com, 'pre_pr_fat'= pre_pr_fat, 'pre_pr_int'= pre_pr_int,
  'pre_ir_com'= pre_ir_com, 'pre_ir_fat'= pre_ir_fat, 'pre_ir_int'= pre_ir_int, 
  'pos_le_com'= pos_le_com, 'pos_le_fat'= pos_le_fat, 'pos_le_int'= pos_le_int, 
  'pos_cr_com'= pos_cr_com, 'pos_cr_fat'= pos_cr_fat, 'pos_cr_int'= pos_cr_int,
  'pos_pr_com'= pos_pr_com, 'pos_pr_fat'= pos_pr_fat, 'pos_pr_int'= pos_pr_int, 
  'pos_ir_com'= pos_ir_com, 'pos_ir_fat'= pos_ir_fat, 'pos_ir_int'= pos_ir_int, 
  'pos_dv_com'= pos_dv_com, 'pos_dv_fat'= pos_dv_fat, 'pos_dv_int'= pos_dv_int, 
  
  'ls_pre_le_com'= ls_pre_le_com, 'ls_pre_le_fat'= ls_pre_le_fat, 'ls_pre_le_int'= ls_pre_le_int, 
  'ls_pre_cr_com'= ls_pre_cr_com, 'ls_pre_cr_fat'= ls_pre_cr_fat, 'ls_pre_cr_int'= ls_pre_cr_int, 
  'ls_pre_pr_com'= ls_pre_pr_com, 'ls_pre_pr_fat'= ls_pre_pr_fat, 'ls_pre_pr_int'= ls_pre_pr_int, 
  'ls_pre_ir_com'= ls_pre_ir_com, 'ls_pre_ir_fat'= ls_pre_ir_fat, 'ls_pre_ir_int'= ls_pre_ir_int,
  'ls_pos_le_com'= ls_pos_le_com, 'ls_pos_le_fat'= ls_pos_le_fat, 'ls_pos_le_int'= ls_pos_le_int, 
  'ls_pos_cr_com'= ls_pos_cr_com, 'ls_pos_cr_fat'= ls_pos_cr_fat, 'ls_pos_cr_int'= ls_pos_cr_int, 
  'ls_pos_pr_com'= ls_pos_pr_com, 'ls_pos_pr_fat'= ls_pos_pr_fat, 'ls_pos_pr_int'= ls_pos_pr_int,
  'ls_pos_ir_com'= ls_pos_ir_com, 'ls_pos_ir_fat'= ls_pos_ir_fat, 'ls_pos_ir_int'= ls_pos_ir_int, 
  'ls_pos_dv_com'= ls_pos_dv_com, 'ls_pos_dv_fat'= ls_pos_dv_fat, 'ls_pos_dv_int'= ls_pos_dv_int, 
  
  'fg_els_com'= fg_els_com, 'fg_els_fat'= fg_els_fat, 'fg_els_int'= fg_els_int)

openxlsx::write.xlsx(suppl, file = file.path(respath, paste0('ALSPAC_SupplementaryResults_',date,'.xlsx')), 
                     overwrite=T)

# ------------------------------------------------------------------------------
# NON-LINEARITY ----------------------------------------------------------------
# ------------------------------------------------------------------------------
fit_spline <- function(outc) {
  # Initiate 
  spl_list <- list()
  for (mod in c('exerc','sleep','mdiet')) {
    # Data placeholder
    d = complete(sample_mids, 1)
    # Create grid (min to max)
    lim = range(d[,mod])
    mod.grid = seq(lim[1], lim[2])
    nd = list(mod.grid); names(nd) = mod
    # Initiate plot
    plot(d[,mod], d[,outc], main = 'Splines', xlab=mod, ylab=outc)
    
    # Initiate dataframe
    spl <- data.frame(row.names = mod.grid)
    
    for (i in 1:30) { 
      imp_i <- complete(sample_mids, i)
      c
      pred = predict(fit, newdata = nd, se=T)
      spl <- cbind(spl, pred$fit)
      lines(mod.grid, pred$fit, col='blue',lwd=1)
      # lines(mod.grid, pred$fit+2*pred$se.fit,lty='dashed',lwd=2)
      # lines(mod.grid, pred$fit-2*pred$se.fit,lty='dashed',lwd=2)
    }
    spl_list[[paste0(substr(outc,1,6),'_',mod)]] <- spl
  }
  return(spl_list)
}

nln_int = fit_spline('intern_z')
nln_fat = fit_spline('adipos_z')

openxlsx::write.xlsx(c(nln_int, nln_fat), 
                     file = file.path(respath, paste0('ALSPAC_splines_',date,'.xlsx')), 
                     overwrite=T)

# Comorbidity outcome ----------------------------------------------------------

covs <- '+ age_child + m_bmi_before_pregnancy + m_smoking + m_drinking' # + sex + ethnicity (factors)

for (mod in c('sleep','mdiet','exerc')) { # 'exerc',
  for (exp in c('ELS_z')) { # 'prenatal_stress_z','postnatal_stress_z'
    # Determine output name and initiate 
    out_name = paste0(substr(exp,1,3),'_',mod,'_pps')
    probpred <- data.frame()
    # loop through imputed datasets 
    for (i in 1:30) {
      message('==== ', mod, ' * ', exp,' ==== imp: ', i,' ====================================')
      imp_i = complete(sample_mids, i)
      fit = nnet::multinom(as.formula(paste0('comorb ~',exp,'*',mod,covs)), data=imp_i, Hess=T, 
                           model=T,trace=F)
      
      lim = range(imp_i[,mod], na.rm=T)
      mod.grid = seq(round(lim[1]), round(lim[2]), 2)
      # loop through moderator levels 
      for (l in mod.grid) {
        pred <- MNLpred::mnl_fd_ova(model = fit, data = imp_i,
                                    x = exp,
                                    by = 0.5,
                                    z = mod,
                                    z_values=c(l,l+1),
                                    seed = 310896, # default
                                    nsim = 1000) # faster
        # Only keep comorbid 
        pps <- pred$plotdata[pred$plotdata$comorb=='multimorbid',]
        pps$.imp <- i
        probpred <- rbind(probpred, pps)
      }
    }
    # Average across imputations 
    output = probpred[probpred$.imp==1, c(exp,'comorb',mod)]
    impmean = impmin = impmax = rep(0, nrow(output))
    for (i in 1:30){
      impres = probpred[probpred$.imp==i, c('mean','lower','upper')]
      impmean = impmean + impres[,'mean']
      impmin =  impmin +  impres[,'lower']
      impmax =  impmax +  impres[,'upper']
    }
    output$mean = impmean/30
    output$lower = impmin/30
    output$upper = impmax/30
    
    assign(out_name, output)
  }
} 

# names = ls()[grepl('_pps', ls())]
# for (n in names) { cat("'",n,"'= ",n,", ", sep='')}

pps <- list('els_exerc_pps'= ELS_exerc_pps, 'els_mdiet_pps'= ELS_mdiet_pps, 'els_sleep_pps'= ELS_sleep_pps, 
            'pos_exerc_pps'= pos_exerc_pps, 'pos_mdiet_pps'= pos_mdiet_pps, 'pos_sleep_pps'= pos_sleep_pps, 
            'pre_exerc_pps'= pre_exerc_pps, 'pre_mdiet_pps'= pre_mdiet_pps, 'pre_sleep_pps'= pre_sleep_pps)

pps <- list('els_exerc_pps'= ELS_exerc_pps, 'els_mdiet_pps'= ELS_mdiet_pps, 'els_sleep_pps'= ELS_sleep_pps)
openxlsx::write.xlsx(pps,
                     file = file.path(respath, paste0('ALSPAC_predprobs_',date,'.xlsx')), 
                     overwrite=T)

# ========================================
library(reghelper)

covs <- '+ sex + ethnicity + age_child + m_bmi_before_pregnancy + m_smoking + m_drinking'

for (mod in c('exerc','sleep','mdiet')) {
  for (out in c('intern_z','adipos_z')) { 
    exp='ELS_z' # 'prenatal_stress_z','postnatal_stress_z'
    # Determine output name and initiate 
    out_name = paste0(substr(out,1,3),'_',mod,'_ssl')
    simplesl <- list()
    # loop through imputed datasets 
    for (i in 1:30) {
      message('==== ', mod, ' * ', exp,' ==== imp: ', i,' ====================================')
      imp_i = complete(sample_mids, i)
      imp_i$sleep = round(imp_i$sleep)
      
      fit = lm(as.formula(paste0(out, '~',exp,'*',mod,covs)), data=imp_i)
      
      lvls = list(levels(as.factor(imp_i[,mod])))
      names(lvls) = mod
      ssl = reghelper::simple_slopes(fit, confint=T, levels=lvls)
      
      ssl$.imp <- i
      simplesl[[i]] <- sapply( ssl[, -1], as.numeric )
    }
    # Average across imputations 
    add <- function(x) Reduce("+", x)
    
    output = data.frame(round(add(simplesl)/30, 4))
    
    names(output) = c(mod, 'mean','se','lower','upper','t','df','P')
    
    assign(out_name, output[, -ncol(output)]) 
  }
}

names = ls()[grepl('_ssl', ls())]
for (n in names) { cat("'",n,"'= ",n,", ", sep='')}

ssls <- list('adi_exerc_ssl'= adi_exerc_ssl, 'adi_mdiet_ssl'= adi_mdiet_ssl, 'adi_sleep_ssl'= adi_sleep_ssl, 
             'int_exerc_ssl'= int_exerc_ssl, 'int_mdiet_ssl'= int_mdiet_ssl, 'int_sleep_ssl'= int_sleep_ssl)
openxlsx::write.xlsx(ssls,
                     file = file.path(respath, paste0('ALSP_simpslopes_',date,'.xlsx')), 
                     overwrite=T)
