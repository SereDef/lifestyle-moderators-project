library(openxlsx)

# lil help
s <- function(var) {summary(var)}; sf <- function(var) {summary(as.factor(var))}

# Data location 
if(!exists('datapath')) { datapath <- dirname(file.choose()) } # volumes = X:/Psychology/ResearchProjects/EWalton/ (EarlyCause/WP3/SD/Data)

data <- readRDS(file.path(datapath, 'GenR_data_raw.rds'))
data <- data[,-grep('IDM',names(data))]
imput <- readRDS(file.path(datapath, 'imp_sample.rds'))

descpath <- file.path(dirname(datapath),'Descriptives')

date <- format(Sys.Date(), "%d%m%y")

# Save histograms of all variables in the full sample
pdf(file.path(descpath, paste0('GENR_0_hist_original_sample_',date,'.pdf')))
for (var in names(data)[-which(names(data) %in% c('IDC','sleep_bedtime', 'sleep_waketime'))]) { # these cannot be converted to numeric
  hist(as.numeric(data[,var]), col = 'blue', main = toupper(var), xlab = var) }
dev.off() 

# Sample selection =============================================================
sink(file.path(descpath, paste0('GENR_0_flowchart_selection_',date,'.txt')))
samp1 <- data[data$pre_percent_missing<50,]; 
cat('Excluding participants with insufficient prenatal information (<50%):\n', nrow(data),'-',nrow(data)-nrow(samp1),'=',nrow(samp1), '\n\n')
samp2 <- samp1[samp1$post_percent_missing<50,]; 
cat('Excluding participants with insufficient postnatal information (<50%):\n', nrow(samp1),'-',nrow(samp1)-nrow(samp2),'=',nrow(samp2), '\n\n')
samp3 <- samp2[samp2$twin==0,]; 
cat('Excluding twins:\n', nrow(samp2),'-',nrow(samp2)-nrow(samp3),'=',nrow(samp3), '\n\n')
sample <- samp3[samp3$IDC %in% imput$data$IDC,];
cat('Excluding siblings:\n', nrow(samp3),'-',nrow(samp3)-nrow(sample),'=',nrow(sample),'\n\n')

complete <- sample[!is.na(sample$intern_score_13) & !is.na(sample$tot_fat_percent_13),]
cat('\nNOTE: Complete outcome data:\n', nrow(sample),'-',nrow(sample)-nrow(complete),'=',nrow(complete),'\n\n')

cat('===== Sensitivity analysis - selected sample ================================\n')

for (var in names(data)[-which(names(data)%in%c('IDC','sibling','twin'))]) {
  # NOTE: I do not treat binary variables as factors for this one, change length condition to modify
  if (length(levels(as.factor(data[,var])))<2) { # this is a factor, need to compare proportions
    fac_data <- as.factor(data[,var])
    fac_samp <- as.factor(sample[,var])
    for (l in levels(fac_data)) {
      if (l != 'NaN') {
        t = prop.test(x = c(as.numeric(table(fac_data)[l]), as.numeric(table(fac_samp)[l])), 
                      n = c(nrow(data), nrow(sample)))
        if (t$estimate[1] > t$estimate[2]) { rel = 'full > selected' } else { rel = 'full < selected' }
        if (t$p.value < 0.05) { cat('\n', var,'-', l,'=',rel, 'P =', format(round(t$p.value,3),nsmall=3)) }
      }
    }
    cat('\n')
  } else {
    num_data <- as.numeric(data[,var])
    num_samp <- as.numeric(sample[,var])
    t = t.test(num_data, num_samp, var.equal=TRUE)
    if (t$estimate[1] > t$estimate[2]) { rel = 'full > selected' } else { rel = 'full < selected' }
    if (t$p.value < 0.05) { cat('\n', var, '=',rel, 'P =', format(round(t$p.value,3),nsmall=3), '\n') }
  }
}
cat('\n\n===== Sensitivity analysis - complete outcome sample ============================\n')

for (var in names(data)[-which(names(data)%in%c('IDC','sibling','twin','sleep_bedtime', 'sleep_waketime'))]) {
  # NOTE: I do not treat binary variables as factors for this one, change length condition to modify
  num_data <- as.numeric(data[,var])
  num_samp <- as.numeric(complete[,var])
  t = t.test(num_data, num_samp, var.equal=TRUE)
  if (t$estimate[1] > t$estimate[2]) { rel = 'full > selected' } else { rel = 'full < selected' }
  if (!is.na(t$p.value) & t$p.value < 0.05) { cat('\n', var, '=',rel, 'P =', format(round(t$p.value,3),nsmall=3), '\n') }
}

sink()

# Save sample as file
saveRDS(sample, file.path(datapath, 'GENR_raw_data_sample.rds'))

# Save histograms of all variables in the selected sample
pdf(file.path(descpath, paste0('GENR_0_hist_selected_sample_',date,'.pdf')))
for (var in names(data)[-which(names(data) %in% c('IDC','sleep_bedtime', 'sleep_waketime'))]) { # these cannot be converted to numeric
  hist(as.numeric(data[,var]), col = 'lightblue', main = toupper(var), xlab = var) }
dev.off() 

# ==============================================================================

# Calculate the percentage missing data ----------------------------------------

miss <- data.frame( names(data),
  paste0(colSums(is.na(data)),' (',round((colSums(is.na(data))/nrow(data))*100, 1), '%)'),
  paste0(colSums(is.na(sample)),' (',round((colSums(is.na(sample))/nrow(sample))*100, 1), '%)'),
  paste0(colSums(is.na(complete)),' (',round((colSums(is.na(complete))/nrow(complete))*100, 1), '%)'))
names(miss)<- c('Var name', paste0('Full sample, N = ',nrow(data)),
                            paste0('Selected sample, N = ',nrow(sample)),
                            paste0('Complete outcome, N = ',nrow(complete)))
# View(miss)

# ==============================================================================

summ_to_df <- function(var) {
  if (length(levels(as.factor(data[,var])))<5) { # treat this is a factor
    sum_full = as.matrix(summary(as.factor(data[,var])))
    sum_samp = as.matrix(summary(as.factor(sample[,var])))
    sum_comp = as.matrix(summary(as.factor(complete[,var])))
  } else { 
    sum_full = rbind(as.matrix(summary(data[,var])), SD = sd(as.numeric(data[,var]),na.rm=T))
    sum_samp = rbind(as.matrix(summary(sample[,var])), SD = sd(as.numeric(data[,var]),na.rm=T))
    sum_comp = rbind(as.matrix(summary(complete[,var])), SD = sd(as.numeric(data[,var]),na.rm=T))
  }
  # handle different dimentions
  diff1 = rep(NA, (nrow(sum_full) - nrow(sum_samp))+1)
  diff2 = rep(NA, (nrow(sum_full) - nrow(sum_comp))+1)
  filler = rep(NA, nrow(sum_full)+2)
  # Transform to dataframe adding variable name and levels and empty lines around it 
  df = data.frame('value' = c(var, row.names(sum_full), NA), 'full' = c(NA, sum_full, NA), '.' = filler,
                  'sample' = c(NA, sum_samp, diff1), '..' = filler, 'complete_outcome' = c(NA, sum_comp, diff2))
  return(df)
}

summary_df <- data.frame(matrix(nrow=0, ncol=6))
names(summary_df)=c('Value', paste0('Full sample, N = ',nrow(data)),'.',
                             paste0('Selected sample, N = ',nrow(sample)),'..',
                             paste0('Complete outcome, N = ',nrow(complete)))

for (v in names(data)[-which(names(data)%in%c('IDC','IDM'))]) { summary_df = rbind(summary_df, summ_to_df(v)) } 

# write.csv(summary_df, file.path(descpath, paste0('1_summary_pre_imp_',date,'.csv')))

# FULL AND SELECTED SAMPLE (AFTER IMPUTATION ) =================================

describe <- function(imput) {
  # determine categorical and continuous vars 
  lvl_length <- lapply(imput$data, function(var) length(levels(as.factor(var)))) 
  # Cutoff 15 levels: consider it categorical
  cat_vars <- names(which(lvl_length < 5))
  # correct one issue: this has too few levels 
  # cat_vars <- setdiff(cat_vars, c('post_direct_victimization','mdiet'))
  
  # Stack imputed datasets in long format, excluding the original data
  impdat <- mice::complete(imput, action="long", include = F)
  # Set to factors/ numeric when appropiate 
  impdat[,cat_vars] <- lapply(impdat[,cat_vars] , as.factor)
  impdat[, !names(impdat)%in%cat_vars] <- lapply(impdat[,!names(impdat)%in%cat_vars] , as.numeric)
  
  pool_descriptives <- function(implist, column_names, categorical=T) {
    summ <- with(implist, by(implist, .imp, function(x) summary(x[, -c(1, 2)],digits=4))) 
    if (categorical==F) {
      # Pool summary 
      num_pool <- lapply(summ, function(m) matrix(as.numeric(sapply(strsplit(m, ":"), "[[", 2)), nrow = dim(m)[1], ncol=dim(m)[2]))
      pool_mean <- Reduce("+",num_pool)/length(num_pool)
      # Pool SDs
      sds <- with(implist, by(implist, .imp, function(x) round(apply(x[, -c(1, 2)], 2, sd, na.rm = T), 4)))
      pool_sds <- Reduce("+",sds)/length(sds)
      # Bind SDs to other metrics
      summ_df <- data.frame(rbind(pool_mean,pool_sds))
      # Define column and row names
      colnames(summ_df) <- colnames(implist[-c(1,2)])
      rownames(summ_df) <- c('Min','1stQ','Median','Mean','3rdQ','Max','SD')
    } else { 
      pool_mean <- Reduce("+",summ)/length(summ) 
      summ_df <- data.frame(pool_mean)
      colnames(summ_df) <- 'counts' # colnames(implist[-c(1,2)])
      rownames(summ_df) <- names(pool_mean)
    }
    return(summ_df)
  }
  # Continuous data
  cnt <- impdat[, -c(which(colnames(impdat) %in% cat_vars))]
  cnt_summ <- cbind(c('Min','1stQ','Median','Mean','3rdQ','Max','SD'), pool_descriptives(cnt, categorical = F))
  
  # Categorical / binary
  cat_summ <- NA
  for (v in cat_vars) {
    sel <- impdat[, c('.imp','.id', v)]
    v_summ <- pool_descriptives(sel)
    v_summ <- cbind(row.names(v_summ), v_summ)
    v_summ$percent <- (v_summ$count / nrow(imput$data))*100
    cat_summ <- rbind(cat_summ, v, v_summ, NA)
  }
  
  # Correlation matrix in the imputed set
  cors_imp <- miceadds::micombine.cor(mi.res = impdat, 
                                      variables = colnames(impdat)[!colnames(impdat) %in% 
                                                                     c('.imp','.id','IDC',cat_vars)]) 
  
  # Export the outputs of summary statistics into an xlsx file with one model per sheet
  stats <- list( # 's_pre_imp' = summ_sample, 
                's_imp_cnt' = cnt_summ, 's_imp_cat' = cat_summ, 
                'cor_imp' = cors_imp)
  
  return(stats)
}

s = describe(imput)


openxlsx::write.xlsx(c(list('missing_pattern'=miss,'summary'=summary_df), s), 
                     file = file.path(descpath, paste0('GENR_Descriptives_',date,'.xlsx')), 
                     overwrite=T)
