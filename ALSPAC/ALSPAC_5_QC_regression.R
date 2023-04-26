# Required packages
invisible(lapply(c('ggplot2','gridExtra','splines','car','gvlma','mice','miceadds'),
                 require, character.only = T));

# Data location 
if(!exists('datapath')) { 
  datapath <- dirname(file.choose()) # volumes = X:/Psychology/ResearchProjects/EWalton/ (EarlyCause/WP3/SD/Data)
  reqcpath <- file.path(dirname(datapath),'QC') } 

date <- format(Sys.Date(), "%d%m%y")

# Load data
data <- complete(readRDS(file.path(datapath, 'imp_sample_merged_220223.rds')), 1)
data$ELS <- rowSums(data[,c('prenatal_stress','postnatal_stress')])
data$exercise_fac <- data$exercise
data$exercise <- as.numeric(data$exercise)
for (v in c('ELS','exercise','sleep_hr','med_diet')) { data[,paste0(v,'_z')] <- scale(data[,v])}

# Define main exposures, moderators and outcomes
expo <- c('ELS_z','prenatal_stress_z','postnatal_stress_z')
outc <- c('intern_score_13_z','tot_fat_percent_13_z')
mode <- c('exercise_z','med_diet_z','sleep_hr_z')
covs <- '+ age_child + sex + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking'

# Quickly fit the lm() with relevant variables 
form <- function(outc, exp, mod) {
  fit <- lm(as.formula(paste0(outc,' ~ ',exp,'*',mod,covs)), data=data)
  print(summary(fit))
  return(fit)
}

form('intern_score_13_z','prenatal_stress_z','sleep_hr')
data <- complete(readRDS(file.path(datapath, 'imp_sample_merged_220223.rds')), 2)

mod = 'med_diet'
outc = 'intern_score_13_z'

lim = range(data[,mod])
mod.grid = seq(lim[1],lim[2])
nd = list(mod=mod.grid); names(nd)=mod
fit = lm(as.formula(paste0(outc,' ~ splines::ns(',mod,', 5)')), data=data)
pred = predict(fit, newdata = nd, se=T)

plot(data[,mod], data[,outc], main = 'Try')
lines(mod.grid, pred$fit, col='blue',lwd=3)
lines(mod.grid, pred$fit+2*pred$se.fit,lty='dashed',lwd=2)
lines(mod.grid, pred$fit-2*pred$se.fit,lty='dashed',lwd=2)

sink(file.path(dirname(datapath), paste0('ALSPAC_regress_raw_',date,'.txt')))

method = 'poly'
pdf(file.path(reqcpath,paste0('ALSPAC_nonlin_graphs_',method,'_',date,'.pdf')), width=25, height=5)
for (e in expo) {
  for (o in outc) {
    for (m in mode) {
      cat('---------------------------------------------------\n', o, '~',e, '*',m, 
          '\n---------------------------------------------------\n')
     fit <- form(o,e,m)
      
      p <- list() # create empty list to store p-values
      # plot linear relation 
      p[[1]] <- ggplot(data, aes_string(m,o) ) + geom_point() +
        ggtitle(paste(toupper(substr(o,1,3)),'~',toupper(substr(m,1,3)), '--> linear model')) +
        labs(y = toupper(substr(o,1,3)), x = toupper(substr(m,1,3))) +
        stat_smooth(method = lm, formula = y ~ x)
      
      # fit polynomials 2-5
      for (test in c(2:5)) {
        if (method =='poly') {
          fit2 <- form(o, e, paste0('poly(',m,', ',test,', raw=T)'))
          line = formula(paste0('y ~ poly(x,',test,', raw = T)'))
        } else if (method=='expont') {
          fit2 <- form(o, e, paste0(m,' + I(',m,'^',test,')'))
          line = formula(paste0('y ~ I(x^',test,')'))
        } else if (method=='splines') {
          fit2 <- form(o, e, paste0('splines::ns(',m,', ',test,')'))
          line = formula(paste0('y ~ ns(x,',test,')'))
        }
        # LRT significance test (compared to linear)
        a <- anova(fit, fit2); sign <- ''
        if (!is.na(a[2,'Pr(>F)']) & a[2,'Pr(>F)'] < 0.05) {
          cat(test, '-', e, '~', o, '--------------------------------------\n')
          sign <- '*'
          print(summary(fit2)); print(a)
        }
        pval <- paste(round(a[2,'Pr(>F)'], 3), sign)
        
        p[[test]] <- ggplot(data, aes_string(m,o) ) + geom_point() +
          ggtitle(paste(toupper(substr(o,1,3)),'~',toupper(substr(m,1,3)), '-->', method,'=',test,', p =', pval)) +
          labs(y = toupper(substr(o,1,3)), x = toupper(substr(m,1,3))) +
          stat_smooth(method = lm, formula = as.formula(line))
      }  
      do.call(grid.arrange, c(p, ncol=5))
    }
  }
}
dev.off()
sink()

sample$risk_groups <- relevel(sample$risk_groups, ref = 'healthy')

sink(file.path(dirname(datapath), paste0('ALSPAC_logregress_raw_',date,'.txt')))
for (e in expo) {
  for (m in mode) {
    f <- as.formula(paste('risk_groups ~', e,'*',m, covs))
    cat('---------------------------------------------------\n',e, '*',m, 
        '\n---------------------------------------------------\n')
    model <- nnet::multinom(f, data=sample)
    print(summary(model))
    
    z <- summary(model)$coefficients/summary(model)$standard.errors
    # 2-tailed Wald z tests to test significance of coefficients
    p <- (1 - pnorm(abs(z), 0, 1)) * 2
    print(round(p,3))
  }
}
sink()


# ==============================================================================
# ======================= 2. Regression diagnostics ============================
# ==============================================================================



# Save output to txt log file 
sink(file.path(respath,'QC-regression',paste0('ALSPAC_diagnostics_',date,'.txt')))

# ==============================================================================
# Remember that distinct problems can interact: if the errors have a skewed distribution, 
# apparent outliers may be produced in the direction of the skew. Transforming the 
# response to make the errors less skewed can solve this problem. Similarly, properly 
# modeling a nonlinear relationship may bring apparently outlying observations in line 
# with the rest of the data.
# ==============================================================================

pdf(file.path(respath,'QC-regression',paste0('Diagnostic_plots_',date,'.pdf')))
for (e in exposures) {
  for (o in outcomes) {
    # if (any(sapply(c('gmv','tbv'), grepl, o))) { data <- smri 
    # } else if (any(sapply(c('fa','md'), grepl, o, ignore.case=T))) { data <- dti } 
    data <- full
    
    fit <- form(o, e, data)
    title = paste(toupper(substr(e,1,3)),'~',toupper(substr(o,1,3)))
    cat('\n====',title,'===================================================\n')
    print(gvlma::gvlma(fit))
    cat('\n --- General fit -------------------------------\n')
    # Systematic features in residual plots, e.g. curvature or apparent non-constant variance
    # require to modify the structure of the model to match the data more closely. 
    car::residualPlots(fit, main = title) 
    # Useful for revealing problems but less for determining the nature of the problem. 
    # These should be null plots, with no systematic features. One non-null plot (curved general trend) 
    # is enough to suggest that the specified model does not match the data, generally implying a failure 
    # of one or more assumptions.
    # A lack-of-fit test is also computed for each numeric predictor (i.e. a t-test for (regressor)^2 added to the model. 
    # Significant p-values indicate lack-of-fit, confirming the nonlinear pattern visible in the graph. 
    # For the plot of residuals vs fitted values, the Tukey’s test for non-additivity is obtained by adding the squares 
    # of the fitted values to the model and refitting. The test confirms the impression of curvature (model is not adequate).
    car::marginalModelPlots(fit, main = title) # sd=T to check variance assumptions
    # A lowess smooth is fit to the data (in blue) and to the fitted values (red dashed)
    # and the two curves should match on each of the plots. Any mismatches are evidence of bad fit.
    cat('\n --- Outliers & Influential observations ------\n')
    # Generic qqPlot: Studentized residuals against corresponding quantiles of t(n-k—2),
    # with 95% pointwise confidence envelope. Also returns observations with the largest absolute Studentized residuals.
    # Note qqplots don't work well with the form function
    car::qqPlot(lm(as.formula(paste(o,'~',e,covs1)), data=data), main=title, ylab='Studentized Redisuals')
    # Bonferroni p-value for most extreme obs (i.e. highest Studentized residuals)
    print(car::outlierTest(fit))
    # index plots of Cook’s distances, Studentized residuals, corresponding 
    # Bonferroni p values for outlier testing, and the hat-values (measure of leverage).
    car::influenceIndexPlot(fit, main = title)
    car::influencePlot(fit, id.method="identify", main=title, sub="Circle size proportial to Cook's Distance" )
    cat('\n --- Collinearity -----------------------------\n')
    # print(vif(fit)) # variance inflation factors
    print(sqrt(vif(fit)) > 2) # problem?
  }
}
dev.off()

# ==============================================================================
cat('\n=============================================================
     \n===================== NON LINEARITY =========================
     \n=============================================================\n')
test_nlin <- function(method){
  pdf(file.path(respath,'QC',paste0('Nonlin_graphs_',method,'_',date,'.pdf')), width=25, height=5)
  for (e in expo) {
    for (o in outc) {
      # Define datasets
      # if (any(sapply(c('gmv','tbv'), grepl, o))) { data <- smri 
      # } else if (any(sapply(c('fa','md'), grepl, o, ignore.case=T))) { data <- dti } 
      data <- full
      # fit standard (linear) model
      fit <- form(o, e, data)
      
      p <- list() # create empty list to store p-values
      # plot linear relation 
      p[[1]] <- ggplot(data, aes_string(e,o) ) + geom_point() +
        ggtitle(paste(toupper(substr(e,1,3)),'~',toupper(substr(o,1,3)), '--> linear model')) +
        labs(y = toupper(substr(o,1,3)), x = toupper(substr(e,1,3))) +
        stat_smooth(method = lm, formula = y ~ x)
      
      # fit polynomials 2-5
      for (test in c(2:5)) {
        if (method =='poly') {
          fit2 <- form(o, paste0('poly(',e,', ',test,', raw=T)'), data)
          line = formula(paste0('y ~ poly(x,',test,', raw = T)'))
        } else if (method=='expont') {
          fit2 <- form(o, paste0(e,' + I(',e,'^',test,')'), data)
          line = formula(paste0('y ~ I(x^',test,')'))
        } else if (method=='splines') {
          fit2 <- form(o, paste0('splines::ns(',e,', ',test,')'), data)
          line = formula(paste0('y ~ ns(x,',test,')'))
        }
        # LRT significance test (compared to linear)
        a <- anova(fit, fit2); sign <- ''
        if (!is.na(a[2,'Pr(>F)']) & a[2,'Pr(>F)'] < 0.05) {
          cat(test, '-', e, '~', o, '--------------------------------------\n')
          sign <- '*'
          print(summary(fit2)); print(a)
        }
        pval <- paste(round(a[2,'Pr(>F)'], 3), sign)
        
        p[[test]] <- ggplot(data, aes_string(e,o) ) + geom_point() +
          ggtitle(paste(toupper(substr(e,1,3)),'~',toupper(substr(o,1,3)), '-->', method,'=',test,', p =', pval)) +
          labs(y = toupper(substr(o,1,3)), x = toupper(substr(e,1,3))) +
          stat_smooth(method = lm, formula = as.formula(line))
      }  
      do.call(grid.arrange, c(p, ncol=5))
    }
  }
  dev.off()
}

test_nlin('poly')
test_nlin('expont')
test_nlin('splines')

sink()


