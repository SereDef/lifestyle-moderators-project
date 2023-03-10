
library(foreign)

datapath <- dirname(file.choose()) 
# datapath <- '/Users/Serena/Desktop/moderator_analysis/Data'

# Read in the files ------------------------------------------------------------
readsav <- function(file, summary=F) {
  d <- foreign::read.spss(file.path(datapath, file), use.value.labels=T, to.data.frame=T)
  return(d)
}
# Load datasets - moderators
sport  <- readsav('GR1084_C10-C11_04042017.sav');
sleep  <- readsav('SLEEPCHILDF9F13_SOLWASOTSTDiary-interview_22102021.sav'); # 'GR1096-E1-14_SHS_09072020.sav'
diet   <- readsav('CHILDNUTRITIONFOODGROUPS_26042015.sav');
# Load datasets for imputation 
sport4 <- readsav('GR1067-C1-15_01072012.sav');
sport6 <- readsav('GR1075-D1-25_17072015.sav'); # GR1075-D5_13062016.sav
sport9 <- readsav('GR1081_H5-8_21122020.sav');
cbcl1  <- readsav('20200111_GR1029_CBCL A1_1.5_incl_Tscores.sav');
cbcl3  <- readsav('CBCL_3_incl_Tscores__GR1065E2_GR1066A1_20201111.sav');
cbcl6  <- readsav('CHILDCBCL_6_incl_Tscores_20201111.sav'); # CHILDCBCL5_17072015.sav
# Load outcome and covariate information 
outcs   <- readRDS(file.path(datapath, 'PCM_allvars.rds'))[,c('IDC','age_child_dxa','age_child_cbcl')]


# merge them
full <- Reduce(function(x,y) merge(x = x, y = y, by = 'IDC', all.x = T),
                      list(sport, sleep, diet, sport4, sport6, sport9, cbcl1, cbcl3, cbcl6, outcs) ) 

# Initiate dataframe and set up unique child ID
data <- data.frame('IDC' = full$IDC)

# Tidy up
rm(sport, sleep, diet, sport4, sport5, sport9, cbcl1, cbcl3, cbcl6, outcs)

s <- function(v) { summary(v) }
sf <- function(v) { summary(as.factor(v)) }
# ==============================================================================
# EXERCISE (9.8 years) ========================================================
# ==============================================================================
# Do you play sports at a sports club or team? How often do you play sports?

data$exercise <- as.factor(ifelse(full[,'C1100184_cleaned'] == 'No' & is.na(full[,'C1100284_cleaned']), 0,
                        full[,'C1100284_cleaned'])) # add 0 instead of NA for no exercise
levels(data$exercise) <- c('none','1 week','2 week','3 week','4 week','>=5 week')
#.                         'none','<1 week','1-3 week','4-6 week','daily'
summary(data$exercise)
  
data$exercise_bin <- as.factor(ifelse(as.numeric(data$exercise) > 4, 1, 0)) # more than 3 times a week
summary(data$exercise_bin)

data$exercise_age <- full[,'ageChildGR1084']

# Exercise auxiliary timepoints ------------------------------------------------

# Repeated measures of the exercise variables are not neat as in ALSPAC. Therefore we use information on: 
# - Frequency of participation in sport per week (at 6 and 9.7 years)
# - Frequency of walking or cycling to school per week (at 6 and 9.7 years)
# - Hours of outdoor play per week (at 4 and 9.7 years)

#### Frequency of participation in sport per week ##############################
# 6 years -------------
data$exercise_sports_6y <- ifelse(!is.na(full[,'playing_sports']) & full[,'playing_sports'] == 'no', 0, # "Does your child play any sport" yes or no
                           ifelse( is.na(full[,'playing_sports']) & rowSums(is.na(full[,c('sp1_days','sp2_days','sp3_days')])) > 2, NA, # missing values
                                   rowSums(full[,c('sp1_days','sp2_days','sp3_days')], na.rm=T))) # sum of hours across 3 sports
data$exercise_sports_6y[data$exercise_sports_6y > 7] <- 7 # Fix two cases coded as 8 times per week

# 10 years -------------
# How many hours per week does your child spend doing sports (training and compete together)?
# Less than 1 hour / 1 to 2 hours / 2 to 4 hours / more than 4 hours per week 
data$exercise_sports_10y <- full[,'H0700281_cleaned']

#### Frequency of walking or cycling to school per week ########################
# 6 years -------------
# On average, how many days per week does your child walk from home to school or back?
# On average, how many days per week does your child cycle him/herself from home to school or back?
wc_6y <- data.frame('walk_6y' = as.numeric(full[,'D0100175'])-1, # (0) Never; (1); (2); (3); (4); (5) days per week
                   'cycle_6y' = as.numeric(full[,'D0200175'])-1)

data$exercise_walk.cycle_6y <- ifelse(rowSums(is.na(wc_6y)) == 2, NA, rowSums(wc_6y, na.rm=T))
data$exercise_walk.cycle_6y[data$exercise_walk.cycle_6y > 5] <- 5

# 10 years -------------
# On average how many days per week does your child walk or cycle from home to school and/or vice/versa?
# (0) Never; (1) 1-2 days; (2) 3-4 days; (3) 5 days per week
levels(full[,'H0500181_cleaned']) <- c(0, 1.5, 3.5, 5)
data$exercise_walk.cycle_10y <- as.numeric(as.character(full[,'H0500181_cleaned']))

####  Hours of outdoor play per week ###########################################
# 4 years -------------
# How many days per week does your child play outside (outside school times) 
# Never or <1 day per week  / 1 day per week / 2 / 3 / 4 / 5 / 6 / 7 days per week
outdoor_days_4y <- as.numeric(full[,'C1400167'])-1
# How long does your child generally play outside each day (outside school times)
# Less than 30 minutes per day / 30 minutes to 1 hour / 1 to 2 hours / 2 to 3 hours / more than 3 hours per day
levels(full[,'C1500167']) <- c(0.15, 0.75, 1.5, 2.5, 3.5)
outdoor_hours_4y <- as.numeric(as.character(full[,'C1500167']))

data$exercise_outdoor_4y <- outdoor_days_4y * outdoor_hours_4y

# 10 years -------------
# On average how many days per week does your child play outside?
# (0) Never; (1) 1-2 days; (2) 3-4 days; (3) 5 or more days per week
levels(full[,'H0800181_cleaned']) <- c(0, 1.5, 3.5, 5.5)
outdoor_days_10y <- as.numeric(as.character(full[,'H0800181_cleaned']))
# Approximately how long does your child approximately play outside per day? Only consider the days that your child plays outside.
# Less than 30 minutes per day / 30 minutes to 1 hour / 1 to 2 hours / 2 to 3 hours / 3 to 4 hours / More than 4 hours per day
levels(full[,'H0800281_cleaned']) <- c(0.15, 0.75, 1.5, 2.5, 3.5, 4.5)
outdoor_hours_10y <- as.numeric(as.character(full[,'H0800281_cleaned']))

data$exercise_outdoor_10y <- outdoor_days_10y * outdoor_hours_10y

# ==============================================================================
# SLEEP (11.8 / ~14.6 years) ===================================================
# ==============================================================================
# Mother reported sleep duration
data$sleep_hr <- full[,'TST_M']
# NOTE: alternative self-reported sleep duration (but lower N)
data$sleep_SR <- rowMeans(full[, paste0('TST', 1:9)], na.rm=T)
summary(data[,c('sleep_hr','sleep_SR')])
# Correlation between self reported and mother reported sleep duration
round(cor(data[,c('sleep_hr','sleep_SR')], use = 'pairwise.complete.obs')[1,2],2)

data$sleep_hr_age <- ifelse(full[,'ageslaapChild']==999, NA, full[,'ageslaapChild'])

# NOTE: data was collected in two waves, round 12 and again around 15 years. 
hist(data$sleep_hr_age, breaks=100)

full$p <- as.factor(ifelse(full$sleep_hr_age > full$age_child_cbcl, 1,0))
# Exclusion of children that are too old (sleep measured after the outcomes)
data$sleep_hr[(data$sleep_hr_age > full$age_child_cbcl | data$sleep_hr_age > full$age_child_dxa)] <- NA
data$sleep_SR[(data$sleep_hr_age > full$age_child_cbcl | data$sleep_hr_age > full$age_child_dxa)] <- NA
summary(data[,c('sleep_hr','sleep_SR')])

# Dichotomize based on recommendations
data$sleep_hr_bin <- as.factor(ifelse(data$sleep_hr < 9 | data$sleep_hr > 12, 0, 1))
data$sleep_SR_bin <- as.factor(ifelse(data$sleep_SR < 9 | data$sleep_SR > 12, 0, 1))
summary(data[,c('sleep_hr_bin','sleep_SR_bin')])

# Sleep auxiliary timepoints ---------------------------------------------------

# Repeated measures of sleep duration are not available. Therefore we use 
# - mother-reported child sleep quality (CBCL – sleep problems scale) at 1.5, 3 and 6 years, 

data$sleep_probs_1.5y <- full[,'sum_sle']
data$sleep_probs_3y   <- full[,'sum_sle_36m'] # sum_sle_36p
data$sleep_probs_6y   <- full[,'sum_sle_5']

# ==============================================================================
# DIET (MEDITERRANEAN) (10.7 years) ============================================
# ==============================================================================

# Nutrients mean weight (g/day) [ FFQ ]
diet <- data.frame( # Extract form main dataset and rename
        # ------ CEREAL ---------------------------------
         'cer_1' = full[,'FG60_04_Bread_whole'], # Brown/wholegrain bread
         'cer_2' = full[,'FG60_06_Grainproducts_whole'], 
         'cer_3' = full[,'FG60_34_Porridge'], 
         # FG60_05_Bread_white / FG60_07_Grainproducts_white
         # FG60_09_Cereals_white_sugar
         # FG60_41_Sweetsnacks # Sweets, sweet snacks, and cookies
         # FG60_42_Savorysnacks # Savory snacks, excl potato crisp
         # FG60_50_Sweettoppings # Sweet toppings or sandwich fillings
         # FG60_51_Potatoes_crisps
         # FG60_52_Pizza_savorypie # Pizza, savory pie, quiche
         # FG60_53_Pastadishes
         # FG60_56_Nasi_bami_dishes # Stir-fried Indonesian rice or noodles
        
        # ------ FISH ---------------------------------
         'fis_1' = full[,'FG60_10_Fish_fat'], # Fatty fish (>10% fat) 
         'fis_2' = full[,'FG60_11_Shellfish'], # Shellfish
         'fis_3' = full[,'FG60_12_Fish_lowfat'], # Lean fish (<2% fat) 
         'fis_4' = full[,'FG60_13_Fish_modfat'], # Moderately fat fish (2-10%)  
         'fis_5' = full[,'FG60_14_Fish_fat_canned'], # Canned fatty fish
         # FG60_08_Fish_breaded
        
        # ------ DAIRY -------------------------------
         'dai_1' = full[,'FG60_24_Milk_lowfat'], # Skimmed or semi skimmed milk, buttermilk, no added sugar
         'dai_2' = full[,'FG60_25_Yogurt_lowfat'], # Skimmed or semi-skimmed yoghurt or quark, no added sugar
         'dai_3' = full[,'FG60_26_Milk_fullfat'], # Full-fat milk
         'dai_4' = full[,'FG60_27_Yoghurt_fullfat'], # Full-fat yoghurt or quark
         'dai_5' = full[,'FG60_28_Milk_beverages_sugar'], # Milk-based beverages with added sugar 
         'dai_6' = full[,'FG60_29_Yogurt_sugar'], # Yoghurt or quark with added sugar
         'dai_7' = full[,'FG60_30_Cheese_lowfat'], # Low-fat cheese (<=30+)
         'dai_8' = full[,'FG60_31_Cheese_fullfat'], # Full-fat cheese (> 30+)
         # FG60_32_Dairy_desserts # Dairy based desserts (incl ice-cream)
         # FG60_33_Soymilks # Soy milk/ soy dessert
        
        # ------ MEAT -------------------------------- 
         'mea_1' = full[,'FG60_15_Meat_red_unprocessed_lowfat'], # Red meat, unprocessed, low-fat (<=5% SFA) 
         'mea_2' = full[,'FG60_16_Meat_red_processed'], # Red meat, processed
         'mea_3' = full[,'FG60_17_Meat_white_processed'], # White meat, processed
         'mea_4' = full[,'FG60_18_Meat_white_unprocessed_lowfat'], # White meat, unprocessed, low-fat (<=5% SFA) 
         'mea_5' = full[,'FG60_35_Fats_hard'], # Hard fats, butter
         'mea_6' = full[,'FG60_45_Meat_red_unprocessed_highfat'], # Red meat, unprocessed, high-fat (>5% SFA) 
         'mea_7' = full[,'FG60_47_Sausagerolls'], # Sausage rolls (worstebroodje, sausijzenbroodje)
         'mea_8' = full[,'FG60_48_Meat_fastfood'], # Fast food meat (kroket frikande)
         # FG60_19_Meatreplacing_products
         # FG60_36_Fats_soft_lowSFA # Soft fats(>=30% saturated fat of total fat), oils
         # FG60_40_Eggs
         # FG60_49_Sauces_fat # Fat- containing sauses
        
        # ------ VEGETABLES -------------------------------- 
         'veg_1' = full[,'FG60_02_Vegetables_cooked'],
         'veg_2' = full[,'FG60_03_Vegetables_raw'],
         # FG60_43_Potatoes / FG60_57_Potatoes_frenchfries
        
        # ------ LEGUMES ----------------------------------
         'leg_1' = full[,'FG60_20_Pulses_canned'], # Pulses, canned
        
        # ------ FRUIT --------------------------------
         'fru_1' = full[,'FG60_01_Fruit_fresh'], # Fresh fruit
         'fru_2' = full[,'FG60_21_Nuts_unsalted'], # Nuts, unsalted
         'fru_3' = full[,'FG60_44_Fruit_dried'] # Dried fruit (raisins) 
         # FG60_37_Fruitjuice / FG60_22_Peanuts_peanutbutter
)

foodgroups <- c('veg','leg','fru','cer','fis','mea','dai')

# Compute a sum score per food group and the corresponding binary variable
for (foodgroup in foodgroups) {
  message('\n', toupper(foodgroup))
  
  name_c <- paste0('med_diet_',foodgroup) # Continuous 
  name_b <- paste0('med_diet_',foodgroup,'_bin') # Binary
  
  # Compute continuous (i.e., sum of food group components)
  if (foodgroup=='leg') { diet[, name_c] <- diet[,'leg_1'] # only one item in the group
  } else { diet[, name_c] <- rowSums(diet[, grepl(foodgroup, names(diet))], na.rm = F) }
  
  print(summary(diet[, name_c]))
  
  # Dichotomize at median intake
  med <- median(diet[, name_c], na.rm = T)
  if (foodgroup %in% c('mea','dai')) { # there are reverse coded (i.e. non beneficial)
    diet[, name_b] <- ifelse(diet[, name_c] < med, 1, 0)
  } else {
    diet[, name_b] <- ifelse(diet[, name_c] >= med, 1, 0)
  }
  message("Dichotomized ", foodgroup, ' at median intake = ', med)
  print(summary(as.factor(diet[, name_b])))
}

# Sum of binary diet components
data$med_diet <- rowSums(diet[, paste0('med_diet_',foodgroups,'_bin')]) 
summary(data$med_diet)

# Dichotomized diet score (at median)
data$med_diet_bin <- as.factor(ifelse(data$med_diet >= median(data$med_diet, na.rm=T), 1, 0))
summary(data$med_diet_bin)

data$med_diet_age <- full[,'agechild_GR1080']

# Include single nutrients weight in the dataset, I also save a separate dataset
data <- cbind(data, diet[, names(diet)[grep('med_diet', names(diet))]])

# No auxiliary variables for diet imputation -----------------------------------
# NOTE: missing dietary components were automatically imputed by SAS vovris, 
# a software package for the processing of FFQs (Dutman et al., 2011)

# ==============================================================================
# merge with the rest of the variables 
pren_stress <- readRDS(file.path(datapath, 'prenatal_stress.rds'))
post_stress <- readRDS(file.path(datapath, 'postnatal_stress.rds'))
outcome_cov <- readRDS(file.path(datapath, 'PCM_allvars.rds')) # PCM_allvars

data_full <- Reduce(function(x,y) merge(x = x, y = y, by = 'IDC', all.x = T),
                    list(pren_stress, post_stress, outcome_cov, data) ) 

data_full$ELS <- rowSums(data_full[,c('prenatal_stress','postnatal_stress')])

# Scale main variables 
for (v in c('prenatal_stress','postnatal_stress','ELS',
            'intern_score_13','tot_fat_percent_13',
            'exercise','sleep_hr','sleep_SR','med_diet')) {
  data_full[,paste0(v,'_z')] <- scale(as.numeric(data_full[,v]))
}
# Construct group outcome
int = as.factor(ifelse(data_full[, 'intern_score_13'] > quantile(data_full[, 'intern_score_13'], probs = 0.8, na.rm = T), 1, 0))
fat = as.factor(ifelse(data_full[,'tot_fat_percent_13'] > quantile(data_full[, 'tot_fat_percent_13'],  probs = 0.8, na.rm = T), 1, 0))
data_full$risk_groups = as.factor(ifelse(int==0 & fat==0, 'healthy',
                                         ifelse(int==1 & fat==0, 'high_intern',
                                                ifelse(int==0 & fat==1, 'high_adipos',
                                                       ifelse(int==1 & fat==1, 'comorbidity', NA)))))

# ==============================================================================
# Save and flee

# Main file
saveRDS(data, file.path(datapath, 'GenR_moderator_vars.rds'))
# Full data file
saveRDS(data_full, file.path(datapath, 'GenR_data_raw.rds'))
# Diet data for additional descriptive / check-ups
saveRDS(diet, file.path(datapath, 'GenR_diet_vars.rds'))
