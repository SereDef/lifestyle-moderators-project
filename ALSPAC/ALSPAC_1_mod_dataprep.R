
library(foreign)

# File name and location 
filename <- 'EarlyCause_AHupdated_CIDB2957_21OCT21.sav' # change if needed
datapath <- dirname(file.choose()) # volumes = X:/Psychology/ResearchProjects/EWalton/ (EarlyCause/WP3/SD/Data)

# Read in the file
full <- foreign::read.spss(file.path(datapath, filename), 
                           use.value.labels=F, to.data.frame=T)
names(full) <- tolower(names(full)) # all column names to lower case

# Initiate dataframe and set up unique child ID
data <- data.frame('IDC' = paste(full$cidb2957, # mother ID
                                 gsub('\\s+', '', full$qlet), # sibling ID (A or B)
                                 sep = "_")) # --> format "1_A"

# ==============================================================================
# EXERCISE (10.7 years) ========================================================
# ==============================================================================
# Average number of times child participated in vigorous physical activity in past month
# (running, dance, gymnastics, netball, swimming or aerobics) # girls
# (running, football, swimming, athletics) # boys

cleanExercise <- function(var) {
  print(summary(full[,var])) # original variable summary
  # Transform to factor and set up answer values
  out <- as.factor(full[,var])
  levels(out) <- c('none','<1 week','1-3 week','4-6 week','daily')
  
  cat('\n'); print(summary(out))
  return(out)
}

data$exercise <- cleanExercise('pub309') 

data$exercise_bin <- as.factor(ifelse(full[,'pub309'] >= 4, 1, 0)) # more than 3 times a week
summary(data$exercise_bin)

data$exercise_age <- full[,'pub397a']/12

# Exercise auxiliary timepoints ------------------------------------------------
# NOTE: the variable name is based on mean age at questionnaire completion
data$exercise_08.2y <- cleanExercise('pub109')
data$exercise_09.7y <- cleanExercise('pub209')
# data$exercise_10.7y <- cleanExercise('pub309') # MAIN VARIABLE!!
data$exercise_11.7y <- cleanExercise('pub409')
data$exercise_13.2y <- cleanExercise('pub509')
data$exercise_14.7y <- cleanExercise('pub609')
data$exercise_15.4y <- cleanExercise('pub709')
data$exercise_16.1y <- cleanExercise('pub809')
data$exercise_17.0y <- cleanExercise('pub909')

# ==============================================================================
# SLEEP (11 years) =============================================================
# ==============================================================================

# Mother report of bed and wake time are combined into sleep duration (in hours)
cleanSleep <- function(bed_var, wak_var, bed_wake_times=F) {
  # Identify bed and wake time variables 
  # - Time child usually goes to sleep on normal school days - hours & minutes
  # - Time child usually wakes up on normal school days - hours & minutes    
  if (length(bed_var) > 1) { # in one case these were not coded as [var]a and [var]b
    bed_h <- full[,bed_var[1]]; bed_m <- full[,bed_var[2]]
    wak_h <- full[,wak_var[1]]; wak_m <- full[,wak_var[2]]
  } else { # in all other cases these are coded as [var]a and [var]b
    bed_h <- full[,paste0(bed_var,'a')]; bed_m <- full[,paste0(bed_var,'b')]
    wak_h <- full[,paste0(wak_var,'a')]; wak_m <- full[,paste0(wak_var,'b')]
  }
  
  # Clean and recode to date-time format
  # NOTE: I add leading 0s for correct time format e.g., 07:05
  bed_h <- ifelse(bed_h < 5 | bed_h > 24, NA, # implausible value (= 16, 25, 26, 97) recoded as NA
           # ifelse(bed_h > 24, paste0(0, bed_h-24), # 2 values coded 25/26, recoded to 01/02
           ifelse(bed_h %in% c(12,24), '24', bed_h+12))#) # NOTE: one case: 24 recoded to '00'
  
  bed_m <- ifelse(is.na(bed_m), '00', # NA recoded as 00 to avoid loosing values in merging
           ifelse(bed_m > 59, NA, # implausible value (=97) recode as NA. 
           ifelse(bed_m < 10, paste0(0,bed_m), bed_m)))
  
  wak_h <- ifelse(wak_h < 4, NA, # implausible values (2 and 3) recoded as NA
           ifelse(wak_h < 10, paste0(0, wak_h), wak_h)) # no recoding required
  wak_m <- ifelse(is.na(wak_m), '00', # NA recoded as 00
           ifelse(wak_m < 10, paste0(0, wak_m), wak_m))
  
  # Paste hour and minutes together. If hour is reported but minutes are NA, consider minutes = 00
  #bedtime <- ifelse(is.na(bed_h), NA, paste(bed_h, bed_m, sep=':'))
  #waketime <- ifelse(is.na(wak_h),NA, paste(wak_h, wak_m, sep=':'))
  
  # Transform into date-time format to compute time difference. 
  # NOTE: I pick two random consecutive days to ensure correct time attribution
  #bed <- strptime(x = ifelse(bed_h<10, paste('2022-08-31',bedtime), paste('2022-08-30',bedtime)), format = "%Y-%m-%d %H:%M")
  #wak <- strptime(x = paste0('2022-08-31 ',waketime), format = "%Y-%m-%d %H:%M")
  
  # Compute time difference between bed and wake time, = continuous sleep variable 
  #sleeptime <- as.vector(difftime(wak, bed, unit='hours'))
  
  # Transform to factors
  #bedtime <- as.factor(bedtime)
  #waketime <- as.factor(waketime)
  
  # Move midnight values at the end for bed times 
  # levels(bedtime) <- c(levels(bedtime)[-grep('00:', levels(bedtime))],levels(bedtime)[grep('00:', levels(bedtime))])
  
  bedtime <- as.numeric(bed_h) + as.numeric(bed_m)/60
  waketime <- as.numeric(wak_h) + as.numeric(wak_m)/60
  sleeptime <- (24 - bedtime + waketime) 
  
  cat('\nBED TIME:\n')
  print(summary(bedtime))
  cat('\nWAKE TIME:\n')
  print(summary(waketime))
  cat('\nSLEEP DURATION:\n')
  print(summary(sleeptime))
  
  return(data.frame(bedtime, waketime, sleeptime))
}

data[,c('sleep_bedtime','sleep_waketime','sleep_hr')] <- cleanSleep('kw4061', 'kw4060')

# Dichotomize based on recommendations
data$sleep_hr_bin <- as.factor(ifelse(data$sleep_hr < 9 | data$sleep_hr > 12, 0, 1))
summary(data$sleep_hr_bin)

data$sleep_hr_age <- full[,'kw9991a']/12

# Sleep auxiliary timepoints ---------------------------------------------------
# NOTE: the variable name is based on mean age at questionnaire completion
data[,c('sleep_bedt_04.8y','sleep_waket_04.8y','sleep_04.8y')] <- cleanSleep('kl222', 'kl223')
data[,c('sleep_bedt_05.8y','sleep_waket_05.8y','sleep_05.8y')] <- cleanSleep('kn2011','kn2020')
data[,c('sleep_bedt_06.8y','sleep_waket_06.8y','sleep_06.8y')] <- cleanSleep(c('kq252','kq253'), c('kq256', 'kq257'))
data[,c('sleep_bedt_09.7y','sleep_waket_09.7y','sleep_09.7y')] <- cleanSleep('ku341','ku342')

# ==============================================================================
# DIET (MEDITERRANEAN) (10.7 years) ============================================
# ==============================================================================

# Nutrients mean weight (g/day) [ FFQ ]
diet <- data.frame( # Extract form main dataset and rename
        # ------ CEREAL ---------------------------------
         'cer_1' = full[,'fddd200'], # High fibre breakfast cereals
         'cer_2' = full[,'fddd227'], # Brown and granary bread
         'cer_3' = full[,'fddd229'], # Wholemeal bread
         # (201) Other breakfast cereals /  (202) Sweet biscuits
         # (226) White bread / (228) Softgrain white bread / (230) Other bread
         # (266) Pasta, rice, pizza etc.
        
        # ------ FISH ---------------------------------
         'fis_1' = full[,'fddd204'], # Other white fish, shellfish, fish dishes
         'fis_2' = full[,'fddd205'], # Oily fish
         # (203) Coated and fried white fish, shellfish
        
        # ------ DAIRY -------------------------------
         'dai_1' = full[,'fddd206'], # Yoghurt and fromage frais
         'dai_2' = full[,'fddd252'], # Cheese (NOTE: typo in variable label)
         'dai_3' = full[,'fddd254'], # Whole milk 
         'dai_4' = full[,'fddd255'], # Semi-skimmed milk
         'dai_5' = full[,'fddd256'], # Skimmed milk
         'dai_6' = full[,'fddd257'], # Goats and sheeps milk
         # (207-212) ... sweets
         # (258) Soya milk / (259) Other milk and cream 
        
        # ------ MEAT -------------------------------- 
         'mea_1' = full[,'fddd215'], # Coated chicken and turkey
         'mea_2' = full[,'fddd216'], # Chicken, turkey and dishes
         'mea_3' = full[,'fddd217'], # Liver and dishes
         'mea_4' = full[,'fddd218'], # Lamb and dishes
         'mea_5' = full[,'fddd219'], # Pork and dishes
         'mea_6' = full[,'fddd220'], # Beef and dishes 
         'mea_7' = full[,'fddd221'], # Burgers and kebabs
         'mea_8' = full[,'fddd222'], # Sausages
         'mea_9' = full[,'fddd223'], # Offal (excluding liver)
        'mea_10' = full[,'fddd224'], # Other meat and meat products
        'mea_11' = full[,'fddd231'], # Butter
        'mea_12' = full[,'fddd238'], # Ham and bacon
         # (214) Meat pies and pastries 
         # (225) Eggs and egg dishes
         # (232) Full-fat polyunsaturated margarine / (233) Low-fat polyunsaturated margarine
         # (234) Full-fat non-polyunsaturated margarine / (235) Low-fat non-polyunsaturated margarine
         # (236) Polyunsaturated cooking fat / (237) Non-polyunsaturated cooking fat
        
        # ------ VEGETABLES -------------------------------- 
         'veg_1' = full[,'fddd241'], # Raw carrots
         'veg_2' = full[,'fddd242'], # Cooked carrots
         'veg_3' = full[,'fddd243'], # Green leafy vegetables
         'veg_4' = full[,'fddd248'], # Other salad and raw vegetables
         'veg_5' = full[,'fddd249'], # Other cooked vegetables
         'veg_6' = full[,'fddd251'], # Vegetable dishes
         # (239) Fried/roast potatoes and chips / (240) Other potatoes
        
        # ------ LEGUMES ----------------------------------
         'leg_1' = full[,'fddd213'], # Baked beans
         'leg_2' = full[,'fddd244'], # Peas
         'leg_3' = full[,'fddd245'], # Green and runner beans
         'leg_4' = full[,'fddd250'], # Legumes
        
        # ------ FRUIT --------------------------------
         'fru_1' = full[,'fddd246'], # Cooked and canned tomatoes
         'fru_2' = full[,'fddd247'], # Raw tomatoes
         'fru_3' = full[,'fddd261'], # Fruit canned in juice
         'fru_4' = full[,'fddd262'], # Citrus fruit
         'fru_5' = full[,'fddd263'], # Apples and pears
         'fru_6' = full[,'fddd264'], # Bananas
         'fru_7' = full[,'fddd265'], # Other fruit 
         'fru_8' = full[,'fddd267'], # Nuts
          # (253) Fruit juice / (260) Fruit canned in syrup
        
        # Optional: total energy intake -----------
        'tot_energy' = full[,'fddd304']
)

foodgroups <- c('veg','leg','fru','cer','fis','mea','dai')

# Compute a sum score per food group and the corresponding binary variable
for (foodgroup in foodgroups) {
  message('\n', toupper(foodgroup))
  
  name_c <- paste0('med_diet_',foodgroup) # Continuous 
  name_b <- paste0('med_diet_',foodgroup,'_bin') # Binary
  
  # Compute continuous (i.e., sum of food group components)
  diet[, name_c] <- rowSums(diet[, grepl(foodgroup, names(diet))], na.rm = F)
  print(summary(diet[, name_c]))
  
  # Dichotomize at median intake
  med <- median(diet[, name_c], na.rm = T)
  if (foodgroup %in% c('mea','dai')) { # there are reverse coded (i.e. non beneficial)
    diet[, name_b] <- ifelse(diet[, name_c] < med, 1, 0)
  } else if (foodgroup == 'fis') { # NOTE: median=0, so adopted a more lenient cutoff 
    diet[, name_b] <- ifelse(diet[, name_c] > med, 1, 0)
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

data$med_diet_age <- full[,'fd003c']/12

# Include single nutrients weight in the dataset, I also save a separate dataset
data <- cbind(data, diet[, names(diet)[grep('med_diet', names(diet))]])

# Supplementary: sum of continuous diet components ----------------------
# NOTE: WHAT FOR? isn't this just a total grams per day, not a diet score!
# data$med_diet_supp_c <- rowSums(diet[, paste0('med_diet_',foodgroups)]) 
# summary(data$med_diet_supp_c)

# Diet auxiliary timepoints ----------------------------------------------------

# Only one other timepoint available: 7.5 years mother report [ FFQ ]
# NOTE: these are largely, but not always the same numbers as 10.7 FFQ
diet_7.5y <- data.frame( # Extract form main dataset and rename
  # ------ CEREAL ---------------------------------
  'cer_1_7.5y' = full[,'f7dd200'], # High fibre breakfast cereals
  'cer_2_7.5y' = full[,'f7dd227'], # Brown and granary bread
  'cer_3_7.5y' = full[,'f7dd229'], # Wholemeal bread
  # ------ FISH ---------------------------------
  'fis_1_7.5y' = full[,'f7dd204'], # Other white fish, shellfish, fish dishes
  'fis_2_7.5y' = full[,'f7dd205'], # Oily fish
  # ------ DAIRY -------------------------------
  'dai_1_7.5y' = full[,'f7dd206'], # Yoghurt and fromage frais
  'dai_2_7.5y' = full[,'f7dd252'], # Cheese
  'dai_3_7.5y' = full[,'f7dd254'], # Whole milk
  'dai_4_7.5y' = full[,'f7dd255'], # Semi-skimmed milk
  'dai_5_7.5y' = full[,'f7dd256'], # Skimmed milk
  'dai_6_7.5y' = full[,'f7dd258'], # Goats and sheeps milk
  # ------ MEAT -------------------------------- 
  'mea_1_7.5y' = full[,'f7dd215'], # Coated chicken and turkey
  'mea_2_7.5y' = full[,'f7dd216'], # Chicken, turkey and dishes
  'mea_3_7.5y' = full[,'f7dd217'], # Liver and dishes
  'mea_4_7.5y' = full[,'f7dd218'], # Lamb and dishes
  'mea_5_7.5y' = full[,'f7dd219'], # Offal (excluding liver) 
  'mea_6_7.5y' = full[,'f7dd220'], # Pork and dishes 
  'mea_7_7.5y' = full[,'f7dd221'], # Beef and dishes 
  'mea_8_7.5y' = full[,'f7dd222'], # Burgers and kebabs
  'mea_9_7.5y' = full[,'f7dd223'], # Sausages
  'mea_10_7.5y' = full[,'f7dd224'], # Other meat and meat products
  'mea_11_7.5y' = full[,'f7dd231'], # Butter
  'mea_12_7.5y' = full[,'f7dd238'], # Ham and bacon
  # ------ VEGETABLES -------------------------------- 
  'veg_1_7.5y' = full[,'f7dd241'], # Raw carrots
  'veg_2_7.5y' = full[,'f7dd242'], # Cooked carrots
  'veg_3_7.5y' = full[,'f7dd243'], # Green leafy vegetables
  'veg_4_7.5y' = full[,'f7dd248'], # Other salad and raw vegetables
  'veg_5_7.5y' = full[,'f7dd249'], # Other cooked vegetables
  'veg_6_7.5y' = full[,'f7dd251'], # Vegetable dishes
  # ------ LEGUMES ----------------------------------
  'leg_1_7.5y' = full[,'f7dd213'], # Baked beans
  'leg_2_7.5y' = full[,'f7dd244'], # Peas
  'leg_3_7.5y' = full[,'f7dd245'], # Green and runner beans
  'leg_4_7.5y' = full[,'f7dd250'], # Legumes
  # ------ FRUIT --------------------------------
  'fru_1_7.5y' = full[,'f7dd246'], # Cooked and canned tomatoes
  'fru_2_7.5y' = full[,'f7dd247'], # Raw tomatoes
  'fru_3_7.5y' = full[,'f7dd261'], # Fruit canned in juice
  'fru_4_7.5y' = full[,'f7dd262'], # Citrus fruit
  'fru_5_7.5y' = full[,'f7dd263'], # Apples and pears
  'fru_6_7.5y' = full[,'f7dd264'], # Bananas
  'fru_7_7.5y' = full[,'f7dd265'], # Other fruit 
  'fru_8_7.5y' = full[,'f7dd267']  # Nuts
)

# Compute a sum score per food group and the corresponding binary variable
for (foodgroup in foodgroups) {
  message('\n', toupper(foodgroup))
  
  name_c <- paste0('med_diet_7.5y_',foodgroup) # Continuous 
  name_b <- paste0('med_diet_7.5y_',foodgroup,'_bin') # Binary
  
  # Compute continuous (i.e., sum of food group components)
  diet_7.5y[, name_c] <- rowSums(diet_7.5y[, grepl(foodgroup, names(diet_7.5y))], na.rm = F)
  print(summary(diet_7.5y[, name_c]))
  
  # Dichotomize at median intake
  med <- median(diet_7.5y[, name_c], na.rm = T)
  if (foodgroup %in% c('mea','dai')) { # there are reverse coded (i.e. non beneficial)
    diet_7.5y[, name_b] <- ifelse(diet_7.5y[, name_c] < med, 1, 0)
  } else if (foodgroup == 'fis') { # NOTE: median=0, so adopted a more lenient cutoff 
    diet_7.5y[, name_b] <- ifelse(diet_7.5y[, name_c] > med, 1, 0)
  } else {
    diet_7.5y[, name_b] <- ifelse(diet_7.5y[, name_c] >= med, 1, 0)
  }
  message("Dichotomized ", foodgroup, ' at median intake = ', med)
  print(summary(as.factor(diet_7.5y[, name_b])))
}

# Compute Mediterranean diet score at 7.5 years
data$med_diet_7.5y <- rowSums(diet_7.5y[, paste0('med_diet_7.5y_',foodgroups,'_bin')]) 
summary(data$med_diet_7.5y)

# Include single nutrients weight in the dataset, I also save a separate dataset
data <- cbind(data, diet_7.5y[, names(diet_7.5y)[grep('med_diet_7.5y_', names(diet_7.5y))]])

# ==============================================================================
# merge with the rest of the variables 
pren_stress <- readRDS(file.path(datapath, 'prenatal_stress.rds'))
post_stress <- readRDS(file.path(datapath, 'postnatal_stress.rds'))
outcome_cov <- readRDS(file.path(datapath, 'PCMout_cov_aux.rds'))

# NOTE the IDC in the outcome file was not stripped of extra spaces, let's fix this
outcome_cov$IDC <- gsub('\\s+', '', outcome_cov$IDC)

data_full <- Reduce(function(x,y) merge(x = x, y = y, by = 'IDC', all.x = T),
                   list(pren_stress, post_stress, outcome_cov, data) ) 

data_full$ELS <- rowSums(data_full[,c('prenatal_stress','postnatal_stress')])

# Scale main variables 
for (v in c('prenatal_stress','postnatal_stress','ELS',
            'intern_score_13','tot_fat_percent_13',
            'exercise','sleep_hr','med_diet')) {
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
data <- data.frame(lapply(data_full , as.numeric))
c = round(cor(data, use='pairwise.complete.obs'),2)

# Save and flee

# Main file
saveRDS(data, file.path(datapath, 'moderator_vars.rds'))
# Full data file
saveRDS(data_full, file.path(datapath, 'data_raw.rds'))
# Diet data for additional descriptive / check-ups
saveRDS(cbind(diet, diet_7.5y), file.path(datapath, 'diet_vars.rds'))
