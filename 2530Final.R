setwd('/Users/zhaoxiangding/Desktop/2530/Final')
library(dplyr)

# population data
################################
pop <- read.csv('MYEB1_detailed_population_estimates_series_UK_(2019).csv', header = T, fill = F)
pop <- pop[pop$country == 'E',]
pop <- pop[,c(4,5,18:24)]
pop <- pop %>% 
  group_by(age, sex) %>%
  summarise(pop_13 = sum(population_2013),
            pop_14 = sum(population_2014),
            pop_15 = sum(population_2015),
            pop_16 = sum(population_2016),
            pop_17 = sum(population_2017),
            pop_18 = sum(population_2018),
            pop_19 = sum(population_2019)) %>%
  mutate(agegroup = case_when(age == '0' ~ 'Under 1',
                              age %in% c(1:4) ~ '1-4',
                              age %in% c(5:9) ~ '5-9',
                              age %in% c(10:14) ~ '10-14',
                              age %in% c(15:19) ~ '15-19',
                              age %in% c(20:24) ~ '20-24',
                              age %in% c(25:29) ~ '25-29',
                              age %in% c(30:34) ~ '30-34',
                              age %in% c(35:39) ~ '35-39',
                              age %in% c(40:44) ~ '40-44',
                              age %in% c(45:49) ~ '45-49',
                              age %in% c(50:54) ~ '50-54',
                              age %in% c(55:59) ~ '55-59',
                              age %in% c(60:64) ~ '60-64',
                              age %in% c(65:69) ~ '65-69',
                              age %in% c(70:74) ~ '70-74',
                              age %in% c(75:79) ~ '75-79',
                              age %in% c(80:84) ~ '80-84',
                              age %in% c(85:89) ~ '85-89',
                              age == '90' ~ '90 and over')) %>%
  ungroup() %>%
  group_by(agegroup, sex) %>%
  summarise(pop_13 = sum(pop_13),
            pop_14 = sum(pop_14),
            pop_15 = sum(pop_15),
            pop_16 = sum(pop_16),
            pop_17 = sum(pop_17),
            pop_18 = sum(pop_18),
            pop_19 = sum(pop_19))

# cancer data
################################

cancer_id <- c('C15',
               'C16',
               'C18-C20',
               'C33-C34',
               'C50',
               'C67',
               'C82-C85',
               'C90',
               'C91-C95')

kidney_id <- c('C64', 'C65', 'C66', 'C68')

df_name <- c('C00', 'Site', 'Under 1', '1-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39',
             '40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80-84','85-89','90 and over',
             'sex', 'year') 

# cancer 17_male and female are duplicated
cancer17_male <- read.csv('cancer17_male.csv')
cancer17_kidney <- cancer17_male[cancer17_male$C00 %in% kidney_id,]
cancer17_kidney <- c('C64 to C66, C68', 'Kidney and Urinary Tract', apply(cancer17_kidney[,-c(1,2)],2,sum))
cancer17_male <- cancer17_male[cancer17_male$C00 %in% cancer_id,]
cancer17_male <- rbind(cancer17_male, cancer17_kidney)
cancer17_male$sex <- 'male'
cancer17_male$year <- 2017

cancer17_female <- read.csv('cancer17_demale.csv')
cancer17_kidney <- cancer17_female[cancer17_female$C00 %in% kidney_id,]
cancer17_kidney <- c('C64 to C66, C68', 'Kidney and Urinary Tract', apply(cancer17_kidney[,-c(1,2)],2,sum))
cancer17_female <- cancer17_female[cancer17_female$C00 %in% cancer_id,]
cancer17_female <- rbind(cancer17_female, cancer17_kidney)
cancer17_female$sex <- 'female'
cancer17_female$year <- 2017
colnames(cancer17_male) <- df_name
colnames(cancer17_female) <- df_name
cancer17 <- rbind(cancer17_male, cancer17_female)
colnames(cancer17) <- df_name
cancer16 <- read.csv('cancer16.csv')
cancer16 <- cancer16[rowSums(is.na(cancer16)) == 0 ,]
cancer16_kidney <- cancer16[cancer16$C00 %in% kidney_id,]
cancer16_kidney_male <- cancer16_kidney[cancer16_kidney$Males == 'Males', ]
cancer16_kidney_female <- cancer16_kidney[cancer16_kidney$Males == 'Females', ]
cancer16_kidney_male <- c('C64 to C66, C68', 'Kidney and Urinary Tract', 'Males', apply(cancer16_kidney_male[,-c(1:3)],2,sum))
cancer16_kidney_female <- c('C64 to C66, C68', 'Kidney and Urinary Tract','Females', apply(cancer16_kidney_female[,-c(1:3)],2,sum))
cancer16 <- cancer16[cancer16$C00 %in% cancer_id,]
cancer16 <- rbind(cancer16, cancer16_kidney_male, cancer16_kidney_female)
cancer16 <- cancer16 %>%
  mutate(sex = ifelse(Males == 'Males', 'male', 'female'),
         year = 2016) %>% 
  select(-Males)
colnames(cancer16) <- df_name
cancer <- rbind(cancer16, cancer17)

cancer15 <- read.csv('cancer15.csv')
cancer15 <- cancer15[rowSums(is.na(cancer15)) == 0 ,]
cancer15 <- cancer15[, -4]
cancer15_kidney <- cancer15[cancer15$C00 %in% kidney_id,]
cancer15_kidney_male <- cancer15_kidney[cancer15_kidney$Males == 'Males', ]
cancer15_kidney_female <- cancer15_kidney[cancer15_kidney$Males == 'Females', ]
cancer15_kidney_male <- c('C64 to C66, C68', 'Kidney and Urinary Tract', 'Males', apply(cancer15_kidney_male[,-c(1:3)],2,sum))
cancer15_kidney_female <- c('C64 to C66, C68', 'Kidney and Urinary Tract','Females', apply(cancer15_kidney_female[,-c(1:3)],2,sum))
cancer15 <- cancer15[cancer15$C00 %in% cancer_id,]
cancer15 <- rbind(cancer15, cancer15_kidney_female, cancer15_kidney_male)
cancer15 <- cancer15 %>%
  mutate(sex = ifelse(Males == 'Males', 'male', 'female'),
         year = 2015) %>%
  select(-Males)
colnames(cancer15) <- df_name
cancer <- rbind(cancer, cancer15)

cancer14 <- read.csv('cancer14.csv')
cancer14 <- cancer14[rowSums(is.na(cancer14)) == 0 ,]
cancer14 <- cancer14[, -4]
cancer14_kidney <- cancer14[cancer14$C00 %in% kidney_id,]
cancer14_kidney_male <- cancer14_kidney[cancer14_kidney$Males == 'Males', ]
cancer14_kidney_female <- cancer14_kidney[cancer14_kidney$Males == 'Females', ]
cancer14_kidney_male <- c('C64 to C66, C68', 'Kidney and Urinary Tract', 'Males', apply(cancer14_kidney_male[,-c(1:3)],2,sum))
cancer14_kidney_female <- c('C64 to C66, C68', 'Kidney and Urinary Tract','Females', apply(cancer14_kidney_female[,-c(1:3)],2,sum))
cancer14 <- cancer14[cancer14$C00 %in% cancer_id,]
cancer14 <- rbind(cancer14, cancer14_kidney_female, cancer14_kidney_male)
cancer14 <- cancer14 %>%
  mutate(sex =  ifelse(Males == 'Males', 'male', 'female'),
         year = 2014) %>%
  select(-Males)
colnames(cancer14) <- df_name
cancer <- rbind(cancer, cancer14)

cancer13 <- read.csv('cancer13.csv')
cancer13 <- cancer13[rowSums(is.na(cancer13)) == 20, ]
cancer13 <- cancer13[,colSums(is.na(cancer13)) == 0]
cancer13 <- cancer13[-1, ]
for (i in 1:nrow(cancer13)){
  if (cancer13[i,1] == ''){
    cancer13[i,1] <- cancer13[i-1,1]
    cancer13[i,2] <- cancer13[i-1,2]
  }
}
cancer13_kidney <- cancer13[cancer13$C00 %in% kidney_id,]
cancer13_kidney_male <- cancer13_kidney[cancer13_kidney$M == 'M', ]
cancer13_kidney_female <- cancer13_kidney[cancer13_kidney$M == 'F', ]
cancer13_kidney_male <- c('C64 to C66, C68', 'Kidney and Urinary Tract', 'M', apply(cancer13_kidney_male[,-c(1:3)],2,sum))
cancer13_kidney_female <- c('C64 to C66, C68', 'Kidney and Urinary Tract','F', apply(cancer13_kidney_female[,-c(1:3)],2,sum))
cancer13 <- cancer13[cancer13$C00 %in% cancer_id,]
cancer13 <- rbind(cancer13, cancer13_kidney_female, cancer13_kidney_male)
cancer13 <- cancer13 %>%
  mutate(sex = ifelse(M == 'M', 'male', 'female'),
         year = 2013) %>%
  select(-M)
colnames(cancer13) <- df_name
cancer <- rbind(cancer, cancer13)


# Preprocess
################################

cancer13 <- cancer %>% 
  filter(year == 2013) %>%
  group_by(C00) %>%
  summarise(`Under 1` = sum(as.numeric(`Under 1`)),
            `1-4` = sum(as.numeric(`1-4`)),
            `5-9` = sum(as.numeric(`5-9`)),
            `10-14` = sum(as.numeric(`10-14`)),
            `15-19` = sum(as.numeric(`15-19`)),
            `20-24` = sum(as.numeric(`20-24`)),
            `25-29` = sum(as.numeric(`25-29`)),
            `30-34` = sum(as.numeric(`30-34`)),
            `35-39` = sum(as.numeric(`35-39`)),
            `40-44` = sum(as.numeric(`40-44`)),
            `45-49` = sum(as.numeric(`45-49`)),
            `50-54` = sum(as.numeric(`50-54`)),
            `55-59` = sum(as.numeric(`55-59`)),
            `60-64` = sum(as.numeric(`60-64`)),
            `65-69` = sum(as.numeric(`65-69`)),
            `70-74` = sum(as.numeric(`70-74`)),
            `75-79` = sum(as.numeric(`75-79`)),
            `80-84` = sum(as.numeric(`80-84`)),
            `85-89` = sum(as.numeric(`85-89`)),
            `90 and over` = sum(as.numeric(`90 and over`))) %>%
  arrange(C00)

cancer14 <- cancer %>%
  filter(year == 2014) %>%
  group_by(C00) %>%
  summarise(`Under 1` = sum(as.numeric(`Under 1`)),
            `1-4` = sum(as.numeric(`1-4`)),
            `5-9` = sum(as.numeric(`5-9`)),
            `10-14` = sum(as.numeric(`10-14`)),
            `15-19` = sum(as.numeric(`15-19`)),
            `20-24` = sum(as.numeric(`20-24`)),
            `25-29` = sum(as.numeric(`25-29`)),
            `30-34` = sum(as.numeric(`30-34`)),
            `35-39` = sum(as.numeric(`35-39`)),
            `40-44` = sum(as.numeric(`40-44`)),
            `45-49` = sum(as.numeric(`45-49`)),
            `50-54` = sum(as.numeric(`50-54`)),
            `55-59` = sum(as.numeric(`55-59`)),
            `60-64` = sum(as.numeric(`60-64`)),
            `65-69` = sum(as.numeric(`65-69`)),
            `70-74` = sum(as.numeric(`70-74`)),
            `75-79` = sum(as.numeric(`75-79`)),
            `80-84` = sum(as.numeric(`80-84`)),
            `85-89` = sum(as.numeric(`85-89`)),
            `90 and over` = sum(as.numeric(`90 and over`))) %>%
  arrange(C00)

cancer15 <- cancer %>%
  filter(year == 2015) %>%
  group_by(C00) %>%
  summarise(`Under 1` = sum(as.numeric(`Under 1`)),
            `1-4` = sum(as.numeric(`1-4`)),
            `5-9` = sum(as.numeric(`5-9`)),
            `10-14` = sum(as.numeric(`10-14`)),
            `15-19` = sum(as.numeric(`15-19`)),
            `20-24` = sum(as.numeric(`20-24`)),
            `25-29` = sum(as.numeric(`25-29`)),
            `30-34` = sum(as.numeric(`30-34`)),
            `35-39` = sum(as.numeric(`35-39`)),
            `40-44` = sum(as.numeric(`40-44`)),
            `45-49` = sum(as.numeric(`45-49`)),
            `50-54` = sum(as.numeric(`50-54`)),
            `55-59` = sum(as.numeric(`55-59`)),
            `60-64` = sum(as.numeric(`60-64`)),
            `65-69` = sum(as.numeric(`65-69`)),
            `70-74` = sum(as.numeric(`70-74`)),
            `75-79` = sum(as.numeric(`75-79`)),
            `80-84` = sum(as.numeric(`80-84`)),
            `85-89` = sum(as.numeric(`85-89`)),
            `90 and over` = sum(as.numeric(`90 and over`))) %>%
  arrange(C00)

cancer16 <- cancer %>%
  filter(year == 2016) %>%
  group_by(C00) %>%
  summarise(`Under 1` = sum(as.numeric(`Under 1`)),
            `1-4` = sum(as.numeric(`1-4`)),
            `5-9` = sum(as.numeric(`5-9`)),
            `10-14` = sum(as.numeric(`10-14`)),
            `15-19` = sum(as.numeric(`15-19`)),
            `20-24` = sum(as.numeric(`20-24`)),
            `25-29` = sum(as.numeric(`25-29`)),
            `30-34` = sum(as.numeric(`30-34`)),
            `35-39` = sum(as.numeric(`35-39`)),
            `40-44` = sum(as.numeric(`40-44`)),
            `45-49` = sum(as.numeric(`45-49`)),
            `50-54` = sum(as.numeric(`50-54`)),
            `55-59` = sum(as.numeric(`55-59`)),
            `60-64` = sum(as.numeric(`60-64`)),
            `65-69` = sum(as.numeric(`65-69`)),
            `70-74` = sum(as.numeric(`70-74`)),
            `75-79` = sum(as.numeric(`75-79`)),
            `80-84` = sum(as.numeric(`80-84`)),
            `85-89` = sum(as.numeric(`85-89`)),
            `90 and over` = sum(as.numeric(`90 and over`))) %>%
  arrange(C00)

cancer17 <- cancer %>%
  filter(year == 2017) %>%
  group_by(C00) %>%
  summarise(`Under 1` = sum(as.numeric(`Under 1`)),
            `1-4` = sum(as.numeric(`1-4`)),
            `5-9` = sum(as.numeric(`5-9`)),
            `10-14` = sum(as.numeric(`10-14`)),
            `15-19` = sum(as.numeric(`15-19`)),
            `20-24` = sum(as.numeric(`20-24`)),
            `25-29` = sum(as.numeric(`25-29`)),
            `30-34` = sum(as.numeric(`30-34`)),
            `35-39` = sum(as.numeric(`35-39`)),
            `40-44` = sum(as.numeric(`40-44`)),
            `45-49` = sum(as.numeric(`45-49`)),
            `50-54` = sum(as.numeric(`50-54`)),
            `55-59` = sum(as.numeric(`55-59`)),
            `60-64` = sum(as.numeric(`60-64`)),
            `65-69` = sum(as.numeric(`65-69`)),
            `70-74` = sum(as.numeric(`70-74`)),
            `75-79` = sum(as.numeric(`75-79`)),
            `80-84` = sum(as.numeric(`80-84`)),
            `85-89` = sum(as.numeric(`85-89`)),
            `90 and over` = sum(as.numeric(`90 and over`))) %>%
  arrange(C00)

pop1 <- pop %>%
  dplyr::select(-sex) %>%
  group_by(agegroup) %>%
  summarise(pop_13 = sum(pop_13),
         pop_14 = sum(pop_14),
         pop_15 = sum(pop_15),
         pop_16 = sum(pop_16),
         pop_17 = sum(pop_17),
         pop_18 = sum(pop_18),
         pop_19 = sum(pop_19)) %>%
  mutate(agegroup = factor(agegroup, levels = c('Under 1', '1-4', '5-9', '10-14', 
                                                '15-19', '20-24', '25-29', '30-34',
                                                '35-39', '40-44', '45-49', '50-54', 
                                                '55-59', '60-64', '65-69', '70-74', 
                                                '75-79', '80-84', '85-89', '90 and over'))) %>%
  arrange(agegroup)

cancer_cor_df <- cancer %>%
  select(-c(2,23,24)) %>%
  group_by(C00) %>%
  summarise(`Under 1` = sum(as.numeric(`Under 1`)),
            `1-4` = sum(as.numeric(`1-4`)),
            `5-9` = sum(as.numeric(`5-9`)),
            `10-14` = sum(as.numeric(`10-14`)),
            `15-19` = sum(as.numeric(`15-19`)),
            `20-24` = sum(as.numeric(`20-24`)),
            `25-29` = sum(as.numeric(`25-29`)),
            `30-34` = sum(as.numeric(`30-34`)),
            `35-39` = sum(as.numeric(`35-39`)),
            `40-44` = sum(as.numeric(`40-44`)),
            `45-49` = sum(as.numeric(`45-49`)),
            `50-54` = sum(as.numeric(`50-54`)),
            `55-59` = sum(as.numeric(`55-59`)),
            `60-64` = sum(as.numeric(`60-64`)),
            `65-69` = sum(as.numeric(`65-69`)),
            `70-74` = sum(as.numeric(`70-74`)),
            `75-79` = sum(as.numeric(`75-79`)),
            `80-84` = sum(as.numeric(`80-84`)),
            `85-89` = sum(as.numeric(`85-89`)),
            `90 and over` = sum(as.numeric(`90 and over`))) %>%
  arrange(C00) %>%
  t() %>%
  as.data.frame()
colnames(cancer_cor_df) <- cancer_cor_df[1,]
cancer_cor_df <- cancer_cor_df[-1,]
cancer_cor_df <- apply(cancer_cor_df, 2, as.numeric)
cancer_cor <- cor(as.matrix(cancer_cor_df))

indicator <- diag(1, nrow = ncol(cancer13)-1)
indicator[10,] <- 1

pop1 <- pop1[,-c(1,7,8)]
################################

#model
################################
library(rstan)
mod = stan_model(file = 'new-model1.stan')
sample = sampling(mod, data = list(N = 5, 
                                   K = ncol(cancer13) - 1, 
                                   M = ncol(cancer_cor), 
                                   y13 = as.matrix(cancer13[,-1]),
                                   y14 = as.matrix(cancer14[,-1]),
                                   y15 = as.matrix(cancer15[,-1]),
                                   y16 = as.matrix(cancer16[,-1]),
                                   y17 = as.matrix(cancer17[,-1]),
                                   pop = as.matrix(pop1),
                                   X = indicator,
                                   tar_sigma1 = 100,
                                   tar_sigma2 = 100,
                                   rho1 = 0.4,
                                   rho2 = 0.4),
                  chains = 2, iter = 1000)



#saveRDS(sample, 'sample_best.rds')
#sample = readRDS('sample_best.rds')
model_p = extract(sample, pars='p_post')$p_post
model_theta = extract(sample, pars='theta_post')$theta_post

plot(sample, plotfun = 'trace', pars = 'alpha2', inc_warmup = F)
plot(sample, plotfun = 'trace', pars = 'beta', inc_warmup = F)


plot(model_p[1:1000,1,5,1], type = "l", col = "blue")
lines(model_p[1001:2000,1,5,1], col = "red")
lines(model_p[2001:3000,1,5,1], col = "green")

plot(model_theta[1:250,5,5,10], type = "l", col = "blue")
lines(model_theta[251:500,5,5,10], col = "red")
lines(model_theta[501:750,5,5,10], col = "green")


plot(model_p[500,1,5,1:20], type = "l", col = "blue")
plot(model_theta[500,1,5,1:20], type = "l", col = "blue")


#saveRDS(sample_y, 'sample_y_best.rds')
sample_y = readRDS('sample_y_best.rds')

plot(sample_y, plotfun = 'trace', pars = 'alpha1', inc_warmup = F)
plot(sample_y, plotfun = 'trace', pars = 'beta', inc_warmup = F)
plot(sample_y, plotfun = 'trace', pars = 'alpha2', inc_warmup = F)


plot(model_p[1:500,1,3,1], type = "l", col = "blue")
lines(model_p[501:1000,1,3,1], col = "red")



plot(model_theta[1:500,1,3,20], type = "l", col = "blue", ylab = c("theta"),
     main = "Trace plot of theta in year 2013, age group: 90+, Cancer type: Oesophagus")
lines(model_theta[501:1000,1,3,20], col = "red")
lines(model_theta[1001:1500,1,3,20], col = "green")

est_y <- matrix(0, nrow = 10, ncol = 20)
r_17 <- matrix(0, nrow = 10, ncol = 20)
for (i in 1:10){
  for (j in 1:20){
    est_y[i,j] = mean(model_p[1:1000,i,5,j])/(1-exp(-mean(model_theta[1:1000,i,5,j])))*mean(model_theta[1:1000,i,5,j])
    r_17[i,j] <- as.numeric((est_y[i,j]/pop1[j,5])*100000)
  }
}
mse_17 <- cancer17[,-1] - est_y


est_y <- matrix(0, nrow = 10, ncol = 20)
r_16 <- matrix(0, nrow = 10, ncol = 20)
for (i in 1:10){
  for (j in 1:20){
    est_y[i,j] = mean(model_p[1:1000,i,4,j])/(1-exp(-mean(model_theta[1:1000,i,4,j])))*mean(model_theta[1:1000,i,4,j])
    r_16[i,j] = as.numeric((est_y[i,j]/pop1[j,4])*100000)
  }
}
mse_16 <- cancer16[,-1] - est_y

est_y <- matrix(0, nrow = 10, ncol = 20)
r_15 <- matrix(0, nrow = 10, ncol = 20)
for (i in 1:10){
  for (j in 1:20){
    est_y[i,j] = mean(model_p[1:1000,i,3,j])/(1-exp(-mean(model_theta[1:1000,i,3,j])))*mean(model_theta[1:1000,i,3,j])
    r_15[i,j] = as.numeric((est_y[i,j]/pop1[j,3])*100000)
  }
}
mse_15 <- cancer15[,-1] - est_y

est_y <- matrix(0, nrow = 10, ncol = 20)
r_14 <- matrix(0, nrow = 10, ncol = 20)
for (i in 1:10){
  for (j in 1:20){
    est_y[i,j] = mean(model_p[1:1000,i,2,j])/(1-exp(-mean(model_theta[1:1000,i,2,j])))*mean(model_theta[1:1000,i,2,j])
    r_14[i,j] = as.numeric((est_y[i,j]/pop1[j,2])*100000)
  }
}
mse_14 <- cancer14[,-1] - est_y

est_y <- matrix(0, nrow = 10, ncol = 20)
r_13 <- matrix(0, nrow = 10, ncol = 20)
for (i in 1:10){
  for (j in 1:20){
    est_y[i,j] = mean(model_p[1:1000,i,1,j])/(1-exp(-mean(model_theta[1:1000,i,1,j])))*mean(model_theta[1:1000,i,1,j])
    r_13[i,j] = as.numeric((est_y[i,j]/pop1[j,1])*100000)
  }
}
mse_13 <- cancer13[,-1] - est_y



hist(c(as.vector(t(mse_17)), as.vector(t(mse_16)),as.vector(t(mse_15)), as.vector(t(mse_14)), as.vector(t(mse_13))), breaks = 20, main = "Histogram of model's prediction compare to observations", xlab = "MSE")

asr_17 <- vector()
for (i in 1:10){
  asr_17[i] <- sum(r_17[i,] *standard_pop[,5])/sum(standard_pop[,5])
}

asr_16 <- vector()
for (i in 1:10){
  asr_16[i] <- sum(r_16[i,] *standard_pop[,4])/sum(standard_pop[,4])
}

asr_15 <- vector()
for (i in 1:10){
  asr_15[i] <- sum(r_15[i,] *standard_pop[,3])/sum(standard_pop[,3])
}

asr_14 <- vector()
for (i in 1:10){
  asr_14[i] <- sum(r_14[i,] *standard_pop[,2])/sum(standard_pop[,2])
}

asr_13 <- vector()
for (i in 1:10){
  asr_13[i] <- sum(r_13[i,] *standard_pop[,1])/sum(standard_pop[,1])
}

asr_matrix <- matrix(c(asr_13, asr_14, asr_15, asr_16, asr_17), nrow = 10, ncol = 5)
asr_matrix <- round(asr_matrix, 2)
write.csv(asr_matrix, 'asr_matrix.csv')
