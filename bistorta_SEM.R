#20200225
#piecewiseSEM on linear mixed-effects models
#script written by M.Wutkowska and D. Ehrich


R.version #R version 3.5.2 (2018-12-20)
#clean the environment
rm(list = ls(all = TRUE))
options(max.print = 1000000)


#loading packages
library(nlme) #v. 3.1-137
library(MASS) #v. 7.3-51.1
library(piecewiseSEM) #v. 2.0.2
library(tidyverse) #v. 1.2.1


#loading environmental data and Bistorta measurements
env_bistorta <- read.delim('env_bistorta_temp.txt', header = TRUE, row.names = 1, sep = '\t')
dim(env_bistorta) #[1] 214  56
summary(env_bistorta)

cor(env_bistorta$ASVs_rarefied, env_bistorta$mean_July_temp_yearofsampling) #[1] -0.1100928

###DISTRIBUTIONS OF VARIABLES
# plants
hist(env_bistorta$rhizome_volume)
hist(log10(env_bistorta$rhizome_volume), breaks=15)
hist(env_bistorta$leaf_length)
hist(log10(env_bistorta$leaf_length), breaks=10)
hist(env_bistorta$ratioInflToStem) 
table(env_bistorta$ratioInflToStem) 

#fungi
hist(env_bistorta$ASVs_rarefied)
hist(env_bistorta$ratio_symbio_sapro, breaks=15) 
hist(log10(env_bistorta$ratio_symbio_sapro), breaks=15) 
hist(env_bistorta$ratio_symbio_sapro_reads) 
hist(log10(env_bistorta$ratio_symbio_sapro_reads), breaks=15)
hist(env_bistorta$composition_abu, breaks=15) 
hist(env_bistorta$composition_pa, breaks=15) 
hist(env_bistorta$Shannon_H) 
hist(env_bistorta$Pielou_J)


# soil
hist(env_bistorta$soil_nitrogen) #mg / g
hist(env_bistorta$soil_carbon) #mg / g
hist((log(env_bistorta$soil_nitrogen)/(1-log(env_bistorta$soil_nitrogen))), breaks=15) 
hist((log(env_bistorta$soil_nitrogen/1000)/(1-log(env_bistorta$soil_nitrogen/1000))), breaks=15) # taking it as permil
hist((log(env_bistorta$soil_carbon)/(1-log(env_bistorta$soil_carbon))), breaks=15) 
hist((log(env_bistorta$soil_carbon/1000)/(1-log(env_bistorta$soil_carbon/1000))), breaks=15) # taking it as permil
hist(env_bistorta$CN_ratio, breaks=20) 
hist(log10(env_bistorta$CN_ratio), breaks=20) 
hist(env_bistorta$pH, breaks=20) 
hist(env_bistorta$organic_matter_percentage, breaks=20) 

# wheather
hist(env_bistorta$prec_sumyear) 
hist(env_bistorta$tempC_meanyear) 
hist(env_bistorta$mean_July_temp_yearpriortosampling) 
hist(env_bistorta$mean_July_temp_yearofsampling) 

# Julian data
env_bistorta$jul <- strptime(env_bistorta$sampling_date, "%d.%m.%Y")$yday


# scaling and renaming variables
#plant
env_bistorta$rv <- scale(log10(env_bistorta$rhizome_volume))
env_bistorta$ll <- scale(log10(env_bistorta$leaf_length))
env_bistorta$i_s <- scale(env_bistorta$ratioInflToStem)

#fungi 
env_bistorta$asv <- scale(env_bistorta$ASVs_rarefied)
env_bistorta$ratt <- scale(log10(env_bistorta$ratio_symbio_sapro))
env_bistorta$ratr <- scale(log10(env_bistorta$ratio_symbio_sapro_reads))
env_bistorta$ccabu <- scale(env_bistorta$composition_abu)
env_bistorta$ccpa <- scale(env_bistorta$composition_pa)
env_bistorta$h <- scale(env_bistorta$Shannon_H)
env_bistorta$j <- scale(env_bistorta$Pielou_J)

#soil
env_bistorta$N <- scale(log(env_bistorta$soil_nitrogen/1000)/(1-log(env_bistorta$soil_nitrogen/1000)))
env_bistorta$C <- scale(log(env_bistorta$soil_carbon/1000)/(1-log(env_bistorta$soil_carbon/1000)))
env_bistorta$CN <- scale(log10(env_bistorta$CN_ratio))
env_bistorta$ph <- scale(env_bistorta$pH)
env_bistorta$om <- scale(env_bistorta$organic_matter_percentage)

#climate
env_bistorta$p <- scale(env_bistorta$prec_sumyear)
env_bistorta$t <- scale(env_bistorta$mean_July_temp_yearofsampling)

cor(env_bistorta$h, env_bistorta$j) #[1,] 0.9496197



################
##rv, is, ll ###
################
#because some values had a NaN/NA status in plant variables:  
env_bistorta_ll_rv_is <- env_bistorta[, c("leaf_length", "rhizome_volume", "ratioInflToStem", "N", "CN", "pH", "p", "t","Locality", "asv", "ratr", "ratt", "ccabu", "ccpa", "h", "j")]
env_bistorta_ll_rv_is <- env_bistorta_ll_rv_is %>%
  drop_na(leaf_length, rhizome_volume, ratioInflToStem)

dim(env_bistorta_ll_rv_is) #188 16 #full cases with all the data present

# SCALING and RENAMING PLANT VARIABLES
env_bistorta_ll_rv_is$rv <- scale(log10(env_bistorta_ll_rv_is$rhizome_volume))
env_bistorta_ll_rv_is$ll <- scale(log10(env_bistorta_ll_rv_is$leaf_length))
env_bistorta_ll_rv_is$i_s <- scale(env_bistorta_ll_rv_is$ratioInflToStem)

# check for correlations between the abiotic variables
cor(env_bistorta_ll_rv_is[ , c("N", "CN", "pH", "p", "t")])



############################################################
### Assesing lme models that are going to be used in SEM ###
############################################################

###Presence-absence
mod_asv <- lme(asv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ll_rv_is)
summary(mod_asv)
plot(mod_asv)
qqnorm(mod_asv) #OK

mod_ratt <- lme(ratt ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ll_rv_is)
summary(mod_ratt)
plot(mod_ratt)
qqnorm(mod_ratt) #OK

mod_ccpa <- lme(ccpa ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ll_rv_is)
summary(mod_ccpa)
plot(mod_ccpa)
qqnorm(mod_ccpa) #OK?

mod_PA_is <- lme(i_s ~ rv + N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_ll_rv_is)
summary(mod_PA_is)
plot(mod_PA_is)
qqnorm(mod_PA_is)
env_bistorta_ll_rv_is[residuals(mod_PA_is) > 6, ] # C7e from snow fences
env_bistorta_PA <- env_bistorta_ll_rv_is[-c(138),]

mod_PA_rv <- lme(rv ~ N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_ll_rv_is)
summary(mod_PA_rv)
plot(mod_PA_rv)
qqnorm(mod_PA_rv) #OK

mod_PA_ll <- lme(ll ~ rv + N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_ll_rv_is)
summary(mod_PA_ll)
plot(mod_PA_ll)
qqnorm(mod_PA_ll) #OK




###Abundance
mod_h <- lme(h ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ll_rv_is)
summary(mod_h)
plot(mod_h)
qqnorm(mod_h)

env_bistorta_ll_rv_is[residuals(mod_h) < -3, ] #B6e
rownames(env_bistorta_ll_rv_is)
env_bistorta_ABU <- env_bistorta_ll_rv_is[-c(134),]


mod_ratr <- lme(ratr ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ll_rv_is)
summary(mod_ratr)
plot(mod_ratr)        
env_bistorta_ll_rv_is[residuals(mod_ratr) > 3, ] #2.1.7
rownames(env_bistorta_ABU)
env_bistorta_ABU <- env_bistorta_ABU[-c(75),]


mod_ccabu <- lme(ccabu ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ll_rv_is)
summary(mod_ccabu)
plot(mod_ccabu) #OK


mod_ABU_is <- lme(i_s ~ N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ll_rv_is)  
summary(mod_ABU_is)
plot(mod_ABU_is)
env_bistorta_ll_rv_is[residuals(mod_ABU_is) > 6, ] #C7e
rownames(env_bistorta_ABU)
env_bistorta_ABU <- env_bistorta_ABU[-c(136),]


mod_ABU_rv <- lme(rv ~ N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ll_rv_is)
summary(mod_ABU_rv)
plot(mod_ABU_rv)
qqnorm(mod_ABU_rv) #OK


mod_ABU_ll <- lme(ll ~ rv + N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ll_rv_is)
summary(mod_ABU_ll)
plot(mod_ABU_ll)
qqnorm(mod_ABU_ll) #OK



#################################################################
#####SEM with mean July temperature from the year of sampling ###
#################################################################

### 1 FULL MODEL - the plant measures depend on both soil and fungi
bist.sem.lme_PA1 <- psem(
  lme(asv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ratt ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ccpa ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(i_s ~ rv + N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_PA),
  lme(rv ~ N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_PA),
  lme(ll ~ rv + N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_PA),
  ccpa %~~% ratt,
  data = env_bistorta_PA)
summary(bist.sem.lme_PA1)  #AIC: 121.229, Fisher's C = 3.229 with P-value = 0.78 and on 6 degrees of freedom

bist.sem.lme_ABU1 <- psem(
  lme(h ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ratr ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ccabu ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(i_s ~ N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ABU),
  lme(rv ~ N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ABU),
  lme(ll ~ rv + N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ABU),
  ratr %~~% h,
  data = env_bistorta_ABU)
summary(bist.sem.lme_ABU1) #AIC: 123.067, Fisher's C = 7.067 with P-value = 0.529 and on 8 degrees of freedom



###2 - the plant measures don't depend on fungi at all.
bist.sem.lme_PA2 <- psem(
  lme(asv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ratt ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ccpa ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(i_s ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(rv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA), 
  lme(ll ~ rv + N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  ratt %~~% ccpa,
  data = env_bistorta_PA)
summary(bist.sem.lme_PA2)  
#AIC: 140.860, Fisher's C = 42.86 with P-value = 0.02 and on 26 degrees of freedom

# And the same based on abundance
bist.sem.lme_ABU2 <- psem(
  lme(h ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ratr ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ccabu ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(i_s ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(rv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ll ~ rv + N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  ratr %~~% h,
  data = env_bistorta_ABU)
summary(bist.sem.lme_ABU2) 
#AIC = 119.276, Fisher's C = 21.276 with P-value = 0.728 and on 26 degrees of freedom



###3 i_s does not depend on fungi
bist.sem.lme_PA3 <- psem(
  lme(asv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ratt ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ccpa ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(i_s ~  N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(rv ~ N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_PA),
  lme(ll ~ rv + N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_PA),
  ccpa %~~% ratt,
  data = env_bistorta_PA)
summary(bist.sem.lme_PA3) 
#AIC: 118.900, Fisher's C = 8.9 with P-value = 0.837 and on 14 degrees of freedom

bist.sem.lme_ABU3 <- psem(
  lme(h ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ratr ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ccabu ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(i_s ~  N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(rv ~ N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ABU),
  lme(ll ~ rv + N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ABU),
  ratr %~~% h,
  data = env_bistorta_ABU)
summary(bist.sem.lme_ABU3) 
# AIC: 121.680, Fisher's C = 11.68 with P-value = 0.632 and on 14 degrees of freedom



###4 ll does not depend on fungi
bist.sem.lme_PA4 <- psem(
  lme(asv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ratt ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ccpa ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(i_s ~  N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_PA),
  lme(rv ~ N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_PA),
  lme(ll ~ rv + N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  ccpa %~~% ratt,
  data = env_bistorta_PA)
summary(bist.sem.lme_PA4) 
#AIC: 132.307, Fisher's C = 22.307 with P-value = 0.073 and on 14 degrees of freedom

bist.sem.lme_ABU4 <- psem(
  lme(h ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ratr ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ccabu ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(i_s ~  N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ABU),
  lme(rv ~ ll + N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ABU),
  lme(ll ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  ratr %~~% h,
  data = env_bistorta_ABU)
summary(bist.sem.lme_ABU4) 
#AIC: 123.442, Fisher's C = 13.442 with P-value = 0.492 and on 14 degrees of freedom



###5 rv does not depend on fungi
bist.sem.lme_PA5 <- psem(
  lme(asv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ratt ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ccpa ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(i_s ~  N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_PA),
  lme(rv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ll ~ rv + N + CN + pH + p + t + asv + ratt + ccpa, random=~1|Locality, data = env_bistorta_PA),
  ccpa %~~% ratt,
  data = env_bistorta_PA)
summary(bist.sem.lme_PA5) 
#AIC: 138.000, isher's C = 28 with P-value = 0.014 and on 14 degrees of freedom

bist.sem.lme_ABU5 <- psem(
  lme(h ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ratr ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ccabu ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(i_s ~  N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ABU),
  lme(rv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ll ~ rv + N + CN + pH + p + t + h + ratr + ccabu, random=~1|Locality, data = env_bistorta_ABU),
  ratr %~~% h,
  data = env_bistorta_ABU)
summary(bist.sem.lme_ABU5) 
#AIC: 123.577, Fisher's C = 13.577 with P-value = 0.482 and on 14 degrees of freedom



###6 community composition does not impact plants
bist.sem.lme_PA6 <- psem(
  lme(asv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ratt ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ccpa ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(i_s ~ N + CN + pH + p + t + asv + ratt, random=~1|Locality, data = env_bistorta_PA),
  lme(rv ~ N + CN + pH + p + t + asv + ratt, random=~1|Locality, data = env_bistorta_PA),
  lme(ll ~ rv + N + CN + pH + p + t + asv + ratt, random=~1|Locality, data = env_bistorta_PA),
  ccpa %~~% ratt,
  data = env_bistorta_PA)
summary(bist.sem.lme_PA6)  
#AIC: 120.317, Fisher's C = 10.317 with P-value = 0.739 and on 14 degrees of freedom

bist.sem.lme_ABU6 <- psem(
  lme(h ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ratr ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ccabu ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(i_s ~ N + CN + pH + p + t + h + ratr, random=~1|Locality, data = env_bistorta_ABU),
  lme(rv ~ N + CN + pH + p + t + h + ratr, random=~1|Locality, data = env_bistorta_ABU),
  lme(ll ~rv + N + CN + pH + p + t + h + ratr, random=~1|Locality, data = env_bistorta_ABU),
  ratr %~~% h,
  data = env_bistorta_ABU)
summary(bist.sem.lme_ABU6) 
#AIC: 123.451, Fisher's C = 13.451 with P-value = 0.491 and on 14 degrees of freedom



###7 community composition does not impact plants and no effect of fungi on IS
bist.sem.lme_PA7 <- psem(
  lme(asv ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ratt ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(ccpa ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(i_s ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_PA),
  lme(rv ~ N + CN + pH + p + t + asv + ratt, random=~1|Locality, data = env_bistorta_PA),
  lme(ll ~ rv + N + CN + pH + p + t + asv + ratt, random=~1|Locality, data = env_bistorta_PA),
  ccpa %~~% ratt,
  data = env_bistorta_PA)
summary(bist.sem.lme_PA7)  
#AIC: 117.973, Fisher's C = 11.973 with P-value = 0.849 and on 18 degrees of freedom

bist.sem.lme_ABU7 <- psem(
  lme(h ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ratr ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(ccabu ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(i_s ~ N + CN + pH + p + t, random=~1|Locality, data = env_bistorta_ABU),
  lme(rv ~ N + CN + pH + p + t + h + ratr, random=~1|Locality, data = env_bistorta_ABU),
  lme(ll ~rv + N + CN + pH + p + t + h + ratr, random=~1|Locality, data = env_bistorta_ABU),
  ratr %~~% h,
  data = env_bistorta_ABU)
summary(bist.sem.lme_ABU7) 
#AIC: 119.532, Fisher's C = 13.532 with P-value = 0.759 and on 18 degrees of freedom
