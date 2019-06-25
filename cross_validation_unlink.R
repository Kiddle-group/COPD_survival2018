# first extract only columns that will be used in modelling
#

#medianOffset_df <- data.frame(Var=c('age','imd','bmi','airflow','egfr'),medianOffset=NA)
#rownames(medianOffset_df) <- medianOffset_df$Var



# survival data, ids and first copd
data <- data.frame(surv_copd[,c('patid','pracid','dead5','censoredAlive','first_copd','region')])

# GENDER
data$gender <- factor(surv_copd$gender)
levels(data$gender) <- c('Male','Female')

# AIRFLOW AND FEV1 NR
data$airflow <- surv_copd$airflow
#medianOffset_df['airflow',2] <-  median(data$airflow[!surv_copd$airflowMiss])
data$airflow[!surv_copd$airflowMiss] <- data$airflow[!surv_copd$airflowMiss] - 64.6 # centre airflow
data$airflowMiss <- surv_copd$airflowMiss

# AGE (CENTRED)

#medianOffset_df['age',2] <- median(surv_copd$age_copd)
data$age <- surv_copd$age_copd - 67.7

# IMD, median impute small number of missing values
#data$imd <- surv_copd$imd2010_20
#medianOffset_df['imd',2] <- median(data$imd,na.rm=T)
#data$imd[is.na(data$imd)] <- median(data$imd,na.rm=T)
#data$imd <- data$imd - medianOffset_df['imd',2]


# SMOKING, random sample to impute small number of missing values
data$smoke <- factor(surv_copd$smokstatus,labels=c('Current','Never','Ex'))
data$smoke[is.na(data$smoke)] <- sample(data$smoke[!is.na(data$smoke)],numNA(data$smoke))

# BMI, centred and NR
data$bmi <- surv_copd$bmi
#medianOffset_df['bmi',2] <- median(data$bmi[!surv_copd$bmiMiss])
data$bmi[!surv_copd$bmiMiss] <- data$bmi[!surv_copd$bmiMiss] - 26 # centre
data$bmiMiss <- surv_copd$bmiMiss

# comorbidity column names for later use
colnames(surv_copd)[cols][1] <- 'c_ckd'
other_comorb <- c(colnames(surv_copd)[cols],'gord')
short_comorb <-  c(unlist(strsplit(colnames(surv_copd)[cols],'_'))[c(seq(from=2,by=2,to=73))],'gord')

# convert co-morbidity data into 'absent/present' factors 
for (i in 1:(length(cols)+1)){
  
  data[,short_comorb[i]] <- factor(data.frame(surv_copd)[,other_comorb[i]])
  levels(data[,short_comorb[i]]) <- c('absent','present')
  
}

# EGFR, centred and NR
#data$egfr <- surv_copd$egfr
#medianOffset_df['egfr',2] <- median(data$egfr[!surv_copd$egfrMiss])
#data$egfr[!surv_copd$egfrMiss] <- data$egfr[!surv_copd$egfrMiss] - median(data$egfr[!surv_copd$egfrMiss])
#data$egfrMiss <- surv_copd$egfrMiss

# ALBUMIN, centred and NR
#data$albumin <- surv_copd$albumin
#data$albuminMiss <- is.na(surv_copd$albumin)
#data$albumin[!data$albuminMiss] <- data$albumin[!data$albuminMiss] - median(data$albumin[!data$albuminMiss])
#data$albumin[data$albuminMiss] <- 0
#
# CRP, centred and NR
#data$crp <- surv_copd$crp
#data$crpMiss <- is.na(surv_copd$crp)
#data$crp[!data$crpMiss] <- data$crp[!data$crpMiss] - median(data$crp[!data$crpMiss])
#data$crp[data$crpMiss] <- 0

# PLATELETS, centred and NR
#data$plate <- surv_copd$plate
#data$plateMiss <- is.na(surv_copd$plate)
#data$plate[!data$plateMiss] <- data$plate[!data$plateMiss] - median(data$plate[!data$plateMiss])
#data$plate[data$plateMiss] <- 0

# CKD coding
#data$ckd_mild <- factor(surv_copd$c_ckd & !surv_copd$cci_ckd)
#levels(data$ckd_mild) <- c('absent','present')
#data$ckd_mod_severe <- factor(surv_copd$cci_ckd_copd)
#levels(data$ckd_mod_severe) <- c('absent','present')
  
# all patients have copd, used in score CMS score below
#data$copd <- 'present'





data$logits <- 1.874370047 + (-0.085820725)*(data$age) + (-0.001089517)*(data$age^2) +
  (0.036567533)*(data$bmi) + (-0.002428287)*(data$bmi^2) + (0.013684613)*(data$airflow) +
  (-0.000277955)*(data$airflow^2) +(0.264175108)*(data$gender=='Female') +
  (0.641164355)*(data$smoke=='Never') + (0.402440002)*(data$smoke=='Ex') +
  (-0.504466266)*data$bmiMiss + (-0.540682601)*data$airflowMiss + (-0.813714263)*(data$can == 'present') +
  (-0.811468347)*(data$hf == 'present') + (-0.720885687)*(data$ap == 'present') + (-0.522215569)*(data$pvd == 'present') +
  (-0.432218107)*(data$atr == 'present') + (0.232582125)*(data$ast == 'present') + (-0.390370747)*(data$dia == 'present') +
  (-0.351968752) *(data$str == 'present') + (-0.379362972) * (data$epi == 'present') + (0.292636023) *(data$ibs == 'present') +
  (-0.412296068) *(data$scz == 'present') + (-0.317463786) *(data$ops == 'present') + (-0.265865416) *(data$ibd == 'present') +
  (-0.286588747) *(data$con == 'present') + (-0.291772449) *(data$dep == 'present') + (-0.245848919) *(data$rhe == 'present')

data$prob <- (exp(data$logits)/(1+exp(data$logits)))

results <- data.frame(region=NA,brier=NA,slope=NA,auc=NA)


results[1,]<- c('All unlinked',val.prob(data$prob[!data$censoredAlive],!data$dead5[!data$censoredAlive])[c('Brier','Slope','C (ROC)')])

results[2,]<- c('England unlinked',val.prob(data$prob[!data$censoredAlive & data$region <=10],!data$dead5[!data$censoredAlive & data$region <=10])[c('Brier','Slope','C (ROC)')])
                
results[3,]<- c('Northern Ireland',val.prob(data$prob[!data$censoredAlive & data$region ==11],!data$dead5[!data$censoredAlive & data$region ==11])[c('Brier','Slope','C (ROC)')])

results[4,]<- c('Scotland',val.prob(data$prob[!data$censoredAlive & data$region ==12],!data$dead5[!data$censoredAlive & data$region ==12])[c('Brier','Slope','C (ROC)')])
                                               
results[5,]<- c('Wales',val.prob(data$prob[!data$censoredAlive & data$region ==13],!data$dead5[!data$censoredAlive & data$region ==13])[c('Brier','Slope','C (ROC)')])

write.csv(results,file='../output/unlinked.csv')
