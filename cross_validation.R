# first extract only columns that will be used in modelling
#

medianOffset_df <- data.frame(Var=c('age','imd','bmi','airflow','egfr'),medianOffset=NA)
rownames(medianOffset_df) <- medianOffset_df$Var

# survival data, ids and first copd
data <- surv_copd[,c('patid','pracid','years','censored','years5','censored5','first_copd')]
data$SurvObj <- with(data, Surv(years, censored == 0))
data$SurvObj2 <- with(data, Surv(years5, censored5 == 0))
data$dead5 <- factor(data$years < 5 & data$censored ==0)
data$surv5 <- factor(!(data$years < 5 & data$censored ==0))

# GENDER
data$gender <- factor(surv_copd$gender)
levels(data$gender) <- c('Male','Female')

# AIRFLOW AND FEV1 NR
data$airflow <- surv_copd$airflow
medianOffset_df['airflow',2] <-  median(data$airflow[!surv_copd$airflowMiss])
data$airflow[!surv_copd$airflowMiss] <- data$airflow[!surv_copd$airflowMiss] - median(data$airflow[!surv_copd$airflowMiss]) # centre airflow
data$airflowMiss <- surv_copd$airflowMiss

# AGE (CENTRED)

medianOffset_df['age',2] <- median(surv_copd$age_copd)
data$age <- surv_copd$age_copd - median(surv_copd$age_copd)

# IMD, median impute small number of missing values
data$imd <- surv_copd$imd2010_20
medianOffset_df['imd',2] <- median(data$imd,na.rm=T)
data$imd[is.na(data$imd)] <- median(data$imd,na.rm=T)
data$imd <- data$imd - medianOffset_df['imd',2]


# SMOKING, random sample to impute small number of missing values
data$smoke <- factor(surv_copd$smokstatus,labels=c('Current','Never','Ex'))
data$smoke[is.na(data$smoke)] <- sample(data$smoke[!is.na(data$smoke)],numNA(data$smoke))

# BMI, centred and NR
data$bmi <- surv_copd$bmi
medianOffset_df['bmi',2] <- median(data$bmi[!surv_copd$bmiMiss])
data$bmi[!surv_copd$bmiMiss] <- data$bmi[!surv_copd$bmiMiss] - median(data$bmi[!surv_copd$bmiMiss]) # centre
data$bmiMiss <- surv_copd$bmiMiss

# comorbidity column names for later use
colnames(surv_copd)[cols][1] <- 'c_ckd'
other_comorb <- c(colnames(surv_copd)[cols],'gord')
short_comorb <-  c(unlist(strsplit(colnames(surv_copd)[cols],'_'))[c(seq(from=2,by=2,to=73))],'gord')

# convert co-morbidity data into 'absent/present' factors 
for (i in 1:(length(cols)+1)){
  
  data[,short_comorb[i]] <- factor(surv_copd[,other_comorb[i]])
  levels(data[,short_comorb[i]]) <- c('absent','present')
  
}

# EGFR, centred and NR
data$egfr <- surv_copd$egfr
medianOffset_df['egfr',2] <- median(data$egfr[!surv_copd$egfrMiss])
data$egfr[!surv_copd$egfrMiss] <- data$egfr[!surv_copd$egfrMiss] - median(data$egfr[!surv_copd$egfrMiss])
data$egfrMiss <- surv_copd$egfrMiss

# ALBUMIN, centred and NR
data$albumin <- surv_copd$albumin
data$albuminMiss <- is.na(surv_copd$albumin)
data$albumin[!data$albuminMiss] <- data$albumin[!data$albuminMiss] - median(data$albumin[!data$albuminMiss])
data$albumin[data$albuminMiss] <- 0

# CRP, centred and NR
data$crp <- surv_copd$crp
data$crpMiss <- is.na(surv_copd$crp)
data$crp[!data$crpMiss] <- data$crp[!data$crpMiss] - median(data$crp[!data$crpMiss])
data$crp[data$crpMiss] <- 0

# PLATELETS, centred and NR
data$plate <- surv_copd$plate
data$plateMiss <- is.na(surv_copd$plate)
data$plate[!data$plateMiss] <- data$plate[!data$plateMiss] - median(data$plate[!data$plateMiss])
data$plate[data$plateMiss] <- 0

# CKD coding
data$ckd_mild <- factor(surv_copd$c_ckd & !surv_copd$cci_ckd)
levels(data$ckd_mild) <- c('absent','present')
data$ckd_mod_severe <- factor(surv_copd$cci_ckd_copd)
levels(data$ckd_mod_severe) <- c('absent','present')
  
# all patients have copd, used in score CMS score below
data$copd <- 'present'

data$silviaM <- silviaScoreM(data)
data$silviaG <- silviaScoreG(data)

# Make reproducible despite random number generation
set.seed(1)

# details of practices, as CV done on practice level
pracs <- sort(unique(data$pracid))
num_pracs <- length(pracs)

# 20% held out test data
num_test <- round(num_pracs*0.2)
test_pracs <- sample(pracs,num_test)
copd_test <- data[data$pracid %in% test_pracs,]

# 80% training data
copd_train <- data[!(data$pracid %in% test_pracs),]

# 5 iterations of 10-fold CV practice indexes
r5cv10_train <- vector('list',5)
r5cv10_test <- vector('list',5)

# single iteration of 10-fold CV practice indexes
cv_train <- vector('list',10)
cv_test <- vector('list',10)

train_pracs <- unique(copd_train$pracid)

# 5 iterations of 10-fold CV
for (r in 1:5){
  
  print(r)

  # first fold, randomly sample practice ids
  cv_test[[1]] <-  sample(train_pracs,30)
  cv_train[[1]] <- train_pracs[!(train_pracs %in% cv_test[[1]])]
  
  #ncv_size <- length(cv10_train[[1]])/10
  #size_int <- as.integer(ncv_size)
  #ncv_sizes <- c(rep(size_int,(ncv_size-size_int)*10),rep(size_int+1,10-(ncv_size-size_int)*10))
  
  # within each iteration, sample without replacement, so remove from pool
  remove <- c(which(train_pracs %in% cv_test[[1]]))
  
  # number practices to add to CV test set
  # total not neat multiple, so later this increases to 31
  to_add <- 30
  
  # for other folds
  for (i in 2:10){
    
    #print(i)
    
    # after certain point increase number to add to 31
    if (i ==7){to_add<-to_add+1}
    
    #randomly sample practice ids
    cv_test[[i]] <-  sample(train_pracs[-remove],to_add)
    cv_train[[i]] <- train_pracs[!(train_pracs %in% cv_test[[i]])]
    
    # within each iteration, sample without replacement, so remove from pool
    remove <- c(remove,which(train_pracs %in% cv_test[[i]]))

  }

  # add each iteration to list
  r5cv10_train[[r]] <- cv_train
  r5cv10_test[[r]] <- cv_test
  
  
}

rm(train_pracs,test_pracs,remove,pracs,num_test,num_pracs,i)
