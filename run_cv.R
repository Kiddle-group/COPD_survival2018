# Stratified by bins of number of comorbidities, association of not-recorded indicators different?
cols <- c(20:55)[-c(3,7,10,17,18,20,34)] # column numbers for Barnett comorbidities

# number of comorbidities per patient
copd_train$comorbid <- apply(copd_train[,cols]=='present',1,'sumBool')
copd_test$comorbid <- apply(copd_test[,cols]=='present',1,'sumBool')


# strings used to describe models later
timeUsed <- c('TTE','nonTTE')
modelApproach <- c('Saturated','Ridge','Lasso','RandomForest')
varSet <- c('Basic (B)','B + CCI','B + CMS','B + co-morbidities (C)','B+C+ C interactions (int) ')

# combining strings above to define modelling approaches
models_df <- data.frame(modelNum=1:20,vars=rep(varSet,c(2,2,2,6,8)),time=c(rep(timeUsed,times=3),rep(timeUsed,each=3),rep(timeUsed,each=4)),model=c(rep(modelApproach[1],6),rep(modelApproach[-4],2),rep(modelApproach,2)),selection=c(rep(NA,6),rep(c(NA,'C','C'),2),rep(c(NA,'int','int','ensemble'),2)))

# hard to get large survival random forests to work, so leave out
models_df <- models_df[-16,]
models_df[,1] <- 1:19

# data frame to collect results of CV
cv_df <- data.frame(reps=rep(1:5,each=10),fold=rep(1:10,times=5),matrix(NA,50,dim(models_df)[1]*4))

# data to collect from each modelling attempt
str2 <- c('brier','slope','auc','runtime')

# populate data frame column names, modelNum corresponds to model_df
colnames(cv_df)[-c(1,2)] <- paste(rep(1:dim(models_df)[1],each=4),rep(str2,dim(models_df)[1]),sep='_')

# works out the prevalence for comorbidity interactions
inter_prev <- function(data,num){as.numeric(table(data)[2]/num)}

set.seed(2)

# row number for results
k <- 0

# random forest tuning parameters in a single fold, as 3 params and long-runtime
r <- 1
i <- 1

# use list to get training folds
cv_train <- copd_train[copd_train$pracid %in% r5cv10_train[[r]][[i]],]

# use list to get test fold
cv_test <- copd_train[copd_train$pracid %in% r5cv10_test[[r]][[i]],]

# for models with interactions, generate interactions
# and appropriately centre and scale data
cv_train_obj <- R2BGLiMS_format(cv_train)
cv_test_obj <- R2BGLiMS_format(cv_test,old=cv_train_obj)#

# for those objects, column names to keep_in (basic)
# or for main or pairwise co-morbidity terms
comorb_int <- colnames(cv_train_obj$data)[c(51:514)]
comorb_main <- colnames(cv_train_obj$data)[c(20:50)]
keep_in <- colnames(cv_train_obj$data)[c(11:19)]

# egfr is not binary so treat differently to other co-morbidity variables
egfrs <- grep(pattern = "egfr",x=comorb_int)
egfr_int <- comorb_int[egfrs][seq(from=1,by = 2,length.out = 29)]

# keep only interactions that are >1% prevalent
comorb_int <- comorb_int[!(comorb_int %in% egfr_int)]
tmp <- apply(cv_train_obj$data[,comorb_int],2,'inter_prev',num=dim(cv_train_obj$data)[1])
comorb_prev <- names(which(tmp>0.01))
comorb_prev <- c(comorb_prev,egfr_int)

# calculate outcome, dead within 5 years
outcome <- dead_by(cv_test$years,cv_test$censored,5)


fit_tree500 <- randomForest(y=cv_train$surv5,x=cv_train_obj$data[,c(keep_in,comorb_main)],ntree=500)
fit_tree1000 <- randomForest(y=cv_train$surv5,x=cv_train_obj$data[,c(keep_in,comorb_main)],ntree=1000)
fit_tree2000 <- randomForest(y=cv_train$surv5,x=cv_train_obj$data[,c(keep_in,comorb_main)],ntree=2000)

fit_mtry6 <- randomForest(y=cv_train$surv5,x=cv_train_obj$data[,c(keep_in,comorb_main)],mtry = 6)
fit_mtry8 <- randomForest(y=cv_train$surv5,x=cv_train_obj$data[,c(keep_in,comorb_main)],mtry = 8)
fit_mtry10 <- randomForest(y=cv_train$surv5,x=cv_train_obj$data[,c(keep_in,comorb_main)],mtry = 10)

fit_size1 <- randomForest(y=cv_train$surv5,x=cv_train_obj$data[,c(keep_in,comorb_main)],nodesize=1)
fit_size2 <- randomForest(y=cv_train$surv5,x=cv_train_obj$data[,c(keep_in,comorb_main)],nodesize=2)
fit_size5 <- randomForest(y=cv_train$surv5,x=cv_train_obj$data[,c(keep_in,comorb_main)],nodesize=5)

forecast_tree500 <- as.numeric(predict(fit_tree500,newdata=cv_test_obj$data[,c(keep_in,comorb_main)],type='prob')[,2])
forecast_tree1000 <- as.numeric(predict(fit_tree1000,newdata=cv_test_obj$data[,c(keep_in,comorb_main)],type='prob')[,2])
forecast_tree2000 <- as.numeric(predict(fit_tree2000,newdata=cv_test_obj$data[,c(keep_in,comorb_main)],type='prob')[,2])

forecast_mtry6 <- as.numeric(predict(fit_mtry6,newdata=cv_test_obj$data[,c(keep_in,comorb_main)],type='prob')[,2])
forecast_mtry8 <- as.numeric(predict(fit_mtry8,newdata=cv_test_obj$data[,c(keep_in,comorb_main)],type='prob')[,2])
forecast_mtry10 <- as.numeric(predict(fit_mtry10,newdata=cv_test_obj$data[,c(keep_in,comorb_main)],type='prob')[,2])

forecast_size1 <- as.numeric(predict(fit_size1,newdata=cv_test_obj$data[,c(keep_in,comorb_main)],type='prob')[,2])
forecast_size2 <- as.numeric(predict(fit_size2,newdata=cv_test_obj$data[,c(keep_in,comorb_main)],type='prob')[,2])
forecast_size5 <- as.numeric(predict(fit_size5,newdata=cv_test_obj$data[,c(keep_in,comorb_main)],type='prob')[,2])

results_tree <- data.frame(value=c(500,1000,2000),brier=NA,slope=NA,auc=NA)
  
results_tree[1,c('value',str2[1:3])] <- c(500,as.numeric(val.prob(1-forecast_tree500,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))
results_tree[2,c('value',str2[1:3])] <- c(1000,as.numeric(val.prob(1-forecast_tree1000,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))
results_tree[3,c('value',str2[1:3])] <- c(2000,as.numeric(val.prob(1-forecast_tree2000,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))

results_mtry <- data.frame(value=c(6,8,10),brier=NA,slope=NA,auc=NA)

results_mtry[1,c('value',str2[1:3])] <- c(6,as.numeric(val.prob(1-forecast_mtry6,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))
results_mtry[2,c('value',str2[1:3])] <- c(8,as.numeric(val.prob(1-forecast_mtry8,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))
results_mtry[3,c('value',str2[1:3])] <- c(10,as.numeric(val.prob(1-forecast_mtry10,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))

results_size <- data.frame(value=c(1,2,5),brier=NA,slope=NA,auc=NA)

results_size[1,c('value',str2[1:3])] <- c(1,as.numeric(val.prob(1-forecast_size1,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))
results_size[2,c('value',str2[1:3])] <- c(2,as.numeric(val.prob(1-forecast_size2,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))
results_size[3,c('value',str2[1:3])] <- c(5,as.numeric(val.prob(1-forecast_size5,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))

formula <- as.formula(paste('surv5~',paste(colnames(cv_train)[c(12:19,21:58)][-c(17,14,10,27,41,24,25)],collapse='+')))
log_basic <- glm(formula,data=cv_train,family=binomial)

forecast <- predict(log_basic,type = "response",newdata=cv_test)

results_tree[4,] <- c('Logistic',as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))
results_mtry[4,] <- c('Logistic',as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))
results_size[4,] <- c('Logistic',as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]))

set.seed(2)

# 5 iterations
for (r in 1:5){

  # 10 folds
  for (i in 1:10){
    
    k <- k+1
    
    print(paste('repetition',r))
    print(paste('fold',i))
    
    # use list to get training folds
    cv_train <- copd_train[copd_train$pracid %in% r5cv10_train[[r]][[i]],]
    
    # use list to get test fold
    cv_test <- copd_train[copd_train$pracid %in% r5cv10_test[[r]][[i]],]
    
    # for models with interactions, generate interactions
    # and appropriately centre and scale data
    cv_train_obj <- R2BGLiMS_format(cv_train)
    cv_test_obj <- R2BGLiMS_format(cv_test,old=cv_train_obj)#
    
    # for those objects, column names to keep_in (basic)
    # or for main or pairwise co-morbidity terms
    comorb_int <- colnames(cv_train_obj$data)[c(51:514)]
    comorb_main <- colnames(cv_train_obj$data)[c(20:50)]
    keep_in <- colnames(cv_train_obj$data)[c(11:19)]
    
    # egfr is not binary so treat differently to other co-morbidity variables
    egfrs <- grep(pattern = "egfr",x=comorb_int)
    egfr_int <- comorb_int[egfrs][seq(from=1,by = 2,length.out = 29)]
    
    # keep only interactions that are >1% prevalent
    comorb_int <- comorb_int[!(comorb_int %in% egfr_int)]
    tmp <- apply(cv_train_obj$data[,comorb_int],2,'inter_prev',num=dim(cv_train_obj$data)[1])
    comorb_prev <- names(which(tmp>0.01))
    comorb_prev <- c(comorb_prev,egfr_int)
    
    # calculate outcome, dead within 5 years
    outcome <- dead_by(cv_test$years,cv_test$censored,5)
    
    print('Basic + CMS')
    
    # Saturated basic TTE

    a <- Sys.time()
    formula <- as.formula(SurvObj2~age+gender+imd+smoke+bmi+airflow+bmiMiss+airflowMiss)
    fit <- coxph(formula,data=cv_train)
    
    forecast <- pred_coxph(fit,cv_test)
    b <- Sys.time()
    
    cv_df[k,paste(1,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    
    # Saturated basic non-TTE

    a <- Sys.time()
    formula <- as.formula(surv5~age+gender+imd+smoke+bmi+airflow+bmiMiss+airflowMiss)
    log_basic <- glm(formula,data=cv_train,family=binomial)
    
    forecast <- predict(log_basic,type = "response",newdata=cv_test)
    b <- Sys.time()
    
    cv_df[k,paste(2,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    # Saturated mendonica TTE
    
    a <- Sys.time()
    formula <- as.formula(SurvObj2~age+gender+imd+smoke+bmi+airflow+bmiMiss+airflowMiss+silviaG)
    fit <- coxph(formula,data=cv_train)
    
    forecast <- pred_coxph(fit,cv_test)
    b <- Sys.time()
    
    cv_df[k,paste(5,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    # Saturated mendonica non-TTE
    
    a <- Sys.time()
    formula <- as.formula(surv5~age+gender+imd+smoke+bmi+airflow+bmiMiss+airflowMiss+silviaG)
    log_basic <- glm(formula,data=cv_train,family=binomial)
    
    forecast <- predict(log_basic,type = "response",newdata=cv_test)
    b <- Sys.time()
    
    cv_df[k,paste(6,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
   
    print('Comorbidity main effects')
    
    # Saturated comorb TTE
    
    a <- Sys.time()
    comorb_formula <- as.formula(paste('SurvObj2~',paste(colnames(cv_train)[c(12:19,21:58)][-c(17,14,10,27,41,24,25)],collapse='+')))
    fit <- coxph(comorb_formula,data=cv_train)
    
    forecast <- pred_coxph(fit,cv_test)
    b <- Sys.time()
    
    cv_df[k,paste(7,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    # Ridge comorb TTE
    
    a <- Sys.time()
      
    fit <- cv.glmnet(x=as.matrix(cv_train_obj$data[,c(keep_in,comorb_main)]),y = cv_train_obj$data$SurvObj2,family='cox',
                           penalty.factor=c(rep(0,length(keep_in)),rep(1,length(comorb_main))),alpha=0)
    
    # find coefficients that were chosen using penalty term selected from nested CV
    ind <- which(fit$lambda == fit$lambda.1se)
    coefs <- fit$glmnet.fit$beta[,ind+1]
    
    comorb_formula2 <- as.formula(paste('SurvObj2~',paste(paste(keep_in,collapse = '+'),paste(comorb_main,collapse = '+'),sep='+')))
    
    # need to convert to coxph format to get baseline hazard, use iterations of zero to not do any more model udpating
    new.fit <- coxph(comorb_formula2,iter.max=0,init=coefs,data=cv_train_obj$data)
    
    forecast <- pred_coxph(new.fit,cv_test_obj$data)
    b <- Sys.time()
    
    cv_df[k,paste(8,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    # Lasso comorb TTE
    
    a <- Sys.time()
    
    fit <- cv.glmnet(x=as.matrix(cv_train_obj$data[,c(keep_in,comorb_main)]),y = cv_train_obj$data$SurvObj2,family='cox',
                     penalty.factor=c(rep(0,length(keep_in)),rep(1,length(comorb_main))),alpha=1)
    
    ind <- which(fit$lambda == fit$lambda.1se)
    
    coefs <- fit$glmnet.fit$beta[,ind+1]
    
    new.fit <- coxph(comorb_formula2,iter.max=0,init=coefs,data=cv_train_obj$data)
    
    forecast <- pred_coxph(new.fit,cv_test_obj$data)
    b <- Sys.time()
    
    cv_df[k,paste(9,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    
    # Saturated comorb non-TTE
    a <- Sys.time()
    formula <- as.formula(paste('surv5~',paste(colnames(cv_train)[c(12:19,21:58)][-c(17,14,10,27,41,24,25)],collapse='+')))
    log_basic <- glm(formula,data=cv_train,family=binomial)
    
    forecast <- predict(log_basic,type = "response",newdata=cv_test)
    b <- Sys.time()
    
    cv_df[k,paste(10,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    
    # Ridge comorb non-TTE
    
    a <- Sys.time()
    
    fit <- cv.glmnet(x=as.matrix(cv_train_obj$data[,c(keep_in,comorb_main)]),y = as.numeric(cv_train$surv5)-1,family='binomial',
                     penalty.factor=c(rep(0,length(keep_in)),rep(1,length(comorb_main))),alpha=0)
    
    forecast <- predict(fit,newx=as.matrix(cv_test_obj$data[,c(keep_in,comorb_main)]),type='response')   
    b <- Sys.time()
    
    cv_df[k,paste(11,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    # Lasso comorb non-TTE
    
    a <- Sys.time()
    
    fit <- cv.glmnet(x=as.matrix(cv_train_obj$data[,c(keep_in,comorb_main)]),y = as.numeric(cv_train$surv5)-1,family='binomial',
                     penalty.factor=c(rep(0,length(keep_in)),rep(1,length(comorb_main))),alpha=1)
    
    forecast <- predict(fit,newx=as.matrix(cv_test_obj$data[,c(keep_in,comorb_main)]),type='response')   
    
    b <- Sys.time()
    
    cv_df[k,paste(12,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    print('Comorbidity interactions')
    
    
    # Saturated interaction TTE
    
    a <- Sys.time()
    formula <- as.formula(paste('SurvObj2~',paste(paste(keep_in,collapse='+'),paste(comorb_main,collapse='+'),paste(comorb_prev,collapse='+'),sep='+')))
    fit <- coxph(formula,data=cv_train_obj$data)
    
    forecast <- pred_coxph(fit,cv_test_obj$data)
    b <- Sys.time()
    
    cv_df[k,paste(13,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a)) 
    
    
    # Ridge interaction TTE
    
    a <- Sys.time()
    
    fit <- cv.glmnet(x=as.matrix(cv_train_obj$data[,c(keep_in,comorb_main,comorb_prev)]),y = cv_train_obj$data$SurvObj2,family='cox',
                     penalty.factor=c(rep(0,length(keep_in)),rep(0,length(comorb_main)),rep(1,length(comorb_prev))),alpha=0)
    
    ind <- which(fit$lambda == fit$lambda.1se)
    
    coefs <- fit$glmnet.fit$beta[,ind+1]
    
    comorb_formula3 <- as.formula(paste('SurvObj2~',paste(paste(keep_in,collapse = '+'),paste(comorb_main,collapse = '+'),paste(comorb_prev,collapse = '+'),sep='+')))
    
    new.fit <- coxph(comorb_formula3,iter.max=0,init=coefs,data=cv_train_obj$data)
    
    forecast <- pred_coxph(new.fit,cv_test_obj$data)
    b <- Sys.time()
    
    cv_df[k,paste(14,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    # Lasso interaction TTE
    
    a <- Sys.time()
    
    fit <- cv.glmnet(x=as.matrix(cv_train_obj$data[,c(keep_in,comorb_main,comorb_prev)]),y = cv_train_obj$data$SurvObj2,family='cox',
                     penalty.factor=c(rep(0,length(keep_in)),rep(0,length(comorb_main)),rep(1,length(comorb_prev))),alpha=1)
    
    ind <- which(fit$lambda == fit$lambda.1se)
    
    coefs <- fit$glmnet.fit$beta[,ind+1]
    
    new.fit <- coxph(comorb_formula3,iter.max=0,init=coefs,data=cv_train_obj$data)
    
    forecast <- pred_coxph(new.fit,cv_test_obj$data)
    b <- Sys.time()
    
    cv_df[k,paste(15,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    # Saturated interaction non-TTE
    
    a <- Sys.time()
    formula <- as.formula(paste('dead~',paste(paste(keep_in,collapse='+'),paste(comorb_main,collapse='+'),paste(comorb_prev,collapse='+'),sep='+')))
    fit <- glm(formula,data=cv_train_obj$data,family=binomial)
    
    forecast <- predict(fit,type = "response",newdata=cv_test_obj$data)
    b <- Sys.time()
    
    cv_df[k,paste(16,str2,sep='_')] <- c(as.numeric(val.prob(forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a)) 
    
    # Ridge interaction non-TTE
    
    a <- Sys.time()
    
    fit <- cv.glmnet(x=as.matrix(cv_train_obj$data[,c(keep_in,comorb_main,comorb_prev)]),y = as.numeric(cv_train$surv5)-1,family='binomial',
                     penalty.factor=c(rep(0,length(keep_in)),rep(0,length(comorb_main)),rep(1,length(comorb_prev))),alpha=0)
    
    forecast <- predict(fit,newx=as.matrix(cv_test_obj$data[,c(keep_in,comorb_main,comorb_prev)]),type='response')  
    b <- Sys.time()
    
    cv_df[k,paste(17,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    # Lasso interaction non-TTE
    
    a <- Sys.time()
    
    fit <- cv.glmnet(x=as.matrix(cv_train_obj$data[,c(keep_in,comorb_main,comorb_prev)]),y = as.numeric(cv_train$surv5)-1,family='binomial',
                     penalty.factor=c(rep(0,length(keep_in)),rep(0,length(comorb_main)),rep(1,length(comorb_prev))),alpha=1)
    
    forecast <- predict(fit,newx=as.matrix(cv_test_obj$data[,c(keep_in,comorb_main,comorb_prev)]),type='response')   
    b <- Sys.time()
    
    cv_df[k,paste(18,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    print('Random forest')
    
    # Random forest
    
    a <- Sys.time()
    
    fit <- randomForest(y=cv_train$surv5,x=cv_train_obj$data[,c(keep_in,comorb_main)])
    
    forecast <- as.numeric(predict(fit,newdata=cv_test_obj$data[,c(keep_in,comorb_main)],type='prob')[,2])
    
    b <- Sys.time()
    
    cv_df[k,paste(19,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    # save as you go along, in case of crash
    save(cv_df,file='../output/proc_data/cv_df.RData')
  
  }
  
}

# no need to above from old code, re-run as set.seed, this code can just add the cci results
#
# load('../output/proc_data/cv_df.RData')
#
# colnames(cv_df)[11:70] <- paste(rep(5:dim(models_df)[1],each=4),rep(str2,dim(models_df)[1]-4),sep='_')
#
# cv_df <- cbind(cv_df[,1:10],matrix(NA,50,8),cv_df[,11:70])
#
# colnames(cv_df)[c(11:18)] <- paste(rep(3:4,each=4),rep(str2,2),sep='_')

set.seed(2)

k <- 0

# 5 iterations
for (r in 1:5){
  
  # 10 folds
  for (i in 1:10){
    
    k <- k+1
    
    print(paste('repetition',r))
    print(paste('fold',i))
    
    # use list to get training folds
    cv_train <- copd_train[copd_train$pracid %in% r5cv10_train[[r]][[i]],]
    
    # use list to get test fold
    cv_test <- copd_train[copd_train$pracid %in% r5cv10_test[[r]][[i]],]

    
    # calculate outcome, dead within 5 years
    outcome <- dead_by(cv_test$years,cv_test$censored,5)
    
    print('Basic + CMS')
    
    
    # Saturated CCI TTE
    
    a <- Sys.time()
    formula <- as.formula(SurvObj2~age + gender+imd+smoke+bmi+airflow+bmiMiss+airflowMiss+cci)
    fit <- coxph(formula,data=cv_train)
    
    forecast <- pred_coxph(fit,cv_test)
    b <- Sys.time()
    
    cv_df[k,paste(3,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    # Saturated CCI non-TTE
    
    a <- Sys.time()
    formula <- as.formula(surv5~age + gender+imd+smoke+bmi+airflow+bmiMiss+airflowMiss+cci)
    log_basic <- glm(formula,data=cv_train,family=binomial)
    
    forecast <- predict(log_basic,type = "response",newdata=cv_test)
    b <- Sys.time()
    
    cv_df[k,paste(4,str2,sep='_')] <- c(as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')]),as.numeric(b-a))
    
    
 
  }
  
}

# correct small formatting error
cv_df[1:50,c(3:10,19:78)] <- cv_df[-1,c(3:10,19:78)]
cv_df <- cv_df[-51,]


save(cv_df,file='../output/proc_data/cv_df.RData')

# need in long format for ggplot, i.e. on model approach on one iteration of one fold CV per row
cv_df_long <- data.frame(reps=NA,fold=NA,brier=NA,slope=NA,auc=NA,runtime=NA,model=NA)

for (i in 1:19){
  
  tmp <- cv_df[,c(1:2,(1:4)+(4*i-2))]
  
  colnames(tmp)[3:6] <- str2
  
  tmp$model <- i
  
  cv_df_long <- rbind(cv_df_long,tmp)
  
}

cv_df_long <- cv_df_long[-1,]

cv_df_long$Model <- factor(cv_df_long$model)

save(cv_df,cv_df_long,file='../output/cv_df.RData')


# model refinement

k <- 0

tweak_df <- data.frame(reps=rep(1:5,each=10),fold=rep(1:10,times=5),matrix(NA,50,5*3))

colnames(tweak_df)[-c(1,2)] <- paste(rep(1:5,each=3),rep(str2[1:3],5),sep='_')

copd_train$yearsSince2004 <-as.numeric(copd_train$first_copd - as.Date('2004-01-01')) /365

set.seed(2)

for (r in 1:5){
  
  for (i in 1:10){
    
    k <- k+1
    
    print(paste('repetition',r))
    print(paste('fold',i))
    
    cv_train <- copd_train[copd_train$pracid %in% r5cv10_train[[r]][[i]],]
    
    cv_test <- copd_train[copd_train$pracid %in% r5cv10_test[[r]][[i]],]
    
    outcome <- as.numeric(dead_by(cv_test$years,cv_test$censored,5))
    
    # Saturated comorb non-TTE - plus COPD diagnosis date
    
    comorb_formula <- as.formula(paste('surv5~yearsSince2004 + ',paste(colnames(cv_train)[c(12:19,21:58)][-c(17,14,10,27,41,24,25)],collapse='+')))
    
    fit <- glm(comorb_formula,data=cv_train,family=binomial)
    
    forecast <- predict(fit,type = "response",newdata=cv_test)
    
    tweak_df[k,paste(1,str2[1:3],sep='_')] <- as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')])
    
    # Saturated comorb non-TTE - plus not recorded interactions
    
    
    comorb_formula <- as.formula(paste('surv5~egfr*bmiMiss+egfr*airflowMiss+egfr*airflowMiss*bmiMiss + bmi*egfrMiss+bmi*airflowMiss+bmi*airflowMiss*egfrMiss + airflow*bmiMiss+airflow*egfrMiss+airflow*egfrMiss*bmiMiss + ',paste(colnames(cv_train)[c(12:19,21:58)][-c(17,14,10,27,41,24,25)],collapse='+')))
  
    fit <- glm(comorb_formula,data=cv_train,family=binomial)
    
    forecast <- predict(fit,type = "response",newdata=cv_test)
    
    tweak_df[k,paste(2,str2[1:3],sep='_')] <- as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')])
    
    
    # Saturated comorb non-TTE - plus other blood tests
    
    
    comorb_formula <- as.formula(paste('surv5~albumin+albuminMiss+crp+crpMiss+plate+plateMiss+',paste(colnames(cv_train)[c(12:19,21:58)][-c(17,14,10,27,41,24,25)],collapse='+')))
    
    fit <- glm(comorb_formula,data=cv_train,family=binomial)
    
    forecast <- predict(fit,type = "response",newdata=cv_test)
    
    tweak_df[k,paste(3,str2[1:3],sep='_')] <- as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')])
    
    
    # Saturated comorb non-TTE - plus poly 2
    
    comorb_formula <- as.formula(paste('surv5~age+I(age^2)+bmi+I(bmi^2)+imd+I(imd^2)+airflow+I(airflow^2)+egfr+I(egfr^2)+',paste(colnames(cv_train)[c(12:19,21:58)][-c(2,4,5,7,17,14,10,27,41,24,25,45)],collapse='+')))
    
    fit <- glm(comorb_formula,data=cv_train,family=binomial)
    
    forecast <- predict(fit,type = "response",newdata=cv_test)
    
    tweak_df[k,paste(4,str2[1:3],sep='_')] <- as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')])
    
    # Saturated comorb non-TTE - plus poly 3
    
    comorb_formula <- as.formula(paste('surv5~age+I(age^2)+I(age^3)+bmi+I(bmi^2)+I(bmi^3)+imd+I(imd^2)+I(imd^3)+airflow+I(airflow^2)+I(airflow^3)+egfr+I(egfr^2)+I(egfr^3)+',paste(colnames(cv_train)[c(12:19,21:58)][-c(2,4,5,7,17,14,10,27,41,24,25,45)],collapse='+')))
    
    fit <- glm(comorb_formula,data=cv_train,family=binomial)
    
    forecast <- predict(fit,type = "response",newdata=cv_test)
    
    tweak_df[k,paste(5,str2[1:3],sep='_')] <- as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')])
    
  }
  
}

set.seed(3)

# extra refinement

k <- 0

tweak_df2 <- data.frame(reps=rep(1:5,each=10),fold=rep(1:10,times=5),matrix(NA,50,13*3))

colnames(tweak_df2)[-c(1,2)] <- paste(rep(6:13,each=3),rep(str2[1:3],5),sep='_')

for (r in 1:5){
  
  for (i in 1:10){
    
    k <- k+1
    
    print(paste('repetition',r))
    print(paste('fold',i))
    
    cv_train <- copd_train[copd_train$pracid %in% r5cv10_train[[r]][[i]],]
    
    cv_test <- copd_train[copd_train$pracid %in% r5cv10_test[[r]][[i]],]
    
    outcome <- as.numeric(dead_by(cv_test$years,cv_test$censored,5))
    
    
    # Saturated comorb non-TTE -  0 interactions
    comorb_formula2 <- as.formula(paste('surv5~age+bmi+imd+airflow+egfr+I(egfrMiss*(comorbid==0))+I(bmiMiss*(comorbid==0))+I(airflowMiss*(comorbid==0))+',paste(colnames(cv_train)[c(12:19,21:58)][-c(2,4,5,7,17,14,10,27,41,24,25,45)],collapse='+')))
    

    fit <- glm(comorb_formula2,data=cv_train,family=binomial)
    
    forecast <- predict(fit,type = "response",newdata=cv_test)
    
    tweak_df2[k,paste(6,str2[1:3],sep='_')] <- as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')])
    
    
    # Saturated comorb non-TTE -  quadratic + 0 interactions
    comorb_formula3 <- as.formula(paste('surv5~age+I(age^2)+bmi+I(bmi^2)+imd+I(imd^2)+airflow+I(airflow^2)+egfr+I(egfr^2)+I(egfrMiss*(comorbid==0))+I(bmiMiss*(comorbid==0))+I(airflowMiss*(comorbid==0))+',paste(colnames(cv_train)[c(12:19,21:58)][-c(2,4,5,7,17,14,10,27,41,24,25,45)],collapse='+')))
    
    
    fit <- glm(comorb_formula3,data=cv_train,family=binomial)
    
    forecast <- predict(fit,type = "response",newdata=cv_test)
    
    tweak_df2[k,paste(7,str2[1:3],sep='_')] <- as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')])
    

    
    
  }
  
}

set.seed(4)

# supp table 2 models

k <- 0



for (r in 1:5){
  
  for (i in 1:10){
    
    k <- k+1
    
    print(paste('repetition',r))
    print(paste('fold',i))
    
    cv_train <- copd_train[copd_train$pracid %in% r5cv10_train[[r]][[i]],]
    
    cv_test <- copd_train[copd_train$pracid %in% r5cv10_test[[r]][[i]],]
    
    outcome <- as.numeric(dead_by(cv_test$years,cv_test$censored,5))
    
    # age,smoking,gender,airflow,bmi
    
    formula <- as.formula(surv5~age+I(age^2)+bmi+I(bmi^2)+airflow+I(airflow^2)+gender+smoke+bmiMiss+airflow+airflowMiss)
    log_basic <- glm(formula,data=cv_train,family=binomial)
    
    forecast <- predict(log_basic,type = "response",newdata=cv_test)
    
    tweak_df2[k,paste(8,str2[1:3],sep='_')] <- as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')])
    
  
   
    
    # age,smoking,gender,airflow,bmi + CCI
    
    formula <- as.formula(surv5~age+I(age^2)+bmi+I(bmi^2)+airflow+I(airflow^2)+gender+smoke+bmiMiss+airflow+airflowMiss + cci)
    

    fit <- glm(formula,data=cv_train,family=binomial)
    
    forecast <- predict(fit,type = "response",newdata=cv_test)
    
    tweak_df2[k,paste(9,str2[1:3],sep='_')] <- as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')])
    
    
    # age,smoking,gender,airflow,bmi + our 16
    
    formula <- as.formula(surv5~age+I(age^2)+bmi+I(bmi^2)+airflow+I(airflow^2)+gender+smoke+bmiMiss+airflowMiss+can+hf+ap+pvd+atr+ast+dia+str+epi+ibs+scz+ops+ibd+con+dep+rhe)
    
    fit <- glm(formula,data=cv_train,family=binomial)
    
    forecast <- predict(fit,type = "response",newdata=cv_test)
    
    tweak_df2[k,paste(10,str2[1:3],sep='_')] <- as.numeric(val.prob(1-forecast,outcome,pl=F)[c('Brier','Slope','C (ROC)')])
    
  
     
  }
  
}


# again, need in long format for plotting
cv_df_long_2 <- data.frame(reps=NA,fold=NA,brier=NA,slope=NA,auc=NA,model=NA)

tmp <- cv_df[,c(1:2,(1:3)+(4*2-2))]

colnames(tmp)[3:5] <- str2[1:3]

tmp$model <- 1

cv_df_long_2 <- rbind(cv_df_long_2,tmp)

tmp <- cv_df[,c(1:2,(1:3)+(4*4-2))]

colnames(tmp)[3:5] <- str2[1:3]

tmp$model <- 2

cv_df_long_2 <- rbind(cv_df_long_2,tmp)

tmp <- cv_df[,c(1:2,(1:3)+(4*6-2))]

colnames(tmp)[3:5] <- str2[1:3]

tmp$model <- 3

cv_df_long_2 <- rbind(cv_df_long_2,tmp)

tmp <- cv_df[,c(1:2,(1:3)+(4*10-2))]

colnames(tmp)[3:5] <- str2[1:3]

tmp$model <- 4

cv_df_long_2 <- rbind(cv_df_long_2,tmp)

for (i in 1:5){
  
  tmp <- tweak_df[,c(1:2,(1:3)+(3*i-1))]
  
  colnames(tmp)[3:5] <- str2[1:3]
  
  tmp$model <- i +4
  
  cv_df_long_2 <- rbind(cv_df_long_2,tmp)
  
}

for (i in 6:10){

  tmp <- tweak_df2[,c(1:2,(1:3)+3*(i-5)-1)]

  colnames(tmp)[3:5] <- str2[1:3]

  tmp$model <- i +4

  cv_df_long_2 <- rbind(cv_df_long_2,tmp)

}


cv_df_long_2 <- cv_df_long_2[-1,]

cv_df_long_short <- cv_df_long_2[cv_df_long_2$model %in% c(5,7:9),]

cv_df_long_short$model[cv_df_long_short$model == 5] <- 20
cv_df_long_short$model[cv_df_long_short$model == 7] <- 21
cv_df_long_short$model[cv_df_long_short$model == 8] <- 22
cv_df_long_short$model[cv_df_long_short$model == 9] <- 23

cv_df_supp1 <- rbind(cv_df_long[,-c(6,8)],cv_df_long_short)

# labels for plot
vars <- c('B','B + CCI','B + CMS','B + C','B + C^2','All')
models <- c('','ridge','lasso','random forest')

models_df <- data.frame(Model_number=1:19,Variables=c(rep(vars[1:3],each=2),rep(vars[4:5],each=6), vars[6]),
                        Penalty=c(rep(models[1],6),rep(models[1:3],times=4),''),
                        Model=c(rep(c('Cox','Logistic'),times=3),rep(c('Cox','Logistic'),times=2,each=3),'Random forest'))

models_df$Approach <- c(as.character(models_df$Model[1:7]),paste(models_df$Model,models_df$Penalty)[8:9],'Logistic',
                        paste(models_df$Model,models_df$Penalty)[11:12],'Cox',
                        paste(models_df$Model,models_df$Penalty)[14:15],'Logistic',
                        paste(models_df$Model,models_df$Penalty)[17:18],'Random forest')

extra_models_df <- models_df[c(10,10,10,10),]
extra_models_df$Model_number <- 20:23
extra_models_df$Approach <- c('+ diagnosis year','+albumin+CRP+platelets','+ quadratic','+ cubic')

models_df <- rbind(models_df,extra_models_df)

models_df$Approach <- factor(models_df$Approach,levels = levels(factor(models_df$Approach))[c(5,8,2,4,3,1,6,9,7,10,11)])
models_df$Variables <- factor(models_df$Variables,levels = levels(factor(models_df$Variables))[c(2,5,6,3,4,1)])


colnames(cv_df_supp1)[6] <- 'Model_number'

long_plus <- merge(cv_df_supp1,models_df,by='Model_number')

long_plus_dt <- data.table(long_plus)

write.csv(long_plus,file='../output/supp_fig1_data.csv')

# make plots, error msg is not a problem, but have to prevent it stopping script

a <- tryCatch(
  {
    ggsave(filename="../output/supp1_brier.pdf", plot=ggplot(long_plus,aes(x=Approach,y=brier,fill=Approach)) + geom_boxplot() +facet_grid(~Variables,scales='free_x',space='free_x') + 
             theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
             #gridExtra::grid.arrange(last_plot(), top="Variables and interactions considered ") +
             guides(fill=FALSE) + ylab('Brier score') +
             geom_hline(yintercept=min(long_plus_dt[,list(median=median(brier)),by=Model_number][,2]),color="red", linetype="dashed") +
             gridExtra::grid.arrange(last_plot(), top="Variables and/or interactions considered")
           ,width=6,height = 4)
  }, warning = function(w) {}, error = function(e) {}, finally = {}
)

a <- tryCatch(
  {
    ggsave(filename="../output/supp1_slope.pdf", plot=ggplot(long_plus,aes(x=Approach,y=slope,fill=Approach)) + geom_boxplot() +facet_grid(~Variables,scales='free_x',space='free_x') + 
             theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
             #gridExtra::grid.arrange(last_plot(), top="Variables and/or interactions considered ") +
             guides(fill=FALSE) + ylab('Calibration slope') +
             geom_hline(yintercept=1,color="red", linetype="dashed") +
             gridExtra::grid.arrange(last_plot(), top="Variables and/or interactions considered")
           ,width=6,height = 4)
  }, warning = function(w) {}, error = function(e) {}, finally = {}
)

a <- tryCatch(
  {
    ggsave(filename="../output/supp1_auc.pdf", plot=ggplot(long_plus,aes(x=Approach,y=auc,fill=Approach)) + geom_boxplot() +facet_grid(~Variables,scales='free_x',space='free_x') + 
             theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
             #gridExtra::grid.arrange(last_plot(), top="Variables and/or interactions considered ") +
             guides(fill=FALSE) + ylab('Area Under Curve (c-index)') +
             geom_hline(yintercept=max(long_plus_dt[,list(median=median(auc)),by=Model_number][,2]),color="red", linetype="dashed") +
             gridExtra::grid.arrange(last_plot(), top="Variables and/or interactions considered") 
           ,width=6,height = 4)
  }, warning = function(w) {}, error = function(e) {}, finally = {}
)





# now lets make supp fig 2

cv_df_long_short_2 <- cv_df_long_2[cv_df_long_2$model %in% c(12:14),]

tmp <- cbind(tweak_df[,1:2],tweak_df[,paste(4,str2[1:3],sep='_')])

colnames(tmp)[3:5] <- str2[1:3]

tmp$model <- 1


cv_df_supp2 <- rbind(tmp,cv_df_long_short_2)


# plot model labels
cv_df_supp2$Model[cv_df_supp2$model == 12] <- 'Basic - IMD (B-I)'
cv_df_supp2$Model[cv_df_supp2$model == 13] <- '(B-I) + CCI'
cv_df_supp2$Model[cv_df_supp2$model == 14] <- '21 variable model'
cv_df_supp2$Model[cv_df_supp2$model == 1] <- '36 variable model'

cv_df_supp2$Model <- factor(cv_df_supp2$Model,levels=levels(factor(cv_df_supp2$Model))[c(4,1,2,3)])

supp2_dt <- data.table(cv_df_supp2)

write.csv(cv_df_supp2,file='../output/supp_fig2_data.csv')

# make plots

a <- tryCatch(
  {
    ggsave(filename="../output/supp2_brier.pdf", plot=
             ggplot(cv_df_supp2,aes(x=Model,y=brier,fill=Model)) + geom_boxplot()  +
             guides(fill=FALSE) + ylab('Brier score') + xlab('Logistic regression model') +
             geom_hline(yintercept= min(supp2_dt[,list(median=median(brier)),by=model][,2]),color="red", linetype="dashed") 
           ,width=6,height = 4)
  }, warning = function(w) {}, error = function(e) {}, finally = {}
)



a <- tryCatch(
  {
    ggsave(filename="../output/supp2_auc.pdf", plot=
             ggplot(cv_df_supp2,aes(x=Model,y=auc,fill=Model)) + geom_boxplot()  +
             guides(fill=FALSE) + ylab('Area Under Curve (c-index)') + xlab('Logistic regression model') +
             geom_hline(yintercept= max(supp2_dt[,list(median=median(auc)),by=model][,2]),color="red", linetype="dashed") 
           ,width=6,height = 4)
    
  }, warning = function(w) {}, error = function(e) {}, finally = {}
)



a <- tryCatch(
  {
    ggsave(filename="../output/supp2_slope.pdf", plot=
             ggplot(cv_df_supp2,aes(x=Model,y=slope,fill=Model)) + geom_boxplot()  +
             guides(fill=FALSE) + ylab('Calibration slope') + xlab('Logistic regression model') +
             geom_hline(yintercept= 1,color="red", linetype="dashed") 
           ,width=6,height = 4)
  }, warning = function(w) {}, error = function(e) {}, finally = {}
)

comorb_formula <- as.formula(paste('surv5~age+I(age^2)+bmi+I(bmi^2)+imd+I(imd^2)+airflow+I(airflow^2)+egfr+I(egfr^2)+',paste(colnames(cv_train)[c(12:19,21:58)][-c(2,4,5,7,17,14,10,27,41,24,25,45)],collapse='+')))



fit <- glm(comorb_formula,data=copd_train,family=binomial)

save(fit,file='../output/final.RData')

coef_df <- data.frame(logOddsRatio=fit$coefficients,confint.default(fit))

#coef_df <- coef_df[order(coef_df$logOddsRatio),]

coef_df$label <- rownames(coef_df)

colnames(coef_df)[2:3] <- c('log_lower','log_upper')

coef_df$oddsRatio <- exp(coef_df$logOddsRatio)
coef_df$upper <- exp(coef_df$log_upper)
coef_df$lower <- exp(coef_df$log_lower)

coef_df2 <- coef_df[-c(1:11),]

percentile_df <- data.frame(Var=c('age','imd','bmi','airflow','egfr'),upperQuartile=NA,lowerQuartile=NA,medianOffset=medianOffset_df$medianOffset)
rownames(percentile_df) <- percentile_df$Var

for (i in as.character(percentile_df$Var)){
  
  if( i %in% c('age','imd')){
    
    percentile_df[i,'upperQuartile'] <-  quantile(copd_train[,i],0.75) 
    percentile_df[i,'lowerQuartile'] <-  quantile(copd_train[,i],0.25)     
    
  } else {
    
    miss_str <- paste(i,'Miss',sep='')
    
    percentile_df[i,'upperQuartile'] <-  quantile(copd_train[,i][!copd_train[,miss_str]],0.75)
    percentile_df[i,'lowerQuartile'] <-  quantile(copd_train[,i][!copd_train[,miss_str]],0.25)  
    
  }
  
  percentile_df[i,'uq_lor'] <- coef_df[i,'logOddsRatio']*percentile_df[i,'upperQuartile'] + coef_df[paste('I(',i,'^2)',sep=''),'logOddsRatio']*percentile_df[i,'upperQuartile']^2
  percentile_df[i,'uq_lower'] <- coef_df[i,'log_lower']*percentile_df[i,'upperQuartile'] + coef_df[paste('I(',i,'^2)',sep=''),'log_lower']*percentile_df[i,'upperQuartile']^2
  percentile_df[i,'uq_upper'] <- coef_df[i,'log_upper']*percentile_df[i,'upperQuartile'] + coef_df[paste('I(',i,'^2)',sep=''),'log_upper']*percentile_df[i,'upperQuartile']^2
  
  percentile_df[i,'lq_lor'] <- coef_df[i,'logOddsRatio']*percentile_df[i,'lowerQuartile'] + coef_df[paste('I(',i,'^2)',sep=''),'logOddsRatio']*percentile_df[i,'lowerQuartile']^2
  percentile_df[i,'lq_lower'] <- coef_df[i,'log_lower']*percentile_df[i,'lowerQuartile'] + coef_df[paste('I(',i,'^2)',sep=''),'log_lower']*percentile_df[i,'lowerQuartile']^2
  percentile_df[i,'lq_upper'] <- coef_df[i,'log_upper']*percentile_df[i,'lowerQuartile'] + coef_df[paste('I(',i,'^2)',sep=''),'log_upper']*percentile_df[i,'lowerQuartile']^2
  
  
  
}


upperQuart <- exp(percentile_df[,5:7])
rownames(upperQuart) <- paste(rownames(upperQuart),round(percentile_df$medianOffset+percentile_df$upperQuartile),sep='')
colnames(upperQuart) <- colnames(coef_df)[5:7]
lowerQuart <- exp(percentile_df[,8:10])
rownames(lowerQuart) <- paste(rownames(lowerQuart),round(percentile_df$medianOffset+percentile_df$lowerQuartile),sep='')
colnames(lowerQuart) <- colnames(coef_df)[5:7]

contOR <- rbind(upperQuart,lowerQuart)

plotORs <- rbind(contOR,coef_df2[,5:7])

wtf <- c('Five_Year_Mortality','Chronic Kidney Disease','Alcohol problems','Anorexia/bulimia','Atrial fibrillation','Blindness/low vision','Bronciectasis','Chronic liver disease','Sinusitis','Coronary heart disease','Dementia','Diabetes','Diverticular disease','Hearing loss','Heart failure','Hypertension',
'Inflammatory bowel disease','Learning disability','Multiple sclerosis','Peripheral vascular disorder','Parkinsons','Prostate disorders','Substance abuse','Connective tissue disorders',
'Stroke','Thyroid disorders','Anxiety','Asthma','Cancer','Constipation','Depression','Epilepsy','Painful condition','Irritable bowel syndrome','Migrane','Psoriasis/eczema','Pyschosis/bipolar','Normalised_MM_burden')

wtf <- sort(wtf)

plotORs <- plotORs[order(plotORs$oddsRatio),]

plotORs$label <- rev(c('Age 59 vs 68','Never smoker','Ex smoker','Female',wtf[c(24,4)],wtf[32],'FEV1 % predicted 79% vs 65%',wtf[17],wtf[20],'Body Mass Index 30 vs 26',wtf[35],'Index of Multiple Deprivation 6 vs 12',
        wtf[22],'eGFR 88 vs 74',wtf[38],'GORD','Index of Multiple Deprivation 16 vs 12',wtf[13],'eGFR 60 vs 74',wtf[c(33,7)],'eGFR NR*',wtf[6],'Body Mass Index 23 vs 26',wtf[c(29,3,11,15,12)],'FEV1 % predicted 50% vs 65%',wtf[c(23,37,18,36,16,34,5,31)],'Body Mass Index NR','FEV1 NR',wtf[c(1,21)],'Age 76 vs 68',wtf[8]))

plotORs$label <- factor(plotORs$label,levels=unique(plotORs$label))


ggsave(filename="../output/finalModel.pdf",ggplot(data=plotORs,aes(x=label,y=oddsRatio,ymin=lower,ymax=upper)) + geom_pointrange() + geom_hline(yintercept=1,lty=2) + coord_flip() + xlab('Variable') + ylab('Survival odds ratio (95% CI)') + theme_bw() +ylim(0.35,2.05),width=5.4,height=8)


coef_df$label <- c('Intercept','Age in years from 67.7','Age in years from 67.7, squared','Body Mass Index from 26','Body Mass Index from 26, squared',
                   'Index of Multiple Deprivation in twentiles from 12','Index of Multiple Deprivation in twentiles from 12, squared',
                   'Airflow obstruction percent from 64.6%','Airflow obstruction percent from 64.6%, squared',
                   'eGFR from 74','eGFR from 74, squared','Female','FEV1 NR','Never smoker','Ex smoker','Body Mass Index NR',
                   wtf[c(1,5,6,7,35,13,16,17,20,21,22,23,31,32,37,11,36,38,3,4,8,12,15,18,29,24,33,34)],'GORD','eGFR not recorded twice')

write.csv(coef_df,file='../output/finalModel.csv')


# Stratified by bins of number of comorbidities, association of not-recorded indicators different?
cols <- c(20:55)[-c(3,7,10,17,18,20,34)] # column numbers for Barnett comorbidities

# number of comorbidities per patient
copd_train$comorbid <- apply(copd_train[,cols]=='present',1,'sumBool')
copd_test$comorbid <- apply(copd_test[,cols]=='present',1,'sumBool')

# bin patients based on that number
none <- copd_train$comorbid==0
onetwo <- copd_train$comorbid>0 & copd_train$comorbid<=2
threemore <- copd_train$comorbid>=3

none_fit <- glm(surv5~age + I(age^2) + bmi + I(bmi^2) + imd + I(imd^2) + airflow + I(airflow^2) + egfr + I(egfr^2) + gender + airflowMiss + 
                  smoke + bmiMiss + egfrMiss,data=copd_train[none,],family=binomial)

onetwo_fit <- glm(surv5~age + I(age^2) + bmi + I(bmi^2) + imd + I(imd^2) + airflow + I(airflow^2) + egfr + I(egfr^2) + gender + airflowMiss + 
                    smoke + bmiMiss + egfrMiss,data=copd_train[onetwo,],family=binomial)

threemore_fit <- glm(surv5~age + I(age^2) + bmi + I(bmi^2) + imd + I(imd^2) + airflow + I(airflow^2) + egfr + I(egfr^2) + gender + airflowMiss + 
                       smoke + bmiMiss + egfrMiss,data=copd_train[threemore,],family=binomial)




airflowMiss_coef_df <- data.frame(label = c('all_patients','0','1 - 2','>2'),logOddsRatio=c(fit$coefficients['airflowMissTRUE'],none_fit$coefficients['airflowMissTRUE'],onetwo_fit$coefficients['airflowMissTRUE'],threemore_fit$coefficients['airflowMissTRUE']),
                                  rbind(confint.default(fit)['airflowMissTRUE',],
                                        confint.default(none_fit)['airflowMissTRUE',],
                                        confint.default(onetwo_fit)['airflowMissTRUE',],
                                        confint.default(threemore_fit)['airflowMissTRUE',]))

colnames(airflowMiss_coef_df)[3:4] <- c('log_lower','log_upper')

airflowMiss_coef_df$oddsRatio <- exp(airflowMiss_coef_df$logOddsRatio)
airflowMiss_coef_df$upper <- exp(airflowMiss_coef_df$log_upper)
airflowMiss_coef_df$lower <- exp(airflowMiss_coef_df$log_lower)  

airflowMiss_coef_df$label <- factor(airflowMiss_coef_df$label,levels=rev(c('all_patients','0','1 - 2','>2')))

ggsave(filename="../output/airflowMiss.pdf",ggplot(data=airflowMiss_coef_df,aes(x=label,y=oddsRatio,ymin=lower,ymax=upper)) + geom_pointrange() + geom_hline(yintercept=1,lty=2) + coord_flip() + xlab('Variable') + ylab('Survival odds ratio (95% CI)') + theme_bw() +ylim(0.35,1.01),width=4,height=4)


bmiMiss_coef_df <- data.frame(label = c('all_patients','0','1 - 2','>2'),logOddsRatio=c(fit$coefficients['bmiMissTRUE'],none_fit$coefficients['bmiMissTRUE'],onetwo_fit$coefficients['bmiMissTRUE'],threemore_fit$coefficients['bmiMissTRUE']),
                                  rbind(confint.default(fit)['bmiMissTRUE',],
                                        confint.default(none_fit)['bmiMissTRUE',],
                                        confint.default(onetwo_fit)['bmiMissTRUE',],
                                        confint.default(threemore_fit)['bmiMissTRUE',]))

colnames(bmiMiss_coef_df)[3:4] <- c('log_lower','log_upper')

bmiMiss_coef_df$oddsRatio <- exp(bmiMiss_coef_df$logOddsRatio)
bmiMiss_coef_df$upper <- exp(bmiMiss_coef_df$log_upper)
bmiMiss_coef_df$lower <- exp(bmiMiss_coef_df$log_lower)  

bmiMiss_coef_df$label <- factor(bmiMiss_coef_df$label,levels=rev(c('all_patients','0','1 - 2','>2')))

ggsave(filename="../output/bmiMiss.pdf",ggplot(data=bmiMiss_coef_df,aes(x=label,y=oddsRatio,ymin=lower,ymax=upper)) + geom_pointrange() + geom_hline(yintercept=1,lty=2) + coord_flip() + xlab('Variable') + ylab('Survival odds ratio (95% CI)') + theme_bw() +ylim(0.35,1.01),width=4,height=4)

egfrMiss_coef_df <- data.frame(label = c('all_patients','0','1 - 2','>2'),logOddsRatio=c(fit$coefficients['egfrMissTRUE'],none_fit$coefficients['egfrMissTRUE'],onetwo_fit$coefficients['egfrMissTRUE'],threemore_fit$coefficients['egfrMissTRUE']),
                              rbind(confint.default(fit)['egfrMissTRUE',],
                                    confint.default(none_fit)['egfrMissTRUE',],
                                    confint.default(onetwo_fit)['egfrMissTRUE',],
                                    confint.default(threemore_fit)['egfrMissTRUE',]))

colnames(egfrMiss_coef_df)[3:4] <- c('log_lower','log_upper')

egfrMiss_coef_df$oddsRatio <- exp(egfrMiss_coef_df$logOddsRatio)
egfrMiss_coef_df$upper <- exp(egfrMiss_coef_df$log_upper)
egfrMiss_coef_df$lower <- exp(egfrMiss_coef_df$log_lower)  

egfrMiss_coef_df$label <- factor(egfrMiss_coef_df$label,levels=rev(c('all_patients','0','1 - 2','>2')))

ggsave(filename="../output/egfrMiss.pdf",ggplot(data=egfrMiss_coef_df,aes(x=label,y=oddsRatio,ymin=lower,ymax=upper)) + geom_pointrange() + geom_hline(yintercept=1,lty=2) + coord_flip() + xlab('Variable') + ylab('Survival odds ratio (95% CI)') + theme_bw() +ylim(0.35,1.31),width=4,height=4)


# scores for final model in train set
forecast <- predict(fit,newdata = copd_train,type='response')
as.numeric(val.prob(forecast,copd_train$surv5=='TRUE',pl=F)[c('Brier','Slope','C (ROC)')])
#[1] 0.1386016 1.0000000 0.7998773

# test final model!!
forecast <- predict(fit,newdata = copd_test,type='response')
as.numeric(val.prob(forecast,copd_test$surv5=='TRUE',pl=F)[c('Brier','Slope','C (ROC)')])
#[1] 0.1377851 0.9913257 0.8046013





# Sensitivity analysis - how well in those with FEV1 and FEV1/FVC < 0.7?
fvc <- read.dta13('../data/steve_fvc.dta')
fvc <- merge(fvc,all_copd,by='patid')
fvc_before <- fvc[fvc[,'first_copd'] - fvc[,'eventdate'] >= 0,]
fvc_before <- fvc_before[!is.na(fvc_before$fvc),]

fvc_before$post_diagnosis <- dateNumeric(fvc_before$eventdate) - dateNumeric(fvc_before$first_copd)
fvc_before <- fvc_before[(fvc_before$post_diagnosis <= 0) & (fvc_before$post_diagnosis >= -365),]

# reorder by patient, days before diagnosis and highest reading for the day
fvc_before <- fvc_before[order(fvc_before$patid,fvc_before$post_diagnosis,fvc_before$fvc,decreasing = c(F,T,T)),]

# nice trick to give per patient order
fvc_before$order <- ave(rep(1,nrow(fvc_before)),fvc_before[,1],FUN=seq_along )

# extract latest observation
fvc_before <- fvc_before[fvc_before$order==1,]

# merge with raw fev1 data

#FEV1 in last 12 months, remove missing FEV1 and deal with units
fev1 <- read.dta13('../data/all_fev1.dta') # data from Quint group
fev1_nona <- fev1[!is.na(fev1$fev1),]
fev1_nona$fev1[fev1_nona$data3 == 89] <- fev1_nona$fev1[fev1_nona$data3 == 89]/1000
fev1_nona <- fev1_nona[fev1_nona$fev1 < 5,] # remove outlying values

# extract just values from eligible patients within year before COPD diagnosis
fev1_eligible <- fev1_nona
fev1_eligible$post_diagnosis <- dateNumeric(fev1_eligible$eventdate) - dateNumeric(fev1_eligible$first_copd)
fev1_before <- fev1_eligible[(fev1_eligible$post_diagnosis <= 0) & (fev1_eligible$post_diagnosis >= -365),]

# reorder by patient, days before diagnosis and highest reading for the day
fev1_before <- fev1_before[order(fev1_before$patid,fev1_before$post_diagnosis,fev1_before$fev1,decreasing = c(F,T,T)),]

# nice trick to give per patient order
fev1_before$order <- ave(rep(1,nrow(fev1_before)),fev1_before[,1],FUN=seq_along )

# extract latest observation
fev1_before <- fev1_before[fev1_before$order==1,]

copd_test <- merge(copd_test,fev1_before[,c('patid','fev1')],by='patid',all.x=T)
copd_test <- merge(copd_test,fvc_before[,c('patid','fvc')],by='patid',all.x=T)

a <- copd_test$fev1 / (copd_test$fvc/1000)
definite_copd <- which(a < 0.7)

#performance in definite COPD

as.numeric(val.prob(forecast,copd_test$surv5=='TRUE',pl=F)[c('Brier','Slope','C (ROC)')])
#[1] 0.1376445 0.9933361 0.8050062

as.numeric(val.prob(forecast[definite_copd],copd_test$surv5[definite_copd]=='TRUE',pl=F)[c('Brier','Slope','C (ROC)')])
#[1] 0.1329983 0.9790026 0.7916849






# PPV NPV
pred <- prediction(predictions = 1-forecast,labels = copd_test$surv5=='FALSE')
ppv <- performance(pred,'ppv')
npv <- performance(pred,'npv')

ppv_npv_df <- data.frame(Probability_cutoff=c(npv@x.values[[1]],npv@x.values[[1]]),Proportion=c(ppv@y.values[[1]],npv@y.values[[1]]))

ppv_npv_df$Measure[1:12097] <- 'PPV'
ppv_npv_df$Measure[-(1:12097)] <- 'NPV'

# downsample to 100 points, excessive to share greater detail on github

ind <- round(seq(from = 2, to = 24192, length.out = 200))

write.csv(ppv_npv_df[ind,],file='../output/36var_PPV_NPV.csv')

ggsave(filename='../output/ppv_npv.pdf',ggplot(data=ppv_npv_df,
       aes(x=Probability_cutoff, y=Proportion, colour=Measure)) +
  geom_line() +theme_bw(),height=4,width=4)


# risk quartiles
h1 <- hoslem.test(copd_test$surv5==TRUE,forecast,g=5)


write.csv(cbind(h1$observed,h1$expected),file='../output/36var_quintiles.csv')


# simpler model


formula <- as.formula(surv5~age+I(age^2)+bmi+I(bmi^2)+airflow+I(airflow^2)+gender+smoke+bmiMiss+airflowMiss+can+hf+ap+pvd+atr+ast+dia+str+epi+ibs+scz+ops+ibd+con+dep+rhe)

fit <- glm(formula,data=copd_train,family=binomial)

# scores for final model in train set
forecast <- predict(fit,newdata = copd_train,type='response')
as.numeric(val.prob(forecast,copd_train$surv5=='TRUE',pl=F)[c('Brier','Slope','C (ROC)')])


# test final model!!
forecast <- predict(fit,newdata = copd_test,type='response')
as.numeric(val.prob(forecast,copd_test$surv5=='TRUE',pl=F)[c('Brier','Slope','C (ROC)')])

as.numeric(val.prob(forecast[definite_copd],copd_test$surv5[definite_copd]=='TRUE',pl=F)[c('Brier','Slope','C (ROC)')])
#0.1335123 0.9872784 0.7906763

# PPV NPV
pred <- prediction(predictions = 1-forecast,labels = copd_test$surv5=='FALSE')
ppv <- performance(pred,'ppv')
npv <- performance(pred,'npv')

ppv_npv_df <- data.frame(Probability_cutoff=c(npv@x.values[[1]],npv@x.values[[1]]),Proportion=c(ppv@y.values[[1]],npv@y.values[[1]]))

ppv_npv_df$Measure[1:12094] <- 'PPV'
ppv_npv_df$Measure[-(1:12094)] <- 'NPV'

# downsample to 100 points, excessive to share greater detail on github

ind <- round(seq(from = 2, to = 24188, length.out = 200))

write.csv(ppv_npv_df[ind,],file='../output/21var_PPV_NPV.csv')

ggsave(filename='../output/ppv_npv_simple.pdf',ggplot(data=ppv_npv_df,
                                               aes(x=Probability_cutoff, y=Proportion, colour=Measure)) +
         geom_line() +theme_bw(),height=4,width=4)


# risk quartiles
h1 <- hoslem.test(copd_test$surv5==TRUE,forecast,g=5)



write.csv(cbind(h1$observed,h1$expected),file='../output/21var_quintiles.csv')


coef_df <- data.frame(logOddsRatio=fit$coefficients,confint.default(fit))

#coef_df <- coef_df[order(coef_df$logOddsRatio),]

coef_df$label <- rownames(coef_df)

colnames(coef_df)[2:3] <- c('log_lower','log_upper')

coef_df$oddsRatio <- exp(coef_df$logOddsRatio)
coef_df$upper <- exp(coef_df$log_upper)
coef_df$lower <- exp(coef_df$log_lower)

write.csv(coef_df,file='../output/21var_model.csv')


