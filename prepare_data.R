# Script to load and merge relevant COPD CPRD data
# Copyright Steven Kiddle 2018

install <- FALSE
first_run <- TRUE

if (install == TRUE){
  
  install.packages(c('ff','ffbase','readstata13','rms','survival','xlsx','randomForest','glmnet','ROCR','ggplot2'))
  install.packages('../CALIBERdatamanage_0.1-14.tar.gz',repos=NULL) 
  
}


# Load required packages and helper functions
library(readstata13) # allows reading in of stata dta files
library(rms) # scoring functions and more
library(survival) # survival 
library(ROCR) # ROC curves and prediction metrics
library(xlsx) # read from excel files
library(CALIBERdatamanage) # allows to store flat file df's in hard drive
library(randomForest)
library(glmnet)
library(ggplot2)
source('kiddle_misc.R') # load helper functions


if (first_run == TRUE){
  
  # identify eligible patients and their demographics
  source('eligible_patients.R')
  save(copd_eligible,eligible,linked,all_copd,file='../output/proc_data/loaded.RData')
  
  # load flat files into flat file data frames (living in hard drive, NOT RAM)
  source('copd_ffdf.R')
  
  source('load_variables.R') 
  save(copd_eligible,file='../output/proc_data/variables.RData')
  
  source('survival_format.R') 
  save(surv_copd,i,file="../output/proc_data/surv_copd.RData")
  
  source('descriptive.R')
  
  source('add_blood.R') # add in platelet, albumin and crp
  save(surv_copd,i,file="../output/proc_data/surv_copd.RData")
  

  
  source('cross_validation.R')
  save(copd_train,copd_test,r5cv10_train,r5cv10_test,file='../output/proc_data/CV.Rdata')
  
  source('run_cv.R')
  
  
  # now group 2
  source('kiddle_misc_nolink.R') # load helper functions
  source('eligible_patients_nolink.R')
  source('copd_ffdf_nolink.R')
  
  source('load_variables_nolink.R') 
  save(copd_eligible,file='../output/proc_data/variables_nolink.RData')
  
  source('descriptive_nolink.R')
  
  source('cross_validation_unlink.R')
  
  
  
  
  
  
} else {
  
  # identify eligible patients and their demographics
  load('../output/proc_data/loaded.RData')
  
  # load flat file data frames (living in hard drive, NOT RAM)
  load.ffdf('../output/proc_data/additional_eligible') 
  load.ffdf('../output/proc_data/clinical_full_before') 
  load.ffdf('../output/proc_data/immunisation_before') 
  load.ffdf('../output/proc_data/referral_before') 
  load.ffdf('../output/proc_data/test_before') 
  load.ffdf('../output/proc_data/therapy_before') 
  
  # 
  load(file='../output/proc_data/variables.RData')
  
  load(file="../output/proc_data/surv_copd.RData")
  
  load(file='../output/proc_data/CV.Rdata')
  
  load.ffdf('../output/proc_data/clinical_before_nolink')
  load.ffdf('../output/proc_data/additional_eligible_nolink') 
  load.ffdf('../output/proc_data/referral_before_nolink') 
  load.ffdf('../output/proc_data/test_before_nolink') 
  load.ffdf('../output/proc_data/therapy_before_nolink') 
  
}






#source('survival_format.R') 
#save(surv_copd,i,file="../output/proc_data/surv_copd.RData")
load("../output/proc_data/surv_copd.RData")

source('descriptive.R')

#source('silvia.R')
#source('cross_validation.R')
#save(copd_train,copd_test,cv10_train,cv10_test,ncv_trains,ncv_tests,file='../output/proc_data/CV.Rdata')
load(file='../output/proc_data/CV.Rdata')
