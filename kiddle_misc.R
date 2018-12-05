dateNumeric <- function(dates,...){ 
  # function to turn dates into numbers, not sure if I still use
  
  as.numeric(as.Date(dates,...))
  
}

find_medcodes <- function(medcodes,lm_time=0){
  # given a list of medcodes, finds them in all tables, ATM focused on COPD diagnosis,
  # later will modify to work with landmark times (lm_time)
  
  # empty fdff to save results in
  relevant <- as.ffdf(data.frame(patid=NA,eventdate=as.Date('1990-01-01'),medcode=NA))
  
  # if no landmark time specified, extract before COPD diagnosis
  if(lm_time==0){
    
    print('from COPD diagnosis')
    
    # for each table, find instances, store them in relevant (ffdf) 
    # and deal with errors if nothing is found
    
    print('Clinical')
    
    a <- tryCatch(
      {
        res <- clinical_before[clinical_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
        relevant <- ffdfappend(relevant,res)
      }, warning = function(w) {}, error = function(e) {}, finally = {}
    )
    
    print('Referral')
    
    a <- tryCatch(
      {
        res <- referral_before[referral_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
        relevant <- ffdfappend(relevant,res)
      }, warning = function(w) {}, error = function(e) {}, finally = {}
    )
    
    print('Test')
    
    a <- tryCatch(
      {
        res <- test_before[test_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
        relevant <- ffdfappend(relevant,res)
      }, warning = function(w) {}, error = function(e) {}, finally = {}
    )
    
    print('Immunisation')
    
    a <- tryCatch(
      {
        res <- immunisation_before[immunisation_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
        relevant <- ffdfappend(relevant,res)
      }, warning = function(w) {}, error = function(e) {}, finally = {}
    )
    
    
  } else {
    
    if(lm_time==2012){
      
      print('2012')
      
      # for each table, find instances, store them in relevant (ffdf) 
      # and deal with errors if nothing is found
      
      print('Clinical')
      
      a <- tryCatch(
        {
          res <- clinical2012_before[clinical2012_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
          relevant <- ffdfappend(relevant,res)
        }, warning = function(w) {}, error = function(e) {}, finally = {}
      )
      
      print('Referral')
      
      a <- tryCatch(
        {
          res <- referral2012_before[referral2012_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
          relevant <- ffdfappend(relevant,res)
        }, warning = function(w) {}, error = function(e) {}, finally = {}
      )
      
      print('Test')
      
      a <- tryCatch(
        {
          res <- test2012_before[test2012_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
          relevant <- ffdfappend(relevant,res)
        }, warning = function(w) {}, error = function(e) {}, finally = {}
      )
      
      print('Immunisation')
      
      a <- tryCatch(
        {
          res <- immunisation2012_before[immunisation2012_before$medcode %in% medcodes$medcode,c('patid','eventdate','medcode')]
          relevant <- ffdfappend(relevant,res)
        }, warning = function(w) {}, error = function(e) {}, finally = {}
      )
      
      
    }

    
  }
  
  return(relevant)
  
}

silviaScoreM <- function(data){
# function to calculate Cambridge Multimorbidity Score for mortality
  
  score <- numeric(dim(data)[1])
  
  score <- score + (data$hyp == 'present') * -2.09
  score <- score + ((data$anx == 'present') | (data$dep == 'present')) * 7.04
  score <- score + (data$pnc == 'present') * 16.46
  score <- score + (data$hel == 'present') * -3.94
  score <- score + (data$ibs == 'present') * -1.33
  score <- score + (data$ast == 'present') * -2.73
  score <- score + (data$dia == 'present') * 10.23
  score <- score + (data$chd == 'present') * 4.22
  score <- score + (data$ckd == 'present') * 16.61
  score <- score + (data$atr == 'present') * 22.14
  score <- score + (data$con == 'present') * 35.42
  score <- score + (data$str == 'present') * 20.63
  score <- score + (data$copd == 'present') * 42.50
  score <- score + (data$rhe == 'present') * -0.39
  score <- score + (data$can == 'present') * 62
  score <- score + (data$ap == 'present') * 12.72
  score <- score + (data$hf == 'present') * 43.47
  score <- score + (data$dem == 'present') * 124.42
  score <- score + (data$scz == 'present') * 7.2
  score <- score + (data$epi == 'present') * 18.26
  
}

silviaScoreG <- function(data){
  # function to calculate general Cambridge Multimorbidity Score
  
  score <- numeric(dim(data)[1])
  
  score <- score + (data$hyp == 'present') * 0.08
  score <- score + ((data$anx == 'present') | (data$dep == 'present')) * 0.5
  score <- score + (data$pnc == 'present') * 0.92
  score <- score + (data$hel == 'present') * 0.09
  score <- score + (data$ibs == 'present') * 0.21
  score <- score + (data$ast == 'present') * 0.19
  score <- score + (data$dia == 'present') * 0.75
  score <- score + (data$chd == 'present') * 0.49
  score <- score + (data$ckd == 'present') * 0.53
  score <- score + (data$atr == 'present') * 1.34
  score <- score + (data$con == 'present') * 1.12
  score <- score + (data$str == 'present') * 0.80
  score <- score + (data$copd == 'present') * 1.46
  score <- score + (data$rhe == 'present') * 0.43
  score <- score + (data$can == 'present') * 1.53
  score <- score + (data$ap == 'present') * 0.65
  score <- score + (data$hf == 'present') * 1.18
  score <- score + (data$dem == 'present') * 2.50
  score <- score + (data$scz == 'present') * 0.64
  score <- score + (data$epi == 'present') * 0.92
  
}

# is a patient known to be dead by X years post COPD diagnosis?
dead_by <- function(years,censored,time){censored == 0 & years < time}

# how many patients in this subset are dead X years post COPD diagnosis?
num_dead <- function(subset,time){length(which(dead_by(surv_copd$years[subset],surv_copd$censored[subset],time)))}

# add up the number of true instances
sumBool <- function(a){sum(as.numeric(a))}

# count number of NA (missing) data
numNA <- function(a){length(which(is.na(a)))}


# given cox model and new data, predict probability of five years survival
pred_coxph <- function(fit,data){
  
  #fit <- coxph(SurvObj2 ~ age + gender + smoke,data=small)
  
  # multiply coefficients by patients variable values to get prognostic scores
  X_new <- model.matrix(fit$formula,data= data)
  X_b <- fit$coef%*%t(X_new[,-1])
  
  # extract baseline hazard and index for five year prediction
  bh <- basehaz(fit)
  five_years <- dim(bh)[1]
  
  # use cox model to extract probability
  exp(-bh[five_years,1])^exp(X_b)
  
  #plot(bh[,2],exp(-bh[,1]))
  #lines(bh[,2],exp(-bh[,1])^exp(X_b[1]),col='blue')
}

R2BGLiMS_format <- function(data,old=vector('list',1)){
  # I was going to use R2BGLiMS, but didn't in the end. But the format is useful for
  # other things. Essentially standardise all data.
  
  # saves me renaming some old code variables
  small <- data
  
  # prepare output
  out <- vector('list',4)
  names(out) <- c('data','scaling_factors','means','year5_scaled')
  
  # columns that should be converted to integers
  cols_int <- c(12:19,21:58)[-c(2,4,5,6,7,45)]
  
  # convert them to integers
  for (i in cols_int){
    
    # missingness indicators are T/F, not factor
    if(i != 14 & i != 19 & i != 58){
      
      small[,i] <- as.integer(small[,i])-1
      
    } else {
      
      small[,i] <- as.integer(small[,i])
      
    }
    
  }
  
  # deal with factors, creating contrasts
  formula <- as.formula(paste('SurvObj~',paste(colnames(small)[c(12:19,21:58)][-c(17,14,10,27,41,24,25)],collapse='+')))
  
  mat <- model.matrix(formula,data= small)[,-1]
  
  std_small <- cbind(small[,1:10],mat)
  
  
  # scale all variables by dividing by standard deviation, from this or a previous run
  for (i in 11:dim(std_small)[2]){
    
    if (length(old) == 4){
      
      out$scaling_factors[i-10] <- old$scaling_factors[i-10]
      
    } else {
      
      out$scaling_factors[i-10] <- sd(mat[,i-10])
      
    }
    
    std_small[,i] <- mat[,i-10]/out$scaling_factors[i-10]
    
  }
  
  
  
  #colnames(std_small)[17:47]
  
  
  # create pairwise interactions between co-morbidities
  for (i in 1:30){
    
    for (j in (i+1):31){
      
      if (i == 30 & j == 31){}
      else {
        
        str <- paste(colnames(std_small)[20:50][i],colnames(std_small)[20:50][j],sep='_')
        
        std_small[,str] <- std_small[,c(20:50)[i]] * std_small[,c(20:50)[j]]
        
      }
      
    }
    
  }
  
  
  # center variables
  for (i in c(11:514)){
    
    
    
    if (length(old) == 4){
      
      out$means[i-10] <- old$means[i-10]
      
    } else {
      
      out$means[i-10] <- mean(std_small[,i])
      
    }
    
    std_small[,i] <- std_small[,i] - out$means[i-10]
    
  }
  
  
  # only need for R2BGLiMs
  out$year5_scaled <- 5/max(std_small$years)
  std_small$years <- std_small$years/max(std_small$years)# Recommend scaling survival times to between 0 and 1
  
  std_small$dead <- 1-std_small$censored
  
  out$data <- std_small
  
  return(out)
  
}

# function to add co-morbidity data for codelists with ever recorded inclusion criteria
ever_recorded <- function(codes,...){
  
  # for each eligible subject, see if code is present
  relevant <- as.data.table.ffdf(find_medcodes(codes,...))
  
  # keep earliest event per patient
  setorderv(relevant,c('patid','eventdate'))
  relevant <- relevant[!duplicated(relevant[,'patid']),]
  
  relevant <- relevant[-1,]
  
  #extract first event per patient
  return(relevant)
  
}

find_latest <- function(relevant,ffdf=T){

  if (ffdf){relevant <- as.data.table.ffdf(relevant)}
  
  # keep latest event per patient
  setorder(relevant,patid,-eventdate)
  relevant <- relevant[!duplicated(relevant[,'patid']),]
 
  return(relevant) 
}

# legacy functions, not needed
surv_weibull <- function(intercept,log_scale,t){exp(-(exp((log(t)-intercept)/exp(log_scale) )))}

pred_weibull <- function(scale,coef,formula,data,t){
  
  mat <- model.matrix(formula,data= data)
  
  surv_weibull(coef%*%t(mat),log(scale),t)
  
}

pred_coxnet <- function(b,data){
  
  X_new <- model.matrix(fit$formula,data= data)
  
  X_b <- b%*%t(data[,names(b)])
  
  bh <- basehaz(fit)
  
  five_years <- dim(bh)[1]
  
  exp(-bh[five_years,1])^exp(X_b)
  
  #plot(bh[,2],exp(-bh[,1]))
  #lines(bh[,2],exp(-bh[,1])^exp(X_b[1]),col='blue')
}