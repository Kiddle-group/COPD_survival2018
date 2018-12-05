print('Survival format data - time from first COPD diagnosis')


print('..load death data from ONS')
#add in ONS death data
death <- read.table(file="../data/Final + new linkage/linkage_download_jan2018/death_patient_16_276R.txt",sep="\t",header=T) 

dod <- as.Date(death$dod,format="%d/%m/%Y") #date of death
censor_date <- as.Date(max(dod,na.rm=T)) # last recorded death date, censor date

# last death date 19th Sept 2017, treat as censor date
# note, laster to do stopped cox regression i censor at 5-years survival

# create survival format data
print('..create time to event format data')
patid_to_use <- copd_eligible$patid
surv_copd <- data.frame(matrix(NA,length(patid_to_use),3))
colnames(surv_copd) <- c('patid','years','censored')
surv_copd[,1] <- patid_to_use

# loop over all patients, recording info
for (i in 1:dim(surv_copd)[1]){
  
  tmp <- data.frame(patid=patid_to_use[i],years=NA,censored=NA)
  
  # find patient in death data and patient dataset
  ind <- which(death$patid == patid_to_use[i])
  ind2 <- which(copd_eligible$patid == patid_to_use[i])
  
  # if patient can be found in both
  if (length(ind)>0 & !all(is.na(death[ind,"dod"])) ){
    
    # calculate years of survival
    tmp[1,'years'] <- as.numeric(as.Date(dod[ind][which(!is.na(dod[ind]))[1]]) - as.Date(copd_eligible[ind2,"first_copd"])) / 365.25
    
    tmp[1,'censored'] <- 0
    
  } else {
    # they are not yet known to have died, mark them as censored
    # note,  later censored at 5-years survival 
    
    tmp[1,'years'] <- as.numeric(censor_date - as.Date(copd_eligible[ind2,"first_copd"])) / 365.25
    
    tmp[1,'censored'] <- 1
    
  }
  
  #add patients data to survival format data frame
  surv_copd[i,] <- tmp
  
  # save and report back every 10,000 patients
  if (i%%10000 == 0){
    
    print(paste(i,'out of',dim(surv_copd)[1]))
    
    save(surv_copd,i,file="../output/proc_data/surv_copd.RData")
    
  }
  
}

# save data
save(surv_copd,i,file="../output/proc_data/surv_copd.RData")

print('..merge with other data')
# merge survival and other data
surv_copd <- merge(surv_copd,copd_eligible,by="patid")

print('..add raw ONS data')
little_death <- death[,c("patid","dod")]
colnames(little_death) <- c("patid","ons_dod")

# add ONS death data
surv_copd <- merge(surv_copd,little_death,by="patid",all.x=TRUE)

save(surv_copd,i,file="../output/proc_data/surv_copd.RData")

print('..final eligibility selection')
#remove deaths before COPD diagnosis and COPD diagnoses before age 35, other exclusion criteria applied in other file
surv_copd <- surv_copd[which(surv_copd$years>0),]
surv_copd <- surv_copd[which(surv_copd$year_copd >= surv_copd$year35),]

# stopped cox data, censor at time horizon of interest to give model fair chance
surv_copd$years5 <- surv_copd$years
surv_copd$censored5 <- surv_copd$censored
surv_copd$censored5[surv_copd$years5 > 5] <- 1
surv_copd$years5[surv_copd$years5 > 5] <- 5


rm(death,little_death,censor_date,dod,patid_to_use,tmp)
