

#AECOPD only happen after COPD diagnosis, so don't worry about for now
#aecopd <- read.dta13('../data/other_files/AECOPD/all_aecopd2.dta')
#aecopd_eligible <- aecopd[aecopd$patid %in% copd_eligible$patid,]

#MRC, dyspnea occurs rarely before diagnosis
#
codes <- read.xlsx('../data/MRC.xlsx',1) # load codes
first <- ever_recorded(codes) # find first instance per patient

# number of eligible patients with dyspnea data
dim(first)

#FEV1 in last 12 months, remove missing FEV1
fev1 <- read.dta13('../data/Steve_fev_pred1.dta') # data from Quint group
fev1_nona <- fev1[!is.na(fev1$percent_fev),]
#fev1_nona$fev1[fev1_nona$data3 == 89] <- fev1_nona$fev1[fev1_nona$data3 == 89]*1000


# extract just values from eligible patients within year before COPD diagnosis
fev1_eligible <- fev1_nona[fev1_nona$patid %in% copd_eligible$patid,]
fev1_eligible$post_diagnosis <- dateNumeric(fev1_eligible$eventdate) - dateNumeric(fev1_eligible$first_copd)
fev1_before_nolink <- fev1_eligible[(fev1_eligible$post_diagnosis <= 0) & (fev1_eligible$post_diagnosis >= -365),]

# reorder by patient, days before diagnosis and highest reading for the day
fev1_before_nolink <- fev1_before_nolink[order(fev1_before_nolink$patid,fev1_before_nolink$post_diagnosis,fev1_before_nolink$percent_fev,decreasing = c(F,T,T)),]

# nice trick to give per patient order
fev1_before_nolink$order <- ave(rep(1,nrow(fev1_before_nolink)),fev1_before_nolink[,1],FUN=seq_along )

# extract latest observation
fev1_before_nolink <- fev1_before_nolink[fev1_before_nolink$order==1,]

# extract FEV1 percent predicted (i.e. airflow obstruction)
fev1_before_nolink$airflow <- fev1_before_nolink$percent_fev

# merge with patient data
copd_eligible <- merge(copd_eligible,fev1_before_nolink[,c('patid','airflow')],by='patid',all.x=T)

# set up not recorded indicator
missing <- is.na(copd_eligible$airflow)
copd_eligible$airflow[missing] <- 0
copd_eligible$airflowMiss <- missing

rm(fev1,fev1_eligible,fev1_nona,fev1_before_nolink,missing) # remove tmp files

# age at fist COPD
copd_eligible$age_copd <- (dateNumeric(copd_eligible[,'first_copd']) - dateNumeric(copd_eligible[,'dob']))/365




#smoking, make an exception, use memory instead of hard drive to avoid strange error

# find instances of medcodes
codes <- read.dta13('../data/smoking_codes.dta')[,c(1,3)]
relevant <- as.data.table.ffdf(find_medcodes(codes))

# recode to make consistent with additional data
codes[codes[,"smokstatus"] == 1,"smokstatus"] <- 3
codes[codes[,"smokstatus"] == 2,"smokstatus"] <- 1

# add data to relevant
relevant <- merge(relevant,all_copd,by='patid')
relevant <- merge(relevant,codes,by='medcode')
relevant$post_diagnosis <- as.numeric(relevant$eventdate - relevant$first_copd)

# find instances use enttype and linked additional data

clinical_relevant2 <- as.data.table.ffdf(clinical_before_nolink[clinical_before_nolink$enttype[] == 4,c('patid','adid','eventdate','first_copd')])
additional_relevant <- as.data.table.ffdf(additional_eligible_nolink[additional_eligible_nolink$enttype[]==4,c('patid','adid','data1')])
colnames(additional_relevant)[3] <- 'smokstatus'

additional_relevant$smokstatus <- as.numeric(as.character(additional_relevant$smokstatus))

# use clinical enttype to identify presence of additional data, then track it down
clinical_relevant2 <- merge(clinical_relevant2,additional_relevant,by=c('patid','adid'))

clinical_relevant2$post_diagnosis <- as.numeric(clinical_relevant2$eventdate - clinical_relevant2$first_copd)

# merge data derived using medcodes and enttype
smoke <- rbind(relevant[,c('patid','smokstatus','post_diagnosis')],clinical_relevant2[,c('patid','smokstatus','post_diagnosis')])

rm(relevant,clinical_relevant2)

# neat trick to extract most recent data for each patient
setorderv(smoke,c('patid','post_diagnosis'),order=c(1,-1))
smoke$order <- ave(rep(1,nrow(smoke)),smoke[,1],FUN=seq_along )
smoke <- smoke[smoke$order==1,]

# add most recent smoking data to patient dataset
copd_eligible <- merge(copd_eligible,smoke[,c('patid','smokstatus')],by='patid',all.x=T)

rm(smoke, codes)

# BMI

# find medcode instances
clinical_relevant <- clinical_before_nolink[clinical_before_nolink$enttype[] == 13,c('patid','eventdate','adid','first_copd','enttype')]
clinical_relevant$before_diagnosis <- clinical_relevant$first_copd - clinical_relevant$eventdate

# sort by patid, then how many days before diagnosis, efficiently using ffdf
clinical_relevant <- as.data.table.ffdf(clinical_relevant)

# neat trick to extract most recent data for each patient
setorderv(clinical_relevant,c('patid','before_diagnosis'),order=c(1,1))

# keep only most recent results
clinical_relevant <- clinical_relevant[!duplicated(clinical_relevant[,'patid']),]

# but actual bmi is in linked additional data, so extract and merge
additional_relevant <- additional_eligible_nolink[additional_eligible_nolink$enttype[]==13,c('patid','adid','data3')]
bmi_both <- merge(clinical_relevant,additional_relevant,by=c('patid','adid'))
bmi_both$bmi <- as.numeric(as.character(bmi_both$data3))
copd_eligible <- merge(copd_eligible,bmi_both[,c('patid','bmi')],by='patid',all.x=T)

# BMI > 70 possible but very unlikely (max val 325000!!)
copd_eligible$bmi[which(copd_eligible$bmi>70)] <- NA

missing <- is.na(copd_eligible$bmi)

# set up not recorded indicator
copd_eligible$bmi[missing] <- 0
copd_eligible$bmiMiss <- missing

rm(clinical_relevant,additional_relevant,bmi_both,missing)



# need to identify co-morbidities 
#
# using codes and algorithms from http://www.phpc.cam.ac.uk/pcu/cprd_cam/codelists/ (plus soon to be released update!)
#

print('identify comorbid')
      
      
# CKD
print('..extract eGFR for CKD')

# find egfr data
egfr <- test_before_nolink[test_before_nolink$enttype[]==466,c('patid','eventdate','data1','data2')] 

#load('../output/proc_data/egfr.RData')

#data1 - Operator (OPR) 0, 1 < , 2 <=,3 =, 4 >, 5 >=, 6 ~
#data2 - eGFR Value
#data3 - Unit of measure (SUM) - 90 is mL/min
#data4 - Qualifier	(TQU) NA
#data5 - Normal range from NA	
#data6 - Normal range to NA
#data7 - Normal range basis (POP) NA

egfr <- egfr[!duplicated(egfr),] # remove duplicates

egfr <- egfr[!is.na(egfr$data2),] # remove missing values

# remove rare confusing cases
ind <- !(egfr$data1==1 & egfr$data2==60) 
egfr <- egfr[ind,] 

# remove unclassifiable cases
ind <- !(egfr$data1==1 & egfr$data2==90) 
egfr <- egfr[ind,] 
ind <- !(egfr$data1==2 & (egfr$data2==90 | egfr$data2==60 )) 
egfr <- egfr[ind,] 

# remove missing data
ind <- !(egfr$data2==0 ) 
egfr <- egfr[ind,] 

# remove unrealistic range data
ind <- !(egfr$data2>200 ) 
egfr <- egfr[ind,] 

# add some info
egfr <- merge(egfr,as.ffdf(all_copd),by='patid',all.x=T)
egfr$before_diagnosis <- egfr$first_copd - egfr$eventdate

# sort by patid, then how many days before diagnosis, efficiently using ffdf

egfr <- as.data.table.ffdf(egfr)
setorderv(egfr,c('patid','before_diagnosis'),order=c(1,1))

# convert to data table only when necessary
# here i can't work out how to easily extract
# first two measures for each person from ffdf

egfr[,'test'] <- ave(rep(1,nrow(egfr)),egfr[,1],FUN=seq_along )
egfr[,'num_test'] <- ave(rep(1,nrow(egfr)),egfr[,1],FUN=sum )

egfr <- egfr[egfr$num_test > 1,]
egfr <- egfr[egfr$test < 3,]

# neat trick to get max of all values for individual
egfr[,'max'] <- ave(egfr$data2,egfr[,1],FUN=max )
egfr <- egfr[egfr$test == 1,]

# extract egfr, and whether above or below relevant cutpoints
egfr$egfr <- egfr$max
egfr$c_ckd_copd <- egfr$max < 60 # used in CMS
egfr$cci_ckd_copd <- egfr$max < 30 # used in CCI

# add egfr/ckd data to patient dataset
copd_eligible <- merge(copd_eligible,egfr[,c('patid','egfr','c_ckd_copd','cci_ckd_copd')],by='patid',all.x=T)
copd_eligible$c_ckd_copd[is.na(copd_eligible$c_ckd_copd)] <- F
copd_eligible$cci_ckd_copd[is.na(copd_eligible$cci_ckd_copd)] <- F

# set up egfr not recorded at least twice indicator
missing <- is.na(copd_eligible$egfr)
copd_eligible$egfr[missing] <- 0
copd_eligible$egfrMiss <- missing

rm(egfr)


# GORD in previous year

# find instances
codes <- read.dta13('../data/gord.dta')
relevant <- as.data.table.ffdf(find_medcodes(codes))

# add info
relevant <- merge(relevant,all_copd,by='patid')
relevant$before_diagnosis <- as.numeric(relevant$first_copd - relevant$eventdate)

# require mention in last year
relevant <- relevant[relevant$before_diagnosis < 365,]

# keep only most recent results
setorderv(relevant,c('patid','before_diagnosis'))
relevant <- relevant[!duplicated(relevant[,'patid']),]
relevant$gord <- T

# add gord data to patient dataset
copd_eligible <- merge(copd_eligible,relevant[,c('patid','gord')],by='patid',all.x=T)
copd_eligible$gord[is.na(copd_eligible$gord)] <- F

# used later to make code robust to added variables
colnum <- dim(copd_eligible)[2]

# in case of crashes, save progress
save(copd_eligible,file='../output/proc_data/1_3.RData')

# for each simple defintion, find first instances and add data to patient dataset
print('..simple definitions')

# alcohol problems

codes <- read.csv('../../Codelists/CPRDCAM_ALC095_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_ap_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# anorexia bulimia

codes <- read.csv('../../Codelists/CPRDCAM_ANO049_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_ab_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# atrial fibrillation

codes <- read.csv('../../Codelists/CPRDCAM_ATR033_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_atr_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# blindness and low vision

codes <- read.csv('../../Codelists/CPRDCAM_BLI098_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_bli_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# bronchiectasis

codes <- read.csv('../../Codelists/CPRDCAM_BRO050_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_bro_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# cancer v1 # not needed for v1.1

#codes <- read.csv('../../Codelists/CPRDCAM_CAN023_MEDCODES.csv') # load codes
#first <- ever_recorded(codes) # find first instance per patient
#colnames(first)[22] <- "c_can_date" # data of first instance per patient
#copd_eligible <- merge(copd_eligible,first[,c(1,22)],by='patid',all.x = TRUE)
#rm(codes,first) #remove tmp data

# chronic liver disease

codes <- read.csv('../../Codelists/CPRDCAM_CLD081_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_cld_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# chronic sinusitis

codes <- read.csv('../../Codelists/CPRDCAM_SIN099_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_sin_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# coronary heart disease 

codes <- read.csv('../../Codelists/CPRDCAM_CHD086_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_chd_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# dementia

codes <- read.csv('../../Codelists/CPRDCAM_DEM083_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_dem_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# diabetes

codes <- read.csv('../../Codelists/CPRDCAM_DIB062_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_dia_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# diverticular disease of intestine

codes <- read.csv('../../Codelists/CPRDCAM_DIV089_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_div_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# hearing loss

codes <- read.csv('../../Codelists/CPRDCAM_HEL076_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_hel_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# Heart failure

codes <- read.csv('../../Codelists/CPRDCAM_HEF035_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_hf_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# Hypertension

codes <- read.csv('../../Codelists/CPRDCAM_HYP002_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_hyp_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# IBD

codes <- read.csv('../../Codelists/CPRDCAM_IBD097_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_ibd_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# Learning disability

codes <- read.csv('../../Codelists/CPRDCAM_LEA102_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_lea_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# Multiple sclerosis

codes <- read.csv('../../Codelists/CPRDCAM_MSC052_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_ms_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# Peripheral vascular disease

codes <- read.csv('../../Codelists/CPRDCAM_PVD074_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_pvd_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# Parkinson's disease

codes <- read.csv('../../Codelists/CPRDCAM_PRK051_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_prk_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# Prostate disorders

codes <- read.csv('../../Codelists/CPRDCAM_PRO092_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_pro_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# Psychoactive substance misuse

codes <- read.csv('../../Codelists/CPRDCAM_OPS096_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_ops_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# Rheumatoid arthritis

codes <- read.csv('../../Codelists/CPRDCAM_RHE088_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_rhe_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# Stroke

codes <- read.csv('../../Codelists/CPRDCAM_STR029_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_str_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

# Thyroid disorders

codes <- read.csv('../../Codelists/CPRDCAM_THY016_MEDCODES.csv') # load codes
first <- ever_recorded(codes) # find first instance per patient
colnames(first)[2] <- "c_thy_date" # date of first instance per patient
copd_eligible <- merge(copd_eligible,first[,c(1,2)],by='patid',all.x = TRUE)
rm(codes,first) #remove tmp data

print('..simple, by COPD diagnosis')

# legacy code

# for ever recorded co-morbidities, turn date into above
#tmp <- vapply(copd_eligible[,33:56],dateNumeric,numeric(dim(copd_eligible)[1]))

# dates after first COPD date?
#tmp2 <- tmp < rep(dateNumeric(copd_eligible[,'first_copd']),24)

# set missing to FALSE as well, means never mention of readcode
#tmp2[which(is.na(tmp2),arr.ind=T)] <- FALSE

# extract co-morbidity labels from column names of tmp2
#tmp_mat <- matrix(unlist(strsplit(colnames(tmp2),'_')),3,24)

# new labels to show indicator of co-morbidity by first COPD date
#colnames(tmp2) <- paste('c',tmp_mat[2,],sep='_')

# identify comorbidities with a missing date, i.e. that haven't been diagnosed
tmp <- vapply(copd_eligible[,32:55],is.na,numeric(dim(copd_eligible)[1]))

# use comorbidity date column names to make comorbidity present column names
tmp_mat <- matrix(unlist(strsplit(colnames(tmp),'_')),3,24)
colnames(tmp) <- paste('c',tmp_mat[2,],sep='_')

# add to dataset if not missing, i.e. if a diagnosis has been made
copd_eligible <- cbind(copd_eligible,!tmp)

rm(tmp,tmp_mat)



print('..complex definitions')
# the rest of the comorbidities are not 'ever recorded', and so require more work to identify
# for each comorbidity we give the definition at the beggining

# Anxiety (read code in last 12 months OR >= 4 anxiolytic/hypnotic presciption in last 12 months)

# load med/prod codes
codes <- read.csv('../../Codelists/CPRDCAM_ANX021_MEDCODES.csv') 
prods <- read.csv('../../Codelists/CPRDCAM_ANX101_PRODCODES.csv')

# find instances of medcodes
relevant <- as.data.table.ffdf(find_medcodes(codes))
relevant <- merge(relevant,all_copd,by='patid')

# num days between first COPD and read code recording
relevant$before_diagnosis <- as.numeric(relevant$first_copd - relevant$eventdate)
relevant <- relevant[relevant$before_diagnosis<=365,] # require to be within 12 months
relevant$c_anx <- T 

# find instances of prodcodes
therapy_relevant <- therapy_before_nolink[therapy_before_nolink$prodcode[] %in% prods$prodcode,c('patid','eventdate','first_copd')]
therapy_relevant$before_diagnosis <- therapy_relevant$first_copd - therapy_relevant$eventdate

 # require to be within 12 months
therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])

# how many per patient? need >=4 
therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
therapy_relevant$c_anx <- T

# combine as can be either, but leave unique patids and then merge
rels <- rbind(relevant[,c('patid','c_anx')],therapy_relevant[,c('patid','c_anx')])
rels <- rels[!duplicated(rels[,'patid']),]

copd_eligible <- merge(copd_eligible,rels,by='patid',all.x=T)
copd_eligible$c_anx[is.na(copd_eligible$c_anx)] <- F


# Asthma (read code in between 24 - 60 months AND presciption in last 12 months)
# Quint codes and rules, basically COPD and asthma meds overlap, and 0-2 yr likely misdiagnosed 
codes <- read.dta13('../data/asthma.dta')
#codes <- read.csv('../../Codelists/CPRDCAM_AST078_MEDCODES.csv') # CPRD @ cambridge codes, not used
prods <- read.csv('../../Codelists/CPRDCAM_AST055_PRODCODES.csv')

# keep only rows with right codes
relevant <- as.data.table.ffdf(find_medcodes(codes))
relevant <- merge(relevant,all_copd,by='patid')

# num days between first COPD and read code recording
relevant$before_diagnosis <- as.numeric(relevant$first_copd - relevant$eventdate)
relevant <- relevant[(relevant$before_diagnosis>=2*365)&(relevant$before_diagnosis<=5*365),]
relevant <- relevant[!duplicated(relevant[,'patid']),]

# treatment in last year
therapy_relevant <- therapy_before_nolink[therapy_before_nolink$prodcode[] %in% prods$prodcode,c('patid','eventdate','first_copd')]
therapy_relevant$before_diagnosis <- therapy_relevant$first_copd - therapy_relevant$eventdate
therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])
therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]


# combine treatment and diagnosis info and add to dataset
rels <- merge(relevant,therapy_relevant,by='patid')
rels$c_ast <- T

copd_eligible <- merge(copd_eligible,rels[,c('patid','c_ast')],by='patid',all.x=T)
copd_eligible$c_ast[is.na(copd_eligible$c_ast)] <- F


# cancer v1.1, any read code in last 5 years

# identify instances, add info
codes <- read.csv('../../Codelists/CPRDCAM_CAN023_MEDCODES.csv')
relevant <- as.data.table.ffdf(find_medcodes(codes))
relevant <- merge(relevant,all_copd,by='patid')

# set to TRUE if recorded in last 5 years
relevant$before_diagnosis <- as.numeric(relevant$first_copd - relevant$eventdate)
relevant <- relevant[(relevant$before_diagnosis<=5*365),]
relevant <- relevant[!duplicated(relevant[,'patid']),]
relevant$c_can <- T

copd_eligible <- merge(copd_eligible,relevant[,c('patid','c_can')],by='patid',all.x=T)
copd_eligible$c_can[is.na(copd_eligible$c_can)] <- F


# Constipation - 4 or more laxative prescriptions in last year
prods <- read.csv('../../Codelists/CPRDCAM_CON084_PRODCODES.csv')

# identify instances within last year, count and require at least 4
therapy_relevant <- therapy_before_nolink[therapy_before_nolink$prodcode[] %in% prods$prodcode,c('patid','eventdate','first_copd')]
therapy_relevant$before_diagnosis <- therapy_relevant$first_copd - therapy_relevant$eventdate
therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])

therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
therapy_relevant$c_con <- T

copd_eligible <- merge(copd_eligible,therapy_relevant[,c('patid','c_con')],by='patid',all.x=T)
copd_eligible$c_con[is.na(copd_eligible$c_con)] <- F


# Depression (read code in last 12 months OR >= 4 presciption in last 12 months)

codes <- read.csv('../../Codelists/CPRDCAM_DEP003_MEDCODES.csv') # load codes
prods <- read.csv('../../Codelists/CPRDCAM_DEP100_PRODCODES.csv')

# keep only rows with right codes
relevant <- as.data.table.ffdf(find_medcodes(codes))
relevant <- merge(relevant,all_copd,by='patid')

# require to be within year of COPD diagnosis
relevant$before_diagnosis <- as.numeric(relevant$first_copd - relevant$eventdate)
relevant <- relevant[relevant$before_diagnosis<=365,]
relevant$c_dep <- T 

# identify treatment instances and require to be within a year
therapy_relevant <- therapy_before_nolink[therapy_before_nolink$prodcode[] %in% prods$prodcode,c('patid','eventdate','first_copd')]
therapy_relevant$before_diagnosis <- therapy_relevant$first_copd - therapy_relevant$eventdate
therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])

# count and require to be at least 4
therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
therapy_relevant$c_dep <- T

# combine as can be either, but leave unique patids
rels <- rbind(relevant[,c('patid','c_dep')],therapy_relevant[,c('patid','c_dep')])
rels <- rels[!duplicated(rels[,'patid']),]

copd_eligible <- merge(copd_eligible,rels,by='patid',all.x=T)
copd_eligible$c_dep[is.na(copd_eligible$c_dep)] <- F


# Epilepsy (read code ever recorded AND presciption in last 12 months)
# and
# painful condition (4 or more prescriptions in last 12 months OR (4 or more anti-epleptics in absence of elilepsy ever recorded) 

codes <- read.csv('../../Codelists/CPRDCAM_EPI069_MEDCODES.csv') # load codes
epi <- read.csv('../../Codelists/CPRDCAM_PNC079_PRODCODES.csv')

# keep only rows with right codes
relevant <- as.data.table.ffdf(find_medcodes(codes))
relevant <- relevant[!duplicated(relevant[,'patid']),]

# identify treatment instances in year before COPD diagnosis
therapy_relevant <- therapy_before_nolink[therapy_before_nolink$prodcode[] %in% epi$prodcode,c('patid','eventdate','first_copd')]
therapy_relevant$before_diagnosis <- therapy_relevant$first_copd - therapy_relevant$eventdate
therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])

# any?
therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]

# merge, i.e. if not in both will have no entry
both <- merge(relevant,therapy_relevant,by='patid')
both$c_epi <- T
copd_eligible <- merge(copd_eligible,both[,c('patid','c_epi')],by='patid',all.x=T)
copd_eligible$c_epi[is.na(copd_eligible$c_epi)] <- F


# PAINFUL CONDITION
prods <- read.csv('../../Codelists/CPRDCAM_PNC004_PRODCODES.csv')

# identify treatment in the preceding year using new list
therapy_relevant2 <- therapy_before_nolink[therapy_before_nolink$prodcode[] %in% prods$prodcode,c('patid','eventdate','first_copd')]
therapy_relevant2$before_diagnosis <- therapy_relevant2$first_copd - therapy_relevant2$eventdate

therapy_relevant2 <- as.data.table.ffdf(therapy_relevant2[therapy_relevant2$before_diagnosis <= 365,])

# any? either in this list or previous
therapy_relevant2$num <- ave(rep(1,nrow(therapy_relevant2)),therapy_relevant2[,1],FUN=sum )
therapy_relevant2 <- therapy_relevant2[!duplicated(therapy_relevant2[,'patid']),]
therapy_relevant2 <- therapy_relevant2[therapy_relevant2$num>=4,]
therapy_relevant2$c_pnc <- T

therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
therapy_relevant <- therapy_relevant[!(therapy_relevant$patid %in% relevant$patid)]
therapy_relevant$c_pnc <- T

# combine as can be either, but leave unique patids
rels <- rbind(therapy_relevant[,c('patid','c_pnc')],therapy_relevant2[,c('patid','c_pnc')])
rels <- rels[!duplicated(rels[,'patid']),]
copd_eligible <- merge(copd_eligible,rels[,c('patid','c_pnc')],by='patid',all.x=T)
copd_eligible$c_pnc[is.na(copd_eligible$c_pnc)] <- F


# IBS (read code ever recorded OR 4 or more presciption in last 12 months)

codes <- read.csv('../../Codelists/CPRDCAM_IBS022_MEDCODES.csv') # load codes
prods <- read.csv('../../Codelists/CPRDCAM_IBS085_PRODCODES.csv')

# keep only rows with right codes
relevant <- as.data.table.ffdf(find_medcodes(codes))
relevant$c_ibs <- T 

# identify treatment in preceding year
therapy_relevant <- therapy_before_nolink[therapy_before_nolink$prodcode[] %in% prods$prodcode,c('patid','eventdate','first_copd')]
therapy_relevant$before_diagnosis <- therapy_relevant$first_copd - therapy_relevant$eventdate
therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])

# four or more?
therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
therapy_relevant$c_ibs <- T

# combine as can be either, but leave unique patids
rels <- rbind(relevant[,c('patid','c_ibs')],therapy_relevant[,c('patid','c_ibs')])
rels <- rels[!duplicated(rels[,'patid']),]

copd_eligible <- merge(copd_eligible,rels,by='patid',all.x=T)
copd_eligible$c_ibs[is.na(copd_eligible$c_ibs)] <- F


# Migrane - 4 or more migrane prescriptions in last year

prods <- read.csv('../../Codelists/CPRDCAM_MIG045_PRODCODES.csv')

# identify tratment instances
therapy_relevant <- therapy_before_nolink[therapy_before_nolink$prodcode[] %in% prods$prodcode,c('patid','eventdate','first_copd')]
therapy_relevant$before_diagnosis <- therapy_relevant$first_copd - therapy_relevant$eventdate

therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])

# four or more?
therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]
therapy_relevant$c_mig <- T

copd_eligible <- merge(copd_eligible,therapy_relevant[,c('patid','c_mig')],by='patid',all.x=T)
copd_eligible$c_mig[is.na(copd_eligible$c_mig)] <- F


# psoriasis or eczema (read code ever recorded AND 4 or more prescriptions in last 12 months) 

codes <- read.csv('../../Codelists/CPRDCAM_PSO090_MEDCODES.csv') # load codes
prods <- read.csv('../../Codelists/CPRDCAM_PSO091_PRODCODES.csv')

# keep only rows with right codes
relevant <- as.data.table.ffdf(find_medcodes(codes))
relevant <- relevant[!duplicated(relevant[,'patid']),]

# identify treatment instances in preceding year
therapy_relevant <- therapy_before_nolink[therapy_before_nolink$prodcode[] %in% prods$prodcode,c('patid','eventdate','first_copd')]
therapy_relevant$before_diagnosis <- therapy_relevant$first_copd - therapy_relevant$eventdate
therapy_relevant <- as.data.table.ffdf(therapy_relevant[therapy_relevant$before_diagnosis <= 365,])

# >= 4 prescriptions
therapy_relevant$num <- ave(rep(1,nrow(therapy_relevant)),therapy_relevant[,1],FUN=sum )
therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
therapy_relevant <- therapy_relevant[therapy_relevant$num>=4,]

# merge diagnoses and treatments, so only patients with both remain
both <- merge(relevant,therapy_relevant,by='patid')
both$c_pso <- T
copd_eligible <- merge(copd_eligible,both[,c('patid','c_pso')],by='patid',all.x=T)
copd_eligible$c_pso[is.na(copd_eligible$c_pso)] <- F


# scizophrenia (read code ever recorded OR lithium ever recorded) 

codes <- read.csv('../../Codelists/CPRDCAM_SCZ094_MEDCODES.csv') # load codes
prods <- read.csv('../../Codelists/CPRDCAM_SCZ108_PRODCODES.csv')

# medcodes
relevant <- as.data.table.ffdf(find_medcodes(codes))
relevant <- relevant[!duplicated(relevant[,'patid']),]
relevant$c_scz <- T

# prodcodes
therapy_relevant <- therapy_before_nolink[therapy_before_nolink$prodcode[] %in% prods$prodcode,c('patid','eventdate','first_copd')]
therapy_relevant <- therapy_relevant[!duplicated(therapy_relevant[,'patid']),]
therapy_relevant$c_scz <- T

# either
rels <- rbind(relevant[,c('patid','c_scz')],therapy_relevant[,c('patid','c_scz')])
rels <- rels[!duplicated(rels[,'patid']),]

copd_eligible <- merge(copd_eligible,rels,by='patid',all.x=T)
copd_eligible$c_scz[is.na(copd_eligible$c_scz)] <- F

rm(all_copd,both,codes,eligible,epi,ind,linked,missing,prods,relevant,rels,therapy_relevant,therapy_relevant2)


