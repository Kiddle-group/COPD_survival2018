# first need to convert files to csv in linux
# cat copd_cohort_march2017_Extract_Therapy_01.txt | tr "\\t" ","  > Therapy_01.csv

print('Making ffdfs, takes a long time but is worth it, only first run ')

copd_eligible <- data.table(copd_eligible)

# turn copd_eligible into ffdf, memory efficient object
eligible <- as.ffdf(copd_eligible)

print('CLINICAL')

# CLINICAL
#
# I have commented out an attempt to go from raw data as i can't correctly format the dates.
# Could be solved by splitting into multiple smaller files, but i went with a Quint group
# pre-existing file instead.

# have trouble reformatting the dates, so just go from dta file
files <- paste('../data/copd_cohort_march2017/Clinical_0',1:7,'.csv',sep='')

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
clinical <- importFFDF(filename = paste('../data/copd_cohort_march2017/Clinical_0',1:7,'.csv',sep=''),datecolnames = '',verbose = T)

save.ffdf(clinical,dir = '../output/proc_data/clinical')

clinical_eligible <- merge(clinical,eligible,by='patid')
save.ffdf(clinical_eligible,dir = '../output/proc_data/clinical_eligible')

# does this work?
clinical_eligible$eventdate <- as.character(clinical_eligible$eventdate) #works fine
clinical_eligible$eventdate <- as.Date(clinical_eligible$eventdate,format='%d/%m/%Y')

clinical_before <- clinical_eligible[clinical_eligible$eventdate < clinical_eligible$first_copd,]
save.ffdf(clinical_before,dir = '../output/proc_data/clinical_before')










#save.ffdf(clinical,dir = '../output/proc_data/clinical')

#copd_eligible$first_copd <- as.Date(copd_eligible$first_copd,format='%Y-%m-%d')



# load Quint group version of clinical with nice date formats
#clinical <- read.dta13('../data/Final + new linkage/Final_cohort/clinical_base.dta')
# clinical <- as.ffdf(clinical)

# extract only data for eligible and save
#clinical_eligible <- merge(clinical[,c('patid','eventdate','medcode','enttype','adid')],copd_eligible,by='patid')
#tmp <- data.table(clinical_eligible)
#tmp_ff <- as.ffdf(tmp)
#clinical_eligible <- tmp_ff
#save.ffdf(clinical_eligible,dir = '../output/proc_data/clinical_eligible')
#rm(tmp,tmp_ff)


# post-hoc turn dates into as.Date (does this work?)
#clinical_eligible$eventdate <- as.character(clinical_eligible$eventdate)
#clinical_eligible$eventdate <- as.Date(clinical_eligible$eventdate,format='%d/%m/%Y')

# extract just the entries that occur before COPD diagnosis and save
#clinical_before <- clinical_eligible[clinical_eligible$eventdate < clinical_eligible$first_copd,]
#save.ffdf(clinical_before,dir = '../output/proc_data/clinical_before')





print('TEST')

# Test
#

# read in raw files to ffdf objet
# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
test <- importFFDF(filename = paste('../data/copd_cohort_march2017/Test_0',1:7,'.csv',sep=''),datecolnames = '',verbose = T)
save.ffdf(test,dir = '../output/proc_data/test')

# extract only data for eligible and save
test_eligible <- merge(test,eligible,by='patid')
save.ffdf(test_eligible,dir = '../output/proc_data/test_eligible')

# post-hoc turn dates into as.Date

test_eligible$eventdate <- as.character(test_eligible$eventdate)
test_eligible$eventdate <- as.Date(test_eligible$eventdate,format='%d/%m/%Y')

# extract just the entries that occur before COPD diagnosis and save
test_before <- test_eligible[test_eligible$eventdate < test_eligible$first_copd,]
save.ffdf(test_before,dir = '../output/proc_data/test_before')

print('REFERRAL')

# Referral
#

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
referral <- importFFDF(filename = paste('../data/copd_cohort_march2017/Referral.csv',sep=''),datecolnames = '',verbose = T)
save.ffdf(referral,dir = '../output/proc_data/referral')

# extract only data for eligible and save
referral_eligible <- merge(referral,eligible,by='patid')
save.ffdf(referral_eligible,dir = '../output/proc_data/referral_eligible')

# post-hoc turn dates into as.Date

referral_eligible$eventdate <- as.character(referral_eligible$eventdate)
referral_eligible$eventdate <- as.Date(referral_eligible$eventdate,format='%d/%m/%Y')

# extract just the entries that occur before COPD diagnosis and save

referral_before <- referral_eligible[referral_eligible$eventdate < referral_eligible$first_copd,]
save.ffdf(referral_before,dir = '../output/proc_data/referral_before')

print('IMMUNISATION')

# Immunisation
#

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
immunisation <- importFFDF(filename = paste('../data/copd_cohort_march2017/Immunisation.csv',sep=''),datecolnames = '',verbose = T)
save.ffdf(immunisation,dir = '../output/proc_data/immunisation')

# extract only data for eligible and save
immunisation_eligible <- merge(immunisation,eligible,by='patid')
save.ffdf(immunisation_eligible,dir = '../output/proc_data/immunisation_eligible')

# post-hoc turn dates into as.Date

immunisation_eligible$eventdate <- as.character(immunisation_eligible$eventdate)
immunisation_eligible$eventdate <- as.Date(immunisation_eligible$eventdate,format='%d/%m/%Y')

# extract just the entries that occur before COPD diagnosis and save
immunisation_before <- immunisation_eligible[immunisation_eligible$eventdate < immunisation_eligible$first_copd,]
save.ffdf(immunisation_before,dir = '../output/proc_data/immunisation_before')

print('ADDITIONAL')

# Additional
# no dates, so no need to filter

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
additional <- read.csv('../data/copd_cohort_march2017/Additional.csv',header=T)
additional <- as.ffdf(additional)
save.ffdf(additional,dir = '../output/proc_data/additional')  

# extract only data for eligible and save
additional_eligible <- merge(additional,eligible,by='patid')
save.ffdf(additional_eligible,dir = '../output/proc_data/additional_eligible')

print('THERAPY')

# Therapy
#

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
therapy <- importFFDF(filename = paste('../data/copd_cohort_march2017/Therapy_0',1:11,'.csv',sep=''),datecolnames = '',verbose = T)
save.ffdf(therapy,dir = '../output/proc_data/therapy')

# extract only data for eligible and save
therapy_eligible <- merge(therapy,eligible,by='patid')
save.ffdf(therapy_eligible,dir = '../output/proc_data/therapy_eligible')

# post-hoc turn dates into as.Date

therapy_eligible$eventdate <- as.character(therapy_eligible$eventdate)
therapy_eligible$eventdate <- as.Date(therapy_eligible$eventdate,format='%d/%m/%Y')

# extract just the entries that occur before COPD diagnosis and save
therapy_before <- therapy_eligible[therapy_eligible$eventdate < therapy_eligible$first_copd,]
save.ffdf(therapy_before,dir = '../output/proc_data/therapy_before')

