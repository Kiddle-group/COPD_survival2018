print('Loading data')

print('..patient info')
all_copd <- read.dta13("../data/other_files/Extract_files/first_copd_222790.dta") # First COPD diagnosis for each patient, from Quint group
patient <- read.dta13('../data/other_files/Extract_files/pat_pract_222790.dta') # data on patient demographics, from Quint group
copd_patient <- merge(patient,all_copd,by='patid') # merge

# generate derived variables used later
copd_patient[,'year35'] <- copd_patient[,'yob'] + 35 # Year at which patient is 35 years old
copd_patient[,'year_copd'] <- as.numeric(substr(copd_patient[,'first_copd'],1,4)) # year of first COPD diagnosis

print('..linkage eligibility')

# Not everybody is eligible for ONS & HES & LSOA
le <- read.table(file="../data/Final + new linkage/Final_cohort/linkage_eligibility.txt",sep="\t",header=T) # check ONS and HES linkage eligibility
le2 <- le[le[,4] & le[,5] & le[,7],] # remove individuals not eligible to link
linked2 <- merge(le2,copd_patient,by="patid") #merge this file with COPD data
linked <- copd_patient$patid %in% linked2$patid #so that you can tell which COPD patients can be linked

rm(le,le2,patient,linked2) #remove tmp files

# Define my own start id, max of current registration date, up to standard date, age 35+

tmp <- copd_patient[,c('dob','crd','uts')]+c(rep(35*365,dim(copd_patient)[1]),rep(0,dim(copd_patient)[1]*2))

copd_patient$startid <- apply(tmp, 1, 'max')

# identify patients who transferred out of practice before start ID

out_before_in <- which(as.Date(copd_patient$startid) > as.Date(copd_patient$tod)) 

copd_patient$earlyout <- 0
copd_patient$earlyout[out_before_in] <- 1

rm(out_before_in,tmp) #remove tmp files

# define eligible COPD patients (eligible for linkage, not dead before COPD diagnosis, eligible for start before end)
#eligible <- linked  & (as.Date(copd_patient$startid) < as.Date('2016-02-01')) & (as.Date(copd_patient$startid) < as.Date(copd_patient$first_copd)) & copd_patient$earlyout==0 & copd_patient$year_copd >= 2004 
eligible <- linked  & (as.Date(copd_patient$first_copd) <= as.Date('2012-09-19')) & (as.Date(copd_patient$startid) < as.Date(copd_patient$first_copd)) & copd_patient$earlyout==0 & copd_patient$year_copd >= 2004 


print('..define cohort - except early death')

copd_eligible <- copd_patient[eligible,] #data on just eligible patients


rm(copd_patient) #remove tmp file



