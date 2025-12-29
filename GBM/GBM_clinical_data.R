
##Extracting the ID's of the patients' in GBM_Test
gbm_list = list.files("./GBM/Features_GBM", pattern = "_pd\\.txt$")
gbm_id = sapply(gbm_list, strsplit, split = "_")
gbm_id = unname(sapply(gbm_id, function(x) x[1]))


#################################################################
#                   Load clinical covariates.
#################################################################
gbm_clinical = read.csv("./GBM/TCGA-GBM_clinical.csv")
#Selected columns (features)
gbm_clinical = gbm_clinical[,c("submitter_id","vital_status","days_to_death","days_to_last_follow_up","gender","age_at_index","treatments_pharmaceutical_treatment_or_therapy")]
gbm_clinical = dplyr::left_join(data.frame(id=gbm_id), gbm_clinical,by=c("id"="submitter_id"))
colnames(gbm_clinical) = c("id","status","os","followup","gender","age","trt")
gbm_clinical$status= (gbm_clinical$status=="Dead")+0
gbm_clinical$trt = (gbm_clinical$trt=="yes")+0
gbm_clinical$gender[gbm_clinical$gender=="not reported"] = NA
#missing 'days to the death' (whether censored or not) values are replaced with the 'days to the last followup'
gbm_clinical$os[is.na(gbm_clinical$os)] = gbm_clinical$followup[is.na(gbm_clinical$os)]
gbm_clinical = subset(gbm_clinical, select = -followup)
gbm_clinical = gbm_clinical[!apply(gbm_clinical,1,function(x) any(is.na(x))),] #Remove subjects with missing values
gbm_clinical = gbm_clinical[gbm_clinical$id!='TCGA-06-0128',] #Vary poor image quality (imperfect image)
gbm_id = gbm_clinical$id

load('./GBM/GBM_ntumor.RData')
gbm_clinical = dplyr::left_join(gbm_clinical,ntumor_df, by=c('id'='id'))

#compute the distance from the center
gbm_tumor_locs = read.table('./GBM/gbm_tumor_locs.txt',header=T)
colnames(gbm_tumor_locs) = c('id','tumor_x','tumor_y','tumor_z')
gbm_tumor_locs$dist = sqrt((gbm_tumor_locs$tumor_x+1)^2+(gbm_tumor_locs$tumor_y+17)^2+(gbm_tumor_locs$tumor_z-18)^2)
gbm_tumor_locs = gbm_tumor_locs[,c('id','dist')]
gbm_clinical = dplyr::left_join(gbm_clinical,gbm_tumor_locs, by=c('id'='id'))

gbm_clinical$gender = (gbm_clinical$gender=='male')+0


######Karnofsky score
load('./GBM/GBM_clinical_data.RData')
tmp_gbm = read.csv('GBM_Table_1.csv')
tmp_gbm = tmp_gbm[,c('patient_ID','OS.time','gender','age','karnofsky_score')]
tmp_gbm = tmp_gbm[tmp_gbm$patient_ID%in%gbm_clinical$id,]

tmp_gbm_check = dplyr::inner_join(gbm_clinical, tmp_gbm, by = c('id' = 'patient_ID' ))
#sanity check
#tmp_gbm_check[tmp_gbm_check$age.x != tmp_gbm_check$age.y,]
#tmp_gbm_check[tmp_gbm_check$OS.time != tmp_gbm_check$os,] #one obs
#tmp_gbm_check[(tmp_gbm_check$gender.y=='MALE') != tmp_gbm_check$gender.x,]

tmp_gbm2 = dplyr::left_join(gbm_clinical, tmp_gbm[,c('patient_ID','karnofsky_score')],by = c('id' = 'patient_ID' ))
gbm_clinical = tmp_gbm2
colnames(gbm_clinical)[ncol(gbm_clinical)] = 'ks'

#save(gbm_clinical,file = './GBM/GBM_clinical_data.RData')
