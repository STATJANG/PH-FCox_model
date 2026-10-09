#################################################################
#                   Load clinical covariates.
#################################################################
lgg_id = sapply(strsplit(list.files('./GBM_Linux(FSL)/LGG/Tumor_img_mni152/'),'_'),function(x) x[1])
lgg_tcga_clinical = read.csv("./GBM_R/GBM_clinical_radiomic_data/TCGA-LGG_clinical.csv")
#Selected columns (features)
lgg_tcga_clinical = lgg_tcga_clinical[,c("submitter_id","vital_status","days_to_death","days_to_last_follow_up","gender","age_at_index","treatments_pharmaceutical_treatment_or_therapy")]
lgg_tcga_clinical = dplyr::left_join(data.frame(id=lgg_id), lgg_tcga_clinical,by=c("id"="submitter_id"))
colnames(lgg_tcga_clinical) = c("id","status","os","followup","gender","age","trt")
lgg_tcga_clinical$status= (lgg_tcga_clinical$status=="Dead")+0
lgg_tcga_clinical$trt = (lgg_tcga_clinical$trt=="yes")+0
lgg_tcga_clinical$gender[lgg_tcga_clinical$gender=="not reported"] = NA
#missing 'days to the death' (whether censored or not) values are replaced with the 'days to the last followup'
lgg_tcga_clinical$os[is.na(lgg_tcga_clinical$os)] = lgg_tcga_clinical$followup[is.na(lgg_tcga_clinical$os)]
lgg_tcga_clinical = subset(lgg_tcga_clinical, select = -followup)
lgg_tcga_clinical = lgg_tcga_clinical[!apply(lgg_tcga_clinical,1,function(x) any(is.na(x))),] #Remove subjects with missing values
lgg_tcga_clinical = lgg_tcga_clinical[lgg_tcga_clinical$id!='TCGA-06-0128',] #Vary poor image quality (imperfect image)
lgg_tcga_clinical$gender = (lgg_tcga_clinical$gender=='male')+0
dim(lgg_tcga_clinical)
lgg_tcga_clinical = unique(lgg_tcga_clinical)
dim(lgg_tcga_clinical)


roi_img = oro.nifti::readNIfTI('./GBM_R/MNI-maxprob-thr0-1mm.nii', reorient = FALSE)
lgg_tcga_clinical$roi = NA

for(ii in 1:nrow(lgg_tcga_clinical)){
  tmp_img <- oro.nifti::readNIfTI(paste0('./GBM_Linux(FSL)/LGG/Tumor_img_mni152/',lgg_tcga_clinical$id[ii],'_mni152_tumor.nii.gz'), reorient = FALSE)
  tmp_tab <- table(roi_img[tmp_img==4])
  tmp_tab <- tmp_tab[names(tmp_tab)%in%c(1:9)] # Do not counter 'background(=0))
  if(length(tmp_tab)>0) lgg_tcga_clinical$roi[ii] = names(which.max(tmp_tab))
  print(tmp_tab);cat(ii,'-th subject has been completed with the label of ',lgg_tcga_clinical$roi[ii],'\n')
}
table(lgg_tcga_clinical$roi)
table(lgg_tcga_clinical$roi)/sum(table(lgg_tcga_clinical$roi))
lgg_tcga_clinical$roi_F  = (lgg_tcga_clinical$roi==3)+0
lgg_tcga_clinical$roi_F[is.na(lgg_tcga_clinical$roi_F)]=0
#saveRDS(lgg_tcga_clinical,'./GBM_R/Revision(SIM)/lgg_tcga_clinical.rds')
survminer::ggsurvplot(survfit(Surv(os,status)~roi_F,data=lgg_tcga_clinical))

###Ntumor ratio variable.
ntumor_df = list(0)
for(ii in 1:nrow(lgg_tcga_clinical)){
  tmp_img <- oro.nifti::readNIfTI(paste0('./GBM_Linux(FSL)/LGG/Tumor_img_mni152/',lgg_tcga_clinical$id[ii],'_mni152_tumor.nii.gz'), reorient = FALSE)
  ntumor_df[[ii]] = table(tmp_img)
  print(ntumor_df[[ii]]);cat('  ii is',ii,'\n')
}
sapply(ntumor_df, length)
lgg_tcga_clinical$ntumor = sapply(ntumor_df, function(x) x['4']/sum(x[-1]))
lgg_tcga_clinical$ntumor[is.na(lgg_tcga_clinical$ntumor)] = 0
#saveRDS(lgg_tcga_clinical,'./GBM_R/Revision(SIM)/lgg_tcga_clinical.rds')