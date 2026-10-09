## Compute SEDT3 values
library(oro.nifti);library(imager)
library(RNifti)

#UPENN
for(ii in 1:nrow(gbm_upenn)){
  tmp_img = readNIfTI(paste0('./GBM_Linux(FSL)/UPENN/Tumor_img_mni152/',strsplit(gbm_upenn$id[ii],'_')[[1]][1],'_mni152_tumor.nii.gz'), reorient = FALSE)
  tmp_img[tmp_img==1] = 1
  tmp_img[tmp_img==2] = 1
  tmp_img[tmp_img==3] = 1
  sedt_img = distance_transform(as.cimg(tmp_img==4), value = 0)
  sedt_img = sedt_img - distance_transform(as.cimg(tmp_img==1), value = 0)
  sedt_img[tmp_img==0]=Inf
  writeNifti(sedt_img, paste0('./GBM_R/Revision(SIM)/GBM_UPENN_sedt/',strsplit(gbm_upenn$id[ii],'_')[[1]][1],'_sedt.nii.gz'))
  cat('\n ',ii,'-th subject has been completed.')
}



#UCSF
for(ii in 1:nrow(gbm_ucsf)){
  tmp_img = readNIfTI(paste0('./GBM_Linux(FSL)/UCSF/Tumor_img_mni152/',gbm_ucsf$id[ii],'_mni152_tumor.nii.gz'), reorient = FALSE)
  tmp_img[tmp_img==1] = 1
  tmp_img[tmp_img==2] = 1
  tmp_img[tmp_img==3] = 1
  sedt_img = distance_transform(as.cimg(tmp_img==4), value = 0)
  sedt_img = sedt_img - distance_transform(as.cimg(tmp_img==1), value = 0)
  sedt_img[tmp_img==0]=Inf
  writeNifti(sedt_img, paste0('./GBM_R/Revision(SIM)/GBM_UCSF_sedt/',gbm_ucsf$id[ii],'_sedt.nii.gz'))
  cat('\n ',ii,'-th subject has been completed.')
}


## (older-version) GBM data sets.
load('./GBM_R/GBM_RData/GBM_clinical.RData')
for(ii in 1:nrow(gbm_clinical)){
  tmp_img = readNIfTI(paste0('./GBM_Linux(FSL)/GBM_tumor/',gbm_clinical$id[ii],'_mni152_tumor.nii.gz'), reorient = FALSE)
  tmp_img[tmp_img==1] = 1
  tmp_img[tmp_img==2] = 1
  tmp_img[tmp_img==3] = 1
  sedt_img = distance_transform(as.cimg(tmp_img==4), value = 0)
  sedt_img = sedt_img - distance_transform(as.cimg(tmp_img==1), value = 0)
  sedt_img[tmp_img==0]=Inf
  writeNifti(sedt_img, paste0('./GBM_R/Revision(SIM)/GBM_sedt/',gbm_clinical$id[ii],'_sedt.nii.gz'))
  cat('\n ',ii,'-th subject has been completed.')
}






##Persistent diagram
gbm_ucsf = readRDS('./GBM_R/Revision(SIM)/gbm_ucsf_clinical.rds')
  pdg = vector(mode='list',length=nrow(gbm_ucsf))
  names(pdg) = gbm_ucsf$id
    for(ii in 1:length(pdg)){ 
      #pdg[[ii]] = read.table(paste0('./GBM_R/Revision(SIM)/GBM_UCSF_PH/',
      #                              gbm_ucsf$id[ii],'_pdg.txt'),header=FALSE) 
      pdg[[ii]] = read.table(unz('./GBM_R/Revision(SIM)/GBM_UCSF_PH.zip',
                                 paste0(gbm_ucsf$id[ii],'_pdg.txt')),header=FALSE)
      colnames(pdg[[ii]]) = c('dim','birth','death')
      pdg[[ii]][is.infinite(pdg[[ii]]$death),'death']= pdg[[ii]][is.infinite(pdg[[ii]]$death),'birth']
    }
  pdg_range0 = sapply(pdg, function(x){ c(min(x$birth[x$dim==0]),max(x$death[x$dim==0])) })
  pdg_range0 = c(min(pdg_range0[1,]),max(pdg_range0[2,]))
  pdg_range1 = sapply(pdg, function(x){ c(min(x$birth[x$dim==1]),max(x$death[x$dim==1])) })
  pdg_range1 = c(min(pdg_range1[1,]),max(pdg_range1[2,]))
  pdg_range2 = sapply(pdg, function(x){ c(min(x$birth[x$dim==2]),max(x$death[x$dim==2])) })
  pdg_range2 = c(min(pdg_range2[1,]),max(pdg_range2[2,]))
  pdg_range = rbind(pdg_range0,pdg_range1,pdg_range2)
  pdg_range[,1] = floor(pdg_range[,1])-1
  pdg_range[,2] = ceiling(pdg_range[,2])+1
  
  #base::save(pdg,pdg_range, file = './GBM_R/Revision(SIM)/GBM_UCSF_pdg.RData')
  
  
  
#Upenn  
  gbm_upenn = readRDS('./GBM_R/Revision(SIM)/gbm_upenn_clinical.rds')
  pdg = vector(mode='list',length=nrow(gbm_upenn))
  names(pdg) = gbm_upenn$id
  for(ii in 1:length(pdg)){ 
    #pdg[[ii]] = read.table(paste0('./GBM_R/Revision(SIM)/GBM_upenn_PH/',
    #                              gbm_upenn$id[ii],'_pdg.txt'),header=FALSE) 
    pdg[[ii]] = read.table(unz('./GBM_R/Revision(SIM)/GBM_UPENN_PH.zip',
                               paste0(gbm_upenn$id[ii],'_pdg.txt')),header=FALSE)
    colnames(pdg[[ii]]) = c('dim','birth','death')
    pdg[[ii]][is.infinite(pdg[[ii]]$death),'death']= pdg[[ii]][is.infinite(pdg[[ii]]$death),'birth']
  }
  pdg_range0 = sapply(pdg, function(x){ c(min(x$birth[x$dim==0]),max(x$death[x$dim==0])) })
  pdg_range0 = c(min(pdg_range0[1,]),max(pdg_range0[2,]))
  pdg_range1 = sapply(pdg, function(x){ c(min(x$birth[x$dim==1]),max(x$death[x$dim==1])) })
  pdg_range1 = c(min(pdg_range1[1,]),max(pdg_range1[2,]))
  pdg_range2 = sapply(pdg, function(x){ c(min(x$birth[x$dim==2]),max(x$death[x$dim==2])) })
  pdg_range2 = c(min(pdg_range2[1,]),max(pdg_range2[2,]))
  pdg_range = rbind(pdg_range0,pdg_range1,pdg_range2)
  pdg_range[,1] = floor(pdg_range[,1])-1
  pdg_range[,2] = ceiling(pdg_range[,2])+1
  
  #base::save(pdg,pdg_range, file = './GBM_R/Revision(SIM)/GBM_UPENN_pdg.RData')
  
  
  
#GBM-TCGA
  gbm_tcga_clinical = readRDS('./GBM_R/Revision(SIM)/gbm_tcga_clinical.rds')
  pdg = vector(mode='list',length=nrow(gbm_tcga_clinical))
  names(pdg) = gbm_tcga_clinical$id
  for(ii in 1:length(pdg)){ 
    pdg[[ii]] = read.table(paste0('./GBM_Python/GBM_tda/Features_AT_NonAT/',
                             gbm_tcga_clinical$id[ii],'_pd.txt'),header=FALSE)
    colnames(pdg[[ii]]) = c('dim','birth','death')
    pdg[[ii]][is.infinite(pdg[[ii]]$death),'death']= pdg[[ii]][is.infinite(pdg[[ii]]$death),'birth']
  }
  pdg_range0 = sapply(pdg, function(x){ c(min(x$birth[x$dim==0]),max(x$death[x$dim==0])) })
  pdg_range0 = c(min(pdg_range0[1,]),max(pdg_range0[2,]))
  pdg_range1 = sapply(pdg, function(x){ c(min(x$birth[x$dim==1]),max(x$death[x$dim==1])) })
  pdg_range1 = c(min(pdg_range1[1,]),max(pdg_range1[2,]))
  pdg_range2 = sapply(pdg, function(x){ c(min(x$birth[x$dim==2]),max(x$death[x$dim==2])) })
  pdg_range2 = c(min(pdg_range2[1,]),max(pdg_range2[2,]))
  pdg_range = rbind(pdg_range0,pdg_range1,pdg_range2)
  pdg_range[,1] = floor(pdg_range[,1])-1
  pdg_range[,2] = ceiling(pdg_range[,2])+1
  
  #base::save(pdg,pdg_range, file = './GBM_R/Revision(SIM)/GBM_TCGA_pdg.RData')
  
  