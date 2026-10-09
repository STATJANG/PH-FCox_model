rm(list=ls())
library(spatstat)
library(glmnet)
library(survival)
library(ggplot2)
library(ggthemes)

zip_path = paste0('./GBM_R/Revision(SIM)/loocv_lgg_tcga.zip')
p_df = NULL
for(ll in 1:1000){
  tmp = read.table(unz(zip_path,paste0('loocv_',ll,'.txt')),header=T)
  tmp[is.na(tmp)] = 1
  p_df = rbind(p_df , c(apply(tmp[,1:3],2,min),apply(tmp[,-(1:3)],2,max)))
}
apply(p_df[,1:3],2,min,na.rm=T)
#apply(p_df[,-(1:3)],2,max,na.rm=T)

optimal_sm_vec = c( apply(p_df[,1:3],2,which.min), apply(p_df[,-(1:3)],2,which.max) )

optimal_ll_vec = c(
  which.min(read.table(unz(zip_path,paste0('loocv_',optimal_sm_vec[1],'.txt')),header=T)[,1]),
  which.min(read.table(unz(zip_path,paste0('loocv_',optimal_sm_vec[2],'.txt')),header=T)[,2]),
  which.min(read.table(unz(zip_path,paste0('loocv_',optimal_sm_vec[3],'.txt')),header=T)[,3]))

lgg_dat = readRDS('./GBM_R/Revision(SIM)/lgg_tcga_clinical.rds')
lgg_df = lgg_dat[,c('id','os','status')]
lgg_df$Xc = I(lgg_dat[,c('age', 'gender', 'ntumor', 'roi_F')])

sm_par_mat = expand.grid( seq(.3, 3, length.out=10), seq(.3, 3, length.out=10), seq(.3, 3, length.out=10))
lam_vec = exp(seq(-6,-2,length.out=30)) #guessted from cv.glmnet
load('./GBM_R/Revision(SIM)/LGG_pdg_TCGA.RData')


lgg_fcox_cv = function(sm,ll){
  
  set.seed(666)
  cv_sample = sample(x = 1:nrow(lgg_dat), size = nrow(lgg_dat), replace = FALSE)
  cv_n = floor(nrow(lgg_dat)/10)
  cv_idx = list(0)
  for(ii in 1:9){ cv_idx[[ii]] = cv_sample[(1:cv_n)+cv_n*(ii-1)]}
  cv_idx[[10]] = cv_sample[(cv_n*9+1):length(cv_sample)]
  
  sm_vector = unname(unlist(sm_par_mat[sm,]))
  
  ##PS construction
  ps_list = list(0)
  for(dd in 0:2){
    ps_tmp = NA
    for(ii in 1:nrow(lgg_dat)){
      pixnum = pdg_range[dd+1,2] - pdg_range[dd+1,1]
      ps_grd = seq(pdg_range[dd+1,1]+0.5, pdg_range[dd+1,2]-0.5, length.out = pixnum)
      ps_grd = cbind( rep(ps_grd,each=pixnum), rep(ps_grd,pixnum) )
      
      PD_bd = pdg[[ii]]
      PD_bd = PD_bd[PD_bd$dim == dd, c("birth","death")]
      
      weight_tmp = PD_bd$death - PD_bd$birth
      weight_tmp = abs(cbind(weight_tmp, PD_bd$birth, PD_bd$death))
      weight_tmp = apply(weight_tmp,1,max) 
      
      surface_z = spatstat.explore::density.ppp(ppp(PD_bd[,1], PD_bd[,2], pdg_range[dd+1,], pdg_range[dd+1,]),
                                                weights = weight_tmp, sigma = sm_vector[dd+1], 
                                                dimyx = c(pixnum,pixnum))$v #kernel = "gaussian" (default)
      surface_z = surface_z[ps_grd[,1]<=ps_grd[,2]]
      surface_z[surface_z<0] = 0
      ps_tmp = rbind(ps_tmp,surface_z)
    }
    ps_list[[dd+1]] = ps_tmp[-1,]
    rownames(ps_list[[dd+1]]) = names(pdg)
  }
  ps_list[[1]] = ps_list[[1]][match(lgg_dat$id, rownames(ps_list[[1]])), ]
  ps_list[[2]] = ps_list[[2]][match(lgg_dat$id, rownames(ps_list[[2]])), ]
  ps_list[[3]] = ps_list[[3]][match(lgg_dat$id, rownames(ps_list[[3]])), ]
  
    tmp_risk1 = tmp_risk2 = tmp_risk3 = NA
    
    ##Cross-validation
    for(cc in 1:length(cv_idx)){
      df_tr = lgg_df[-cv_idx[[cc]],]
      colnames(df_tr$Xc) = c('age', 'gender', 'ntumor', 'roi_F')
      df_te = lgg_df[(cv_idx[[cc]]),]
      colnames(df_te$Xc) = c('age', 'gender', 'ntumor', 'roi_F')
      
      #df_te$Xc[,c(1,3)] = I(sweep(df_te$Xc[,c(1,3)],2,colMeans(df_tr$Xc[,c(1,3)]),'-'))
      #df_te$Xc[,c(1,3)] = I(sweep(df_te$Xc[,c(1,3)],2,apply(df_tr$Xc[,c(1,3)],2,sd),'/'))
      #df_tr$Xc[,c(1,3)] = I(scale(df_tr$Xc[,c(1,3)], scale = TRUE, center = TRUE))
      
      df_tr$X0 = I(ps_list[[1]][-cv_idx[[cc]],])
      df_tr$X1 = I(ps_list[[2]][-cv_idx[[cc]],])
      df_tr$X2 = I(ps_list[[3]][-cv_idx[[cc]],])
      df_te$X0 = I(ps_list[[1]][(cv_idx[[cc]]),,drop=F])
      df_te$X1 = I(ps_list[[2]][(cv_idx[[cc]]),,drop=F])
      df_te$X2 = I(ps_list[[3]][(cv_idx[[cc]]),,drop=F])
      
      #Centering for FPC
      df_te$X0 = I(sweep(df_te$X0,2,colMeans(df_tr$X0),'-'))
      df_te$X1 = I(sweep(df_te$X1,2,colMeans(df_tr$X1),'-'))
      df_te$X2 = I(sweep(df_te$X2,2,colMeans(df_tr$X2),'-'))
      
      df_tr$X0 = I(scale(df_tr$X0, center=TRUE, scale=FALSE ))
      df_tr$X1 = I(scale(df_tr$X1, center=TRUE, scale=FALSE ))
      df_tr$X2 = I(scale(df_tr$X2, center=TRUE, scale=FALSE ))
      
      #FPC
      num_pc=20
      svd0 = svd(as.matrix(df_tr$X0/sqrt(nrow(df_tr$X0))), nv=num_pc, nu=num_pc)  
      svd1 = svd(as.matrix(df_tr$X1/sqrt(nrow(df_tr$X1))), nv=num_pc, nu=num_pc)
      svd2 = svd(as.matrix(df_tr$X2/sqrt(nrow(df_tr$X2))), nv=num_pc, nu=num_pc)
      eigenval0 = cumsum(svd0$d^2)/sum(svd0$d^2)  
      eigenval1 = cumsum(svd1$d^2)/sum(svd1$d^2) 
      eigenval2 = cumsum(svd2$d^2)/sum(svd2$d^2) 
      hh0 = min(min(which(eigenval0>0.9)),num_pc)
      hh1 = min(min(which(eigenval1>0.9)),num_pc)
      hh2 = min(min(which(eigenval2>0.9)),num_pc)
      eigenfun0 = svd0$v[,1:hh0]
      eigenfun1 = svd1$v[,1:hh1]
      eigenfun2 = svd2$v[,1:hh2]
      
      df_tr$X0 =  I(as.matrix(df_tr$X0)%*%eigenfun0)
      df_tr$X1 =  I(as.matrix(df_tr$X1)%*%eigenfun1)
      df_tr$X2 =  I(as.matrix(df_tr$X2)%*%eigenfun2)
      colnames(df_tr$X0) = paste0('dim0.',1:hh0)
      colnames(df_tr$X1) = paste0('dim1.',1:hh1)
      colnames(df_tr$X2) = paste0('dim2.',1:hh2)
      
      df_tr$X0i =  I(df_tr$X0*matrix(rep(lgg_dat$roi_F[-cv_idx[[cc]]],hh0),ncol=hh0))
      df_tr$X1i =  I(df_tr$X1*matrix(rep(lgg_dat$roi_F[-cv_idx[[cc]]],hh1),ncol=hh1))
      df_tr$X2i =  I(df_tr$X2*matrix(rep(lgg_dat$roi_F[-cv_idx[[cc]]],hh2),ncol=hh2))
      colnames(df_tr$X0i) = paste0('dim0.',1:hh0,'i')
      colnames(df_tr$X1i) = paste0('dim1.',1:hh1,'i')
      colnames(df_tr$X2i) = paste0('dim2.',1:hh2,'i')
      
      df_te$X0 =  I(as.matrix(df_te$X0)%*%eigenfun0)
      df_te$X1 =  I(as.matrix(df_te$X1)%*%eigenfun1)
      df_te$X2 =  I(as.matrix(df_te$X2)%*%eigenfun2)
      colnames(df_te$X0) = paste0('dim0.',1:hh0)
      colnames(df_te$X1) = paste0('dim1.',1:hh1)
      colnames(df_te$X2) = paste0('dim2.',1:hh2)
      
      df_te$X0i =  I(df_te$X0*matrix(rep(lgg_dat$roi_F[(cv_idx[[cc]])],hh0),ncol=hh0))
      df_te$X1i =  I(df_te$X1*matrix(rep(lgg_dat$roi_F[(cv_idx[[cc]])],hh1),ncol=hh1))
      df_te$X2i =  I(df_te$X2*matrix(rep(lgg_dat$roi_F[(cv_idx[[cc]])],hh2),ncol=hh2))
      colnames(df_te$X0i) = paste0('dim0.',1:hh0,'i')
      colnames(df_te$X1i) = paste0('dim1.',1:hh1,'i')
      colnames(df_te$X2i) = paste0('dim2.',1:hh2,'i')
      
      #Standardization
      df_te$X0 = I(sweep(df_te$X0,2,colMeans(df_tr$X0),'-'))
      df_te$X1 = I(sweep(df_te$X1,2,colMeans(df_tr$X1),'-'))
      df_te$X2 = I(sweep(df_te$X2,2,colMeans(df_tr$X2),'-'))
      df_te$X0 = I(sweep(df_te$X0,2,apply(df_tr$X0,2,sd),'/'))
      df_te$X1 = I(sweep(df_te$X1,2,apply(df_tr$X1,2,sd),'/'))
      df_te$X2 = I(sweep(df_te$X2,2,apply(df_tr$X2,2,sd),'/'))
      
      df_te$X0i = I(sweep(df_te$X0i,2,colMeans(df_tr$X0i),'-'))
      df_te$X1i = I(sweep(df_te$X1i,2,colMeans(df_tr$X1i),'-'))
      df_te$X2i = I(sweep(df_te$X2i,2,colMeans(df_tr$X2i),'-'))
      df_te$X0i = I(sweep(df_te$X0i,2,apply(df_tr$X0i,2,sd),'/'))
      df_te$X1i = I(sweep(df_te$X1i,2,apply(df_tr$X1i,2,sd),'/'))
      df_te$X2i = I(sweep(df_te$X2i,2,apply(df_tr$X2i,2,sd),'/'))
      
      
      df_tr$X0 = I(scale(df_tr$X0, center=TRUE, scale=TRUE ))
      df_tr$X1 = I(scale(df_tr$X1, center=TRUE, scale=TRUE ))
      df_tr$X2 = I(scale(df_tr$X2, center=TRUE, scale=TRUE ))
      df_tr$X0i = I(scale(df_tr$X0i, center=TRUE, scale=TRUE ))
      df_tr$X1i = I(scale(df_tr$X1i, center=TRUE, scale=TRUE ))
      df_tr$X2i = I(scale(df_tr$X2i, center=TRUE, scale=TRUE ))
      
      
      penal_vec = c(rep(0,ncol(df_tr$Xc)-1),rep(1,hh0+hh1+hh2))
      lam.fit1 = glmnet(y=Surv(df_tr$os, df_tr$status),
                        x=as.matrix(df_tr[,c('Xc','X0','X1','X2')])[,-c(4)], #age,gender,ntumor
                        family='cox',alpha=1,standardize = FALSE,
                        penalty.factor = penal_vec, 
                        lambda=lam_vec[ll])
      
      penal_vec = c(rep(0,ncol(df_tr$Xc)),rep(1,hh0+hh1+hh2))
      lam.fit2 = glmnet(y=Surv(df_tr$os, df_tr$status),
                            x=as.matrix(df_tr[,c('Xc','X0','X1','X2')]), #age,gender,ntumor,location
                            family='cox',alpha=1,standardize = FALSE,
                            penalty.factor = penal_vec, 
                            lambda=lam_vec[ll])
      
      penal_vec = c(rep(0,ncol(df_tr$Xc)),rep(1,(hh0+hh1+hh2)*2))
      lam.fit3 = glmnet(y=Surv(df_tr$os, df_tr$status),
                              x=as.matrix(df_tr[,c('Xc','X0','X0i','X1','X1i','X2','X2i')]), #age,gender,ntumor,location 
                              family='cox',alpha=1,standardize = FALSE,
                              penalty.factor = penal_vec, 
                              lambda=lam_vec[ll])
      
      tmp_risk1[(cv_idx[[cc]])] = as.matrix(df_te[,c('Xc','X0','X1','X2')])[,-c(4)]%*%coef(lam.fit1)
      tmp_risk2[(cv_idx[[cc]])] = as.matrix(df_te[,c('Xc','X0','X1','X2')])%*%coef(lam.fit2)
      tmp_risk3[(cv_idx[[cc]])] = as.matrix(df_te[,c('Xc','X0','X0i','X1','X1i','X2','X2i')])%*%coef(lam.fit3)
    }
    
    km_df1 = data.frame(lgg_dat[,c('id','os','status')],est_risk = tmp_risk1)
    km_df2 = data.frame(lgg_dat[,c('id','os','status')],est_risk = tmp_risk2)
    km_df3 = data.frame(lgg_dat[,c('id','os','status')],est_risk = tmp_risk3)
    
    km_df1$risk_group =  km_df1$est_risk<median(km_df1$est_risk)
    km_df2$risk_group =  km_df2$est_risk<median(km_df2$est_risk)
    km_df3$risk_group =  km_df3$est_risk<median(km_df3$est_risk)

  return(list(km_df1,km_df2,km_df3))
}

km_df3 = lgg_fcox_cv(sm = optimal_sm_vec[3],ll = optimal_ll_vec[3])[[3]]
pval3 = survdiff(Surv(time=os,event=status)~risk_group, data = km_df3)$pvalue
km_df1 = lgg_fcox_cv(sm = optimal_sm_vec[1],ll = optimal_ll_vec[1])[[1]]
pval1 = survdiff(Surv(time=os,event=status)~risk_group, data = km_df1)$pvalue
km_df2 = lgg_fcox_cv(sm = optimal_sm_vec[2],ll = optimal_ll_vec[2])[[2]]
pval2 = survdiff(Surv(time=os,event=status)~risk_group, data = km_df2)$pvalue

km_surv1 = survminer::ggsurvplot(survfit(Surv(time=os,event=status)~risk_group, data = km_df1),
                                 conf.int = TRUE,  pval = sprintf("p-value = %.0e", pval1), 
                                 pval.coord = c(1150, 0.9), 
                                 xlim = c(0, max(km_df1$os)),
                                 #legend.labs = c("High","Low"),title="ROI-specific PS")
                                 legend.labs = c("High","Low"),title=NULL)
km_surv2 = survminer::ggsurvplot(survfit(Surv(time=os,event=status)~risk_group, data = km_df2),
                                 conf.int = TRUE,  pval = sprintf("p-value = %.0e", pval2), 
                                 pval.coord = c(1150, 0.9), 
                                 xlim = c(0, max(km_df2$os)),
                                 #legend.labs = c("High","Low"),title="ROI-specific PS")
                                 legend.labs = c("High","Low"),title=NULL)
km_surv3 = survminer::ggsurvplot(survfit(Surv(time=os,event=status)~risk_group, data = km_df3),
                                 conf.int = TRUE,  pval = sprintf("p-value = %.0e", pval3), 
                                 pval.coord = c(1150, 0.9), 
                                 xlim = c(0, max(km_df3$os)),
                                 #legend.labs = c("High","Low"),title="ROI-specific PS")
                                 legend.labs = c("High","Low"),title=NULL)

concordance(coxph(Surv(time=os,event=status)~est_risk, data = km_df1))
concordance(coxph(Surv(time=os,event=status)~est_risk, data = km_df2))
concordance(coxph(Surv(time=os,event=status)~est_risk, data = km_df3))


##Clinical model for TCGA
set.seed(666)
cv_sample = sample(x = 1:nrow(lgg_dat), size = nrow(lgg_dat), replace = FALSE)
cv_n = floor(nrow(lgg_dat)/10)
cv_idx = list(0)
for(ii in 1:9){ cv_idx[[ii]] = cv_sample[(1:cv_n)+cv_n*(ii-1)]}
cv_idx[[10]] = cv_sample[(cv_n*9+1):length(cv_sample)]

risk_vec = NA
for(cc in 1:length(cv_idx)){
  fit.cox = coxph(Surv(os,status)~age+gender+ntumor+roi_F,data=lgg_dat[-cv_idx[[cc]],])
  risk_vec[cv_idx[[cc]]] = predict(fit.cox, lgg_dat[cv_idx[[cc]],,drop=F])
}
km_df_tcga = data.frame(lgg_dat[,c('id','os','status')],risk_group = risk_vec<median(risk_vec)+0)
pval_clin =survdiff(Surv(time=os,event=status)~risk_group, data = km_df_tcga)$pvalue 

km_surv_clin = survminer::ggsurvplot(survfit(Surv(time=os,event=status)~risk_group, data = km_df_tcga),
                                     conf.int = TRUE,  pval = sprintf("p-value = %.0e", pval_clin), 
                                     pval.coord = c(1150, 0.9), 
                                     xlim = c(0, max(km_df_tcga$os)),
                                     #legend.labs = c("High","Low"),title="ROI-specific PS")
                                     legend.labs = c("High","Low"),title=NULL)
# Survival plots
#gbm_surv = ggpubr::ggarrange(
#  km_surv1$plot,km_surv2$plot,km_surv3$plot,
#  ncol = 3, nrow = 1,
#  labels = c("(a)", "(b)", "(c)"),
#  font.label = list(size = 20, face = "bold")
#)

#ggsave("./GBM_R/Revision(SIM)/SIM_lgg_surv_revision.eps", gbm_surv,
#       device = cairo_ps,
#       fallback_resolution = 600, width=12,height=3.5)




sm = optimal_sm_vec[3]
ll = optimal_ll_vec[3]
sm_vector = unname(unlist(sm_par_mat[sm,]))
ps_list = list(0)
for(dd in 0:2){
  ps_tmp = 0
  for(ii in 1:nrow(lgg_dat)){
    pixnum = pdg_range[dd+1,2] - pdg_range[dd+1,1]
    ps_grd = seq(pdg_range[dd+1,1]+0.5, pdg_range[dd+1,2]-0.5, length.out = pixnum)
    ps_grd = cbind( rep(ps_grd,each=pixnum), rep(ps_grd,pixnum) )
    
    PD_bd = pdg[[ii]]
    PD_bd = PD_bd[PD_bd$dim == dd, c("birth","death")]
    
    weight_tmp = PD_bd$death - PD_bd$birth
    weight_tmp = abs(cbind(weight_tmp, PD_bd$birth, PD_bd$death))
    weight_tmp = apply(weight_tmp,1,max) 
    
    surface_z = spatstat.explore::density.ppp(ppp(PD_bd[,1], PD_bd[,2], pdg_range[dd+1,], pdg_range[dd+1,]),
                                              weights = weight_tmp, sigma = sm_vector[dd+1], 
                                              dimyx = c(pixnum,pixnum))$v #kernel = "gaussian" (default)
    surface_z = surface_z[ps_grd[,1]<=ps_grd[,2]]
    surface_z[surface_z<0] = 0
    ps_tmp = rbind(ps_tmp,surface_z)
  }
  ps_list[[dd+1]] = ps_tmp[-1,]
  rownames(ps_list[[dd+1]]) = names(pdg)
}
ps_list[[1]] = ps_list[[1]][match(lgg_dat$id, rownames(ps_list[[1]])), ]
ps_list[[2]] = ps_list[[2]][match(lgg_dat$id, rownames(ps_list[[2]])), ]
ps_list[[3]] = ps_list[[3]][match(lgg_dat$id, rownames(ps_list[[3]])), ]

lgg_df = lgg_dat[,c('id','os','status')]
lgg_df$Xc = I(lgg_dat[,c('age', 'gender',  'ntumor', 'roi_F')])
lgg_df$Xc = I(scale(lgg_df$Xc, scale = TRUE, center = TRUE))

#Centering for FPC
lgg_df$X0 = I(scale(ps_list[[1]], center=TRUE, scale=FALSE ))
lgg_df$X1 = I(scale(ps_list[[2]], center=TRUE, scale=FALSE ))
lgg_df$X2 = I(scale(ps_list[[3]], center=TRUE, scale=FALSE ))

#FPC
num_pc=20
svd0 = svd(as.matrix(lgg_df$X0/sqrt(nrow(lgg_df$X0))), nv=num_pc, nu=num_pc)  
svd1 = svd(as.matrix(lgg_df$X1/sqrt(nrow(lgg_df$X1))), nv=num_pc, nu=num_pc)
svd2 = svd(as.matrix(lgg_df$X2/sqrt(nrow(lgg_df$X2))), nv=num_pc, nu=num_pc)
eigenval0 = cumsum(svd0$d^2)/sum(svd0$d^2)  
eigenval1 = cumsum(svd1$d^2)/sum(svd1$d^2) 
eigenval2 = cumsum(svd2$d^2)/sum(svd2$d^2) 
hh0 = min(min(which(eigenval0>0.9)),num_pc)
hh1 = min(min(which(eigenval1>0.9)),num_pc)
hh2 = min(min(which(eigenval2>0.9)),num_pc)
eigenfun0 = svd0$v[,1:hh0]
eigenfun1 = svd1$v[,1:hh1]
eigenfun2 = svd2$v[,1:hh2]

lgg_df$X0 =  I(as.matrix(lgg_df$X0)%*%eigenfun0)
lgg_df$X1 =  I(as.matrix(lgg_df$X1)%*%eigenfun1)
lgg_df$X2 =  I(as.matrix(lgg_df$X2)%*%eigenfun2)
colnames(lgg_df$X0) = paste0('dim0.',1:hh0)
colnames(lgg_df$X1) = paste0('dim1.',1:hh1)
colnames(lgg_df$X2) = paste0('dim2.',1:hh2)

lgg_df$X0i =  I(lgg_df$X0*matrix(rep(lgg_dat$roi_F,hh0),ncol=hh0))
lgg_df$X1i =  I(lgg_df$X1*matrix(rep(lgg_dat$roi_F,hh1),ncol=hh1))
lgg_df$X2i =  I(lgg_df$X2*matrix(rep(lgg_dat$roi_F,hh2),ncol=hh2))
colnames(lgg_df$X0i) = paste0('dim0.',1:hh0,'i')
colnames(lgg_df$X1i) = paste0('dim1.',1:hh1,'i')
colnames(lgg_df$X2i) = paste0('dim2.',1:hh2,'i')

lgg_df$X0 = I(scale(lgg_df$X0, center=TRUE, scale=TRUE ))
lgg_df$X1 = I(scale(lgg_df$X1, center=TRUE, scale=TRUE ))
lgg_df$X2 = I(scale(lgg_df$X2, center=TRUE, scale=TRUE ))
lgg_df$X0i = I(scale(lgg_df$X0i, center=TRUE, scale=TRUE ))
lgg_df$X1i = I(scale(lgg_df$X1i, center=TRUE, scale=TRUE ))
lgg_df$X2i = I(scale(lgg_df$X2i, center=TRUE, scale=TRUE ))

penal_vec = c(rep(0,ncol(lgg_df$Xc)),rep(1,(hh0+hh1+hh2)*2))
lam.fit = glmnet(y=Surv(lgg_df$os, lgg_df$status),
                 x=as.matrix(lgg_df[,c('Xc','X0','X0i','X1','X1i','X2','X2i')]),
                 family='cox',alpha=1,standardize = FALSE,
                 penalty.factor = penal_vec, 
                 lambda=lam_vec[ll])

X0_coef = coef(lam.fit)[grep("X0\\.", rownames(coef(lam.fit)))]      
X1_coef = coef(lam.fit)[grep("X1\\.", rownames(coef(lam.fit)))]      
X2_coef = coef(lam.fit)[grep("X2\\.", rownames(coef(lam.fit)))] 
X0i_coef = coef(lam.fit)[grep("X0i\\.", rownames(coef(lam.fit)))]      
X1i_coef = coef(lam.fit)[grep("X1i\\.", rownames(coef(lam.fit)))]      
X2i_coef = coef(lam.fit)[grep("X2i\\.", rownames(coef(lam.fit)))]


coef_list = list(X0=eigenfun0%*%(X0_coef/svd0$d[1:hh0]),
                 X1=eigenfun1%*%(X1_coef/svd1$d[1:hh1]),
                 X2=eigenfun2%*%(X2_coef/svd2$d[1:hh2]))
coef_list_i = list(X0i=eigenfun0%*%(X0i_coef/svd0$d[1:hh0]),
                   X1i=eigenfun1%*%(X1i_coef/svd1$d[1:hh1]),
                   X2i=eigenfun2%*%(X2i_coef/svd2$d[1:hh2]))      

fcoef_list = list(0)      
for(dd in 0:2){
  fcoef_list[[dd+1]] = vector(mode='list',length = 2)
  names(fcoef_list[[dd+1]]) = c('non_frontal','frontal')
  
  pixnum = pdg_range[dd+1,2] - pdg_range[dd+1,1]
  ps_grd = seq(pdg_range[dd+1,1]+0.5, pdg_range[dd+1,2]-0.5, length.out = pixnum)
  ps_grd = cbind( rep(ps_grd,each=pixnum), rep(ps_grd,pixnum) )
  
  coef_df = data.frame(ps_grd[ps_grd[,1]<=ps_grd[,2],],
                       non_frontal = coef_list[[dd+1]], 
                       frontal = coef_list[[dd+1]]+coef_list_i[[dd+1]] )
  colnames(coef_df) = c('x','y','non_frontal','frontal')
  coef_lim = max(abs(coef_df[,3:4]))
  coef_lim = c(-1,1)*coef_lim  
  
  fcoef_list[[dd+1]]$'non_frontal' = ggplot(coef_df,aes(x,y,fill = non_frontal))+
    geom_raster()+
    theme_tufte()+
    labs(x="birth",y="death",title = NULL,fill=NULL) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(slope=1,intercept=0)+
    scale_fill_gradient2(limits = coef_lim ) +
    theme(legend.position = c(0.9, 0.25))+
    theme(plot.title = element_text(hjust = 0.4))
  
  fcoef_list[[dd+1]]$'frontal' = ggplot(coef_df,aes(x,y,fill = frontal))+
    geom_raster()+
    theme_tufte()+
    labs(x="birth",y="death",title = NULL, fill=NULL) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(slope=1,intercept=0)+
    scale_fill_gradient2(limits = coef_lim ) +
    theme(legend.position = c(0.9, 0.25))+
    theme(plot.title = element_text(hjust = 0.4))
}      

gg_fcoef <- ggpubr::ggarrange(
  fcoef_list[[1]]$non_frontal, fcoef_list[[2]]$non_frontal,fcoef_list[[3]]$non_frontal, 
  fcoef_list[[1]]$frontal,fcoef_list[[2]]$frontal,fcoef_list[[3]]$frontal,
  labels = c("(a)", "(c)", "(e)","(b)", "(d)", "(f)"),
  font.label = list(size = 20, face = "bold"),
  nrow = 2,ncol = 3)

#ggplot2::ggsave("./GBM_R/Revision(SIM)/lgg_fcoef_tcga.eps", gg_fcoef,        
#                device = cairo_ps,
#                fallback_resolution = 600, width=10,height=6.6)



#base::save.image('./GBM_R/Revision(SIM)/LGG_result_TCGA(cv).RData')
