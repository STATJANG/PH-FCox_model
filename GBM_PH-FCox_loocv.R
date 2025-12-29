gbm_fitting = function(sm){
  sm_vector = unname(unlist(sm_par_mat[sm,]))
  gbm_ps = persistent_surface_list(PD_list = gbm_pd, minmax_table = persistent_range(gbm_pd), weight_type='mdw', sd_vector = sm_vector)
  gbm_fp = E_score(gbm_ps$surface, max_dim = 2)
  hh.whole =  sapply(gbm_fp$propEV, function(x){ min(min(which(x>0.9)),fpc.max.order) })
  
  #   Functional features
  fp_mat = data.frame(gbm_fp$Eigenscore[[1]][,1:hh.whole[1]], gbm_fp$Eigenscore[[2]][,1:hh.whole[2]], gbm_fp$Eigenscore[[3]][,1:hh.whole[3]])
  fp_mat_F  = data.frame(fp_mat, fp_mat*as.data.frame( matrix( rep(gbm_dat$'roi_F' ,ncol(fp_mat)), ncol = ncol(fp_mat), nrow = nrow(fp_mat))) )
  
  fp_mat_F$id  = rownames(fp_mat_F)
  
  tmp_dat_F  = dplyr::left_join(gbm_dat[,c('id','status','os','gender','age','ntumor','ks','roi_F' )], fp_mat_F  , by=c('id'='id'))
  
  
  
  
  glmnet.fcox.cv_F <- glmnet::cv.glmnet(
    x=as.matrix(tmp_dat_F[,-(1:3)]),
    y=Surv(time = tmp_dat_F$os, event = tmp_dat_F$status),
    penalty.factor = c(0,0,0,0,0,rep(1,ncol(tmp_dat_F)-8)),
    family = 'cox', nfolds = 10, #default=10
    standardize = FALSE, #type.measure = "C"( default = partial likelihood based deviance )
    alpha=1) #Lasso
  
  glmnet.fcox.fit_F = glmnet::glmnet(x=as.matrix(tmp_dat_F[,-(1:3)]),
                                     y=Surv(time = tmp_dat_F$os, event = tmp_dat_F$status),
                                     penalty.factor = c(0,0,0,0,0,rep(1,ncol(tmp_dat_F)-8)), #roi2 is added now.
                                     family = 'cox',
                                     standardize = FALSE,
                                     lambda = glmnet.fcox.cv_F$lambda.min,
                                     alpha=1)  ##Lasso        
  
  
  ####   LOOCV for selecting smoothing parameters.        ####   
  tmp.risk_F = 0
  
  for(ll in 1:nrow(tmp_dat_F)){
    
    ### Construct functional predictor (common for FT, F and T)  ####
    gbm_ps.train = gbm_ps
    gbm_ps.test = vector(mode='list',length = 3)
    for(dd in 1:3){
      gbm_ps.train$surface[[dd]]  = gbm_ps$surface[[dd]][-ll,]
      gbm_ps.test[[dd]]  = gbm_ps$surface[[dd]][ll,,drop=F] - colMeans(gbm_ps$surface[[dd]][-ll,])
    }
    gbm_fp.train = E_score(gbm_ps.train$surface, max_dim = 2) #, trunc_order = c(5,5,5)
    hh.whole =  sapply(gbm_fp.train$propEV, function(x){ min(min(which(x>0.9)),fpc.max.order) })
    
    for(dd in 1:3){
      gbm_ps.test[[dd]] = gbm_ps.test[[dd]] %*% gbm_fp.train$Eigenfunc[[dd]][,1:hh.whole[dd]]
    }
    
    #  functional features (train)
    fp_mat_train = data.frame(
      gbm_fp.train$Eigenscore[[1]][,1:hh.whole[1]],
      gbm_fp.train$Eigenscore[[2]][,1:hh.whole[2]],
      gbm_fp.train$Eigenscore[[3]][,1:hh.whole[3]])
    fp_mat_train_F  = data.frame(fp_mat_train, fp_mat_train*as.data.frame( matrix( rep(gbm_dat$'roi_F'[-ll] ,ncol(fp_mat_train)), ncol = ncol(fp_mat_train), nrow = nrow(fp_mat_train))) )
    
    fp_mat_train_F$id  = rownames(fp_mat_train_F)
    
    #  functional feature data frame
    train_dat_F  = dplyr::right_join(gbm_dat[,c('id','status','os','gender','age','ntumor','ks','roi_F')],  fp_mat_train_F , by=c('id'='id'))
    
    
    loocv.fit_F = glmnet::glmnet(x=as.matrix(train_dat_F[,-(1:3)]),
                                 y=Surv(time = train_dat_F$os, event = train_dat_F$status),
                                 penalty.factor = c(0,0,0,0,0,rep(1,ncol(train_dat_F)-8)),
                                 family = 'cox',
                                 standardize = FALSE,
                                 lambda = glmnet.fcox.cv_F$lambda.min,
                                 alpha=1  )##Lasso
    
    test_dat_F   = cbind( gbm_dat[ll, c('gender','age','ntumor','ks','roi_F' ),drop=F] ,matrix(c(unlist(lapply(gbm_ps.test, unlist)),unlist(lapply(gbm_ps.test, unlist))*gbm_dat$'roi_F'[ll] ), nrow=1))
    
    tmp.risk_F[ll]    = predict(loocv.fit_F ,  newx = as.matrix(test_dat_F) ,type = 'link') 
    
    #print(ll)
  }
  
  dat_km_F  = data.frame(tmp_dat_F[ ,c('os','status')],group = tmp.risk_F  < median(tmp.risk_F ))
  
  
  tryCatch(
    {
      loocv.pvalue_F  = survival::survdiff(Surv(os, status) ~ group, data = dat_km_F )$pvalue
      
      return(c( loocv.pvalue_F , 
                glmnet.fcox.cv_F$lambda.min))
    },
    error = function(err) {
      message("Caught an error:", conditionMessage(err))
      return(c( NA,
                glmnet.fcox.cv_F$lambda.min)) }
  )
  #### The end of LOOCV
}




source('./GBM/GBM_functions.R')
load('./GBM/GBM_clinical_data.RData')
load('./GBM/GBM_pd.RData')

gbm_dat = gbm_clinical[,c('id','status','os','gender','age','ntumor_AT','ks')]
colnames(gbm_dat)[6] = 'ntumor'
gbm_dat$ntumor = gbm_dat$ntumor/max(gbm_dat$ntumor)

tumor_locs = read.table('./GBM/gbm_tumor_locs.txt',header=T)
tumor_locs = dplyr::left_join(gbm_dat, tumor_locs, by=c('id'='id') )
gbm_dat$'roi_FT' = as.numeric(tumor_locs$roi2 %in% c(3,8))
gbm_dat$'roi_T' = as.numeric(tumor_locs$roi2 %in% c(8))
gbm_dat$'roi_F' = as.numeric(tumor_locs$roi2 %in% c(3))
tumor_locs = dplyr::left_join(gbm_dat, tumor_locs, by=c('id'='id') )


fpc.max.order = 10   
sm_par_mat = expand.grid( seq(.3, 3, length.out=10), seq(.3, 3, length.out=10), seq(.3, 3, length.out=10))
prescv = 0 
for(ss in 1:nrow(sm_par_mat)){
  tmp_fit = gbm_fitting(sm = ss)
  prescv = rbind(prescv , tmp_fit);print(ss)
}
prescv = prescv[-1,]

#save(gbm_dat, file = './GBM/GBM_dataset.RData')
#save(prescv,sm_par_mat, file = './GBM/GBM_PH_FCox_loocv.RData')