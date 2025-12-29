rm(list=ls())

simul_size = 300
n_size = 70 # grp_A 70, grp_B 70, a total of 140
source('./GBM/GBM_functions.R')

load_pd = function(gbm_path, type="binary"){
  if(!type%in%c('binary','ternary')){ stop(' type should be either binary or ternary. ') }
  gbm_pd_files = list.files(gbm_path, pattern = "_pd\\.txt$")
  
  #The size of MNI152 coordinate system.
  x_size = 182 ;  y_size = 218 ;  z_size = 182
  gbm_pd_list = vector(mode='list', length=simul_size) #repetition size
  names(gbm_pd_list) = paste0('iter',1:simul_size)
  
  for(iter in 1:simul_size){ # simulation size (number of repetition)
    gbm_pd_list[[iter]] = vector(mode='list', length=n_size*2) #control 100, treatment 100
    names(gbm_pd_list[[iter]]) = c( paste0('control',1:n_size), paste0('feature',1:n_size) )
    
    
    pd_df_iter = read.table( paste0(gbm_path,'/A_pd',iter,'.txt'), header=T)
    colnames(pd_df_iter) = c('iter','dim','birth','death')
    for(num in 1:n_size){ #sample size ( number of patients)
      pd_df = pd_df_iter[pd_df_iter$iter==num,-1]
      #pd_df = read.table( paste0(gbm_path,'/control_iter_',iter-1,'_num_',num-1,'_pd.txt'), col.names = c('dim','birth','death'))
      if(type=="binary"){ 
        pd_df$death[is.infinite(pd_df$death)] = pd_df$birth[is.infinite(pd_df$death)] 
        }
      if(type=="ternary"){ 
        pd_df = pd_df[-min(which(pd_df$dim==0)),]
        pd_df$death[is.infinite(pd_df$death)] = pd_df$birth[is.infinite(pd_df$death)]
      }
      gbm_pd_list[[iter]][[num]] = pd_df
    }
    
    pd_df_iter = read.table( paste0(gbm_path,'/B_pd',iter,'.txt'), header=T)
    colnames(pd_df_iter) = c('iter','dim','birth','death')
    for(num in 1:n_size){ #sample size ( number of trt and number of ctr)
      pd_df = pd_df_iter[pd_df_iter$iter==num,-1]
      #pd_df = read.table( paste0(gbm_path,'/feature_iter_',iter-1,'_num_',num-1,'_pd.txt'), col.names = c('dim','birth','death'))
      if(type=="binary"){ 
        pd_df$death[is.infinite(pd_df$death)] = pd_df$birth[is.infinite(pd_df$death)] 
        }
      if(type=="ternary"){ 
        pd_df = pd_df[-min(which(pd_df$dim==0)),]
        pd_df$death[is.infinite(pd_df$death)] = pd_df$birth[is.infinite(pd_df$death)]
      }
      gbm_pd_list[[iter]][[n_size+num]] = pd_df
    }
    
  }
  return(gbm_pd_list)
}

#gbm_pd[[iter]][[patients]]
gbm_pd = load_pd(gbm_path = "./simulation/simul1_pd", type = 'binary')





#To take an average of functional coef plots later, the same range of ps's should be unified.
range.0 = t(sapply(lapply(gbm_pd,persistent_range),function(x) x[1,]))
range.1 = t(sapply(lapply(gbm_pd,persistent_range),function(x) x[2,]))
simul.range = rbind(c(min(range.0[,1]),max(range.0[,2])),c(min(range.1[,1]),max(range.1[,2])))

##Persistent surface for each subject
gbm_ps = vector(mode='list', length = simul_size)
for(ii in 1:simul_size){ 
  tmp = persistent_surface_list(PD_list = gbm_pd[[ii]], dim_vector = c(0,1), minmax_table = simul.range, weight_type='mdw', sd_vector = c(2,2))
  save(tmp, file = paste0('./simulation/simul1_ps/ps',ii,'.RData'))
  print(ii)
}

##Eigen decomposition for making functional predictors
for(ii in 1:simul_size){ 
  env_ps = new.env()
  tmp_ps = load(file = paste0('./simulation/simul1_ps/ps',ii,'.RData') ,envir = env_ps)
  tmp_ps = env_ps$tmp
  tmp = E_score(tmp_ps$surface, max_dim = 1)
  
  save(tmp, file = paste0('./simulation/simul1_ps/eigen_',ii,'.RData'))
  print(ii)
}
rm(env_ps)








# Generating survival times (Exponential dist assumed!)
rm(list=ls())
simul_size = 300
n_size = 70 # grp_A 70, grp_B 70, a total of 140

f_ind = rep( c(rep(1, round(n_size*0.3) ),rep(0,n_size - round(n_size*0.3) )), 2) #frontal lobe
simul.B = rep(c(0,1),c(n_size,n_size)) # 1 if a subject belongs to the shape group P

true_coef =  c(0.8,0.4, -0.2) #\beta_{B,frontal} , \beta_{b,non-frontal}, \beta_{A,frontal} ( \beta_{A,nf} = 0)
true_coef2 = c(0.2, 0.5) #beta_{frontal}=0.2, beta_{B}=0.5 (beta_{A} = 0, beta_{non_frontal} = 0)

simul_df = vector(mode='list', length=simul_size)
delta.rate = 0.85 # % of not being censored.
set.seed(7)
for(simul in 1:simul_size){
  simul.eta = cbind(simul.B*f_ind, simul.B*(!f_ind), (!simul.B)*f_ind)%*%true_coef
  simul.eta2 = cbind(f_ind , simul.B)%*%true_coef2
  
  rand.u = runif(n_size*2)
  simul.t = -log(rand.u)/exp(simul.eta)*3000 #lambda=1/3000
  simul.t2 = -log(rand.u)/exp(simul.eta2)*3000
  
  ## status = status2
  rand.u.c = runif(n_size*2)
  simul.censor  = -log(rand.u.c)/( delta.rate/(1-delta.rate)*exp(simul.eta ))*3000
  simul.censor2 = -log(rand.u.c)/( delta.rate/(1-delta.rate)*exp(simul.eta2))*3000
  simul.status  = (simul.censor<=simul.t)+0
  simul.status2 = (simul.censor2<=simul.t2)+0 
  
  simul.t[simul.status==0] = simul.censor [simul.status==0]
  simul.t2[simul.status2==0] = simul.censor2[simul.status2==0]
  
  #os.null = survival times generated without an interaction between shape and frontal_lobe
  simul_df[[simul]] = data.frame(os = simul.t, os.null = simul.t2, status= simul.status, status2= simul.status2,
                            f.B = f_ind*simul.B, nf.B = (!f_ind)*simul.B, f.A = f_ind*(!simul.B), nf.A = (!f_ind)*(!simul.B))
}
mean(sapply(simul_df, function(x) mean(x$status)))  #censoring rate.
#save(simul_df, simul_size, n_size, file = './simulation/simul1_simuldf.RData')



###Making matrices of functional predictors.
f_ind = c(rep(1,n_size), c(rep(1, round(n_size*0.3) ),rep(0,n_size - round(n_size*0.3) )) ) #frontal lobe
fp_mat = vector(mode='list', length = simul_size)

for(ii in 1:simul_size){ 
  env_ps = new.env()
  tmp = load(paste0('./GBM_R/gbm_simul/simul1_ps/eigen_',ii,'.RData'), envir = env_ps)
  tmp = env_ps$tmp
  
  hh.whole_tmp =  sapply(tmp$propEV, function(x){ min(which(x>0.9)) } )
  fp_mat_tmp = data.frame(tmp$Eigenscore[[1]][,1:hh.whole_tmp[1]], 
                          tmp$Eigenscore[[2]][,1:hh.whole_tmp[2]])
  fp.colnames = c( paste0('dim0.',1:hh.whole_tmp[1]), paste0('dim1.',1:hh.whole_tmp[2]) )
  colnames(fp_mat_tmp) = fp.colnames
  
  ##Cell-means model
  fp_mat_tmp.f = fp_mat_tmp*as.data.frame( matrix( rep( f_ind ,ncol(fp_mat_tmp)), ncol = ncol(fp_mat_tmp), nrow = nrow(fp_mat_tmp)))
  colnames(fp_mat_tmp.f) = paste0(fp.colnames,'.f')
  
  fp_mat[[ii]] = data.frame(fp_mat_tmp, fp_mat_tmp.f)
  
  print(ii)
}
rm(env_ps)
save(fp_mat,file='./GBM_R/gbm_simul/gbm_simul1_fp.RData')


