### Calculating and drawing functional coefficients. 
rm(list=ls())

simul_size = 300
n_size = 70 # grp_A 70, grp_B 70, a total of 140

##  Drawing Functional coefficients  ##
load(paste0('./GBM_R/gbm_simul/gbm_simul1_fp.RData'))
load(paste0('./GBM_R/gbm_simul/gbm_simul1.RData'))

env_ps = new.env()
ps_grid = load(file = paste0('./GBM_R/gbm_simul/simul1_ps/ps',1,'.RData') ,envir = env_ps)
ps_grid = env_ps$tmp$grid
rm(env_ps)

fcoef  = vector(mode = 'list',length=simul_size)
fcoef.df  = vector(mode = 'list',length=simul_size)
fcoef2  = vector(mode = 'list',length=simul_size)
fcoef2.df  = vector(mode = 'list',length=simul_size)

for(simul in 1:simul_size){
  
  fcoef[[simul]]  = vector(mode = 'list',length=2) # 0th and 1st dim.
  fcoef.df[[simul]]  = vector(mode = 'list',length=2)
  names(fcoef[[simul]]) = names(fcoef.df[[simul]])  = c('dim0','dim1')
  
  fcoef2[[simul]]  = vector(mode = 'list',length=2) # 0th and 1st dim.
  fcoef2.df[[simul]]  = vector(mode = 'list',length=2)
  names(fcoef2[[simul]]) = names(fcoef2.df[[simul]])  = c('dim0','dim1')
  
  
  
  for( dd in 0:1){  # topological dim
    fcoef[[simul]][[dd+1]]  = vector(mode = 'list',length=3) # c('frontal', 'non_frontal', 'diff')
    fcoef.df[[simul]][[dd+1]]  = vector(mode = 'list',length=3)
    names(fcoef[[simul]][[dd+1]]) = names(fcoef.df[[simul]][[dd+1]]) = c('frontal', 'non_frontal', 'diff')
    
    fcoef.df[[simul]][[dd+1]]$frontal  = data.frame(x= ps_grid[[dd+1]][,1] , y =ps_grid[[dd+1]][,2], z = 0 )
    fcoef.df[[simul]][[dd+1]]$non_frontal  = data.frame(x= ps_grid[[dd+1]][,1] , y =ps_grid[[dd+1]][,2], z = 0 )
    fcoef.df[[simul]][[dd+1]]$diff  = data.frame(x= ps_grid[[dd+1]][,1] , y =ps_grid[[dd+1]][,2], z = 0 )
    
    fcoef2[[simul]][[dd+1]]  = vector(mode = 'list',length=3) # c('frontal', 'non_frontal', 'diff')
    fcoef2.df[[simul]][[dd+1]]  = vector(mode = 'list',length=3)
    names(fcoef2[[simul]][[dd+1]]) = names(fcoef2.df[[simul]][[dd+1]]) = c('frontal', 'non_frontal', 'diff')
    
    fcoef2.df[[simul]][[dd+1]]$frontal  = data.frame(x= ps_grid[[dd+1]][,1] , y =ps_grid[[dd+1]][,2], z = 0 )
    fcoef2.df[[simul]][[dd+1]]$non_frontal  = data.frame(x= ps_grid[[dd+1]][,1] , y =ps_grid[[dd+1]][,2], z = 0 )
    fcoef2.df[[simul]][[dd+1]]$diff  = data.frame(x= ps_grid[[dd+1]][,1] , y =ps_grid[[dd+1]][,2], z = 0 )
    
    
    tmp.env = new.env()
    load(paste0('./GBM_R/gbm_simul/simul1_ps/eigen_',simul,'.RData'), envir = tmp.env)
    gbm_eigen = tmp.env$tmp
    hh.whole = sapply(tmp.env$tmp$propEV, function(x){ min(which(x>0.9)) } )
    rm(tmp.env)
    
    tmp_colname  = apply(cbind( rep(paste0('dim',dd),hh.whole[dd+1]) , 1:hh.whole[dd+1]  ), 1, paste0,collapse = '.')
    ps.area = fcoef.df[[simul]][[dd+1]]$non_frontal$x<=fcoef.df[[simul]][[dd+1]]$non_frontal$y #y>x area.
    
      tmp_colname.f = paste0(tmp_colname,'.f')
      
      coef.nf = as.vector(  gbm_eigen$Eigenfunc[[dd+1]][ ,1:hh.whole[dd+1]]  %*% store_coef_Fcox[[simul]] [tmp_colname,,drop=F]  )
      coef.f = as.vector(  gbm_eigen$Eigenfunc[[dd+1]][ ,1:hh.whole[dd+1]]  %*% store_coef_Fcox[[simul]] [tmp_colname.f,,drop=F]  )
      fcoef.df[[simul]][[dd+1]]$frontal$z[ps.area]    = coef.nf + coef.f 
      fcoef.df[[simul]][[dd+1]]$non_frontal$z[ps.area] = coef.nf
      fcoef.df[[simul]][[dd+1]]$diff$z[ps.area] = coef.f
      
      coef2.nf = as.vector(  gbm_eigen$Eigenfunc[[dd+1]][ ,1:hh.whole[dd+1]]  %*% store_coef_Fcox2[[simul]] [tmp_colname,,drop=F]  )
      coef2.f = as.vector(  gbm_eigen$Eigenfunc[[dd+1]][ ,1:hh.whole[dd+1]]  %*% store_coef_Fcox2[[simul]] [tmp_colname.f,,drop=F]  )
      fcoef2.df[[simul]][[dd+1]]$frontal$z[ps.area]    = coef2.nf + coef2.f 
      fcoef2.df[[simul]][[dd+1]]$non_frontal$z[ps.area]   = coef2.nf
      fcoef2.df[[simul]][[dd+1]]$diff$z[ps.area]  = coef2.f
  }
}

#### Averaged functional coef.
fcoef.ave  = vector(mode = 'list',length=2) #dim=0,1
fcoef.ave.df  = vector(mode = 'list',length=2) 
fcoef2.ave  = vector(mode = 'list',length=2) #dim=0,1
fcoef2.ave.df  = vector(mode = 'list',length=2) 


for( dd in 0:1){  # topological dim
  fcoef.ave[[dd+1]]  = vector(mode = 'list',length=3) # c('frontal', 'non_frontal', 'diff')
  fcoef.ave.df[[dd+1]]  = vector(mode = 'list',length=3)
  names(fcoef.ave[[dd+1]]) = names(fcoef.ave.df[[dd+1]]) = c('frontal', 'non_frontal', 'diff')
  
  fcoef2.ave[[dd+1]]  = vector(mode = 'list',length=3) # c('frontal', 'non_frontal', 'diff')
  fcoef2.ave.df[[dd+1]]  = vector(mode = 'list',length=3)
  names(fcoef2.ave[[dd+1]]) = names(fcoef2.ave.df[[dd+1]]) = c('frontal', 'non_frontal', 'diff')
  
  fcoef.ave.df[[dd+1]]$frontal  = Reduce('+', lapply(fcoef.df, function(x) x[[dd+1]]$frontal))/simul_size
  fcoef.ave.df[[dd+1]]$non_frontal  = Reduce('+', lapply(fcoef.df, function(x) x[[dd+1]]$non_frontal))/simul_size
  fcoef.ave.df[[dd+1]]$diff  = Reduce('+', lapply(fcoef.df, function(x) x[[dd+1]]$diff))/simul_size
  
  fcoef2.ave.df[[dd+1]]$frontal  = Reduce('+', lapply(fcoef2.df, function(x) x[[dd+1]]$frontal))/simul_size
  fcoef2.ave.df[[dd+1]]$non_frontal  = Reduce('+', lapply(fcoef2.df, function(x) x[[dd+1]]$non_frontal))/simul_size
  fcoef2.ave.df[[dd+1]]$diff  = Reduce('+', lapply(fcoef2.df, function(x) x[[dd+1]]$diff))/simul_size
}




#z.limit  = max(sapply(fcoef.ave.df,function(x){ max( c(max(abs(x$frontal$z)),max(abs(x$non_frontal$z)),max(abs(x$diff$z)))  )}  ) )
for( dd in 0:1){  # topological dim
  z.limit  = max(sapply(fcoef.ave.df[[dd+1]],function(x){ max(abs(x$z)) }  ) )
  
  #z.limit = max(abs(fcoef.ave.df[[dd+1]]$frontal$z) )
  fcoef.ave[[paste0('dim',dd)]]$frontal <- ggplot(fcoef.ave.df[[dd+1]]$frontal,aes(x,y,fill=z))+
    geom_raster()+
    theme_tufte()+
    #labs(x="birth",y="death",title = paste0('Entire surface of dim=',dd)) +
    labs(x="birth",y="death",title = NULL) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(slope=1,intercept=0)+
    scale_fill_gradient2(limits=c( -z.limit , z.limit ) ) +
    theme(legend.position = c(0.9, 0.25))+
    theme(plot.title = element_text(hjust = 0.4))
  
  #z.limit = max(abs(fcoef.ave.df[[dd+1]]$non_frontal$z) )
  fcoef.ave[[paste0('dim',dd)]]$non_frontal<- ggplot(fcoef.ave.df[[dd+1]]$non_frontal,aes(x,y,fill=z))+
    geom_raster()+
    theme_tufte()+
    #labs(x="birth",y="death",title = paste0('Entire surface of dim=',dd)) +
    labs(x="birth",y="death",title = NULL) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(slope=1,intercept=0)+
    scale_fill_gradient2(limits=c( -z.limit , z.limit ) ) +
    theme(legend.position = c(0.9, 0.25))+
    theme(plot.title = element_text(hjust = 0.4))
  
  #z.limit = max(abs(fcoef.ave.df[[dd+1]]$diff$z) )
  fcoef.ave[[paste0('dim',dd)]]$diff<- ggplot(fcoef.ave.df[[dd+1]]$diff,aes(x,y,fill=z))+
    geom_raster()+
    theme_tufte()+
    #labs(x="birth",y="death",title = paste0('Entire surface of dim=',dd)) +
    labs(x="birth",y="death",title = NULL) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(slope=1,intercept=0)+
    scale_fill_gradient2(limits=c( -z.limit , z.limit ) ) +
    theme(legend.position = c(0.9, 0.25))+
    theme(plot.title = element_text(hjust = 0.4))
}


#z.limit  = max(sapply(fcoef2.ave.df,function(x){ max( c(max(abs(x$frontal$z)),max(abs(x$non_frontal$z)),max(abs(x$diff$z)))  )}  ) )
for( dd in 0:1){  # topological dim
  z.limit  = max(sapply(fcoef2.ave.df[[dd+1]],function(x){ max(abs(x$z)) }  ) )
  
  #z.limit = max(abs(fcoef2.ave.df[[dd+1]]$frontal$z) )
  fcoef2.ave[[paste0('dim',dd)]]$frontal <- ggplot(fcoef2.ave.df[[dd+1]]$frontal,aes(x,y,fill=z))+
    geom_raster()+
    theme_tufte()+
    #labs(x="birth",y="death",title = paste0('Entire surface of dim=',dd)) +
    labs(x="birth",y="death",title = NULL) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(slope=1,intercept=0)+
    scale_fill_gradient2(limits=c( -z.limit , z.limit ) ) +
    theme(legend.position = c(0.9, 0.25))+
    theme(plot.title = element_text(hjust = 0.4))
  
  #z.limit = max(abs(fcoef2.ave.df[[dd+1]]$non_frontal$z) )
  fcoef2.ave[[paste0('dim',dd)]]$non_frontal<- ggplot(fcoef2.ave.df[[dd+1]]$non_frontal,aes(x,y,fill=z))+
    geom_raster()+
    theme_tufte()+
    #labs(x="birth",y="death",title = paste0('Entire surface of dim=',dd)) +
    labs(x="birth",y="death",title = NULL) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(slope=1,intercept=0)+
    scale_fill_gradient2(limits=c( -z.limit , z.limit ) ) +
    theme(legend.position = c(0.9, 0.25))+
    theme(plot.title = element_text(hjust = 0.4))
  
  #z.limit = max(abs(fcoef2.ave.df[[dd+1]]$diff$z) )
  fcoef2.ave[[paste0('dim',dd)]]$diff<- ggplot(fcoef2.ave.df[[dd+1]]$diff,aes(x,y,fill=z))+
    geom_raster()+
    theme_tufte()+
    #labs(x="birth",y="death",title = paste0('Entire surface of dim=',dd)) +
    labs(x="birth",y="death",title = NULL) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(slope=1,intercept=0)+
    scale_fill_gradient2(limits=c( -z.limit , z.limit ) ) +
    theme(legend.position = c(0.9, 0.25))+
    theme(plot.title = element_text(hjust = 0.4))
}


save(fcoef.df, fcoef.ave.df, fcoef2.df, fcoef2.ave.df,
     file=paste0('./GBM_R/gbm_simul/gbm_simul1_fcoef_df.RData'))
save(fcoef.ave, fcoef2.ave,
     file=paste0('./GBM_R/gbm_simul/gbm_simul1_fcoef.RData'))