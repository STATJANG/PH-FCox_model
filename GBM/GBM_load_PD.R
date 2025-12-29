rm(list=ls())

load_pd_locs = function(gbm_path){
  gbm_locs_files = list.files(gbm_path, pattern = "_locs\\.txt$")
  gbm_pd_files = list.files(gbm_path, pattern = "_pd\\.txt$")
  
  #list[patient][dim][death or birth] = cbind(x,y,z)
  locs_name = unname(sapply(gbm_locs_files , function(s) strsplit(s,split = "_")[[1]][1]))
  pd_name = unname(sapply(gbm_pd_files, function(s) strsplit(s,split = "_")[[1]][1]))
  if(all(locs_name==pd_name)==FALSE) {stop("The orders of patients is not matched.")}
  
  #The size of MNI152 coordinate system.
  x_size = 182
  y_size = 218 
  z_size = 182
  
  gbm_pd_locs_list = vector(mode='list', length=length(gbm_locs_files))
  
  for(patient in 1:length(gbm_pd_locs_list)){
    locs_df = read.table(paste(gbm_path,"/" ,gbm_locs_files[patient], sep = "") , col.names = c("dim","birth","death"))
    pd_df   = read.table(paste(gbm_path,"/" ,gbm_pd_files[patient]  , sep = "") , col.names = c("dim","birth","death"))
    
    #Birth Cartesian coordinate
    x_coor_b = locs_df$birth%%x_size
    y_coor_b = (locs_df$birth%%(x_size*y_size))%/%x_size
    z_coor_b = locs_df$birth%/%(x_size*y_size)
    #Death Cartesian coordinate.
    x_coor_d = locs_df$death%%x_size
    y_coor_d = (locs_df$death%%(x_size*y_size))%/%x_size
    z_coor_d = locs_df$death%/%(x_size*y_size)
    ## 0th dim => birth, (1st, 2nd) dim=> death
    # orientation = c(90,126,72)
    tmp.coord = cbind(90-x_coor_b,y_coor_b-126,z_coor_b-72, 90-x_coor_d,y_coor_d-126,z_coor_d-72)
    colnames(tmp.coord) = c('x.birth','y.birth','z.birth','x.death','y.death','z.death')
    
    
    ##Extracting PDs.
    pd_df = pd_df[-min(which(pd_df$dim==0)),]
    pd_df$death[is.infinite(pd_df$death)] = pd_df$birth[is.infinite(pd_df$death)]
    if(dim(pd_df)[1] != dim(locs_df)[1]) { cat('Iter stops at ',patient,'\n');stop("Dims of pd_df and locs_df are different.")  }
    
    colnames(pd_df) = c('dim','pd.birth','pd.death')
    tmp_df = cbind(pd_df, locs_df[,c('birth','death')], tmp.coord)
    colnames(tmp_df) = c('dim','birth','death','lex.birth','lex.death',colnames(tmp.coord))
    gbm_pd_locs_list[[patient]] = tmp_df
  }
  names(gbm_pd_locs_list) = locs_name
  #gbm_ids = sapply(strsplit(list.files(gbm_path_AT_NonAT, pattern = "_locs\\.txt$"),"_"), function(x) x[1] )
  return(gbm_pd_locs_list)
}


#Ternary images with combined labels
gbm_path = "./GBM/Features_GBM" #active vs inactive
gbm_pd   = load_pd_locs(gbm_path)
gbm_pd = gbm_pd[!names(gbm_pd)%in%c("TCGA-06-0128","TCGA-12-1601")] #excluding images of a poor quality.

#Save data as workspace
#save(gbm_pd, file= "./GBM/GBM_pd.RData")  
