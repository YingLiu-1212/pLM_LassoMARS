## _forAnalysis.R
min_comb_aa_dist<-function(aa_data,var_list){
  for(a in 1:length(var_list)){
    print(var_list[a])
    aa_data[which(aa_data[,var_list[a]]=="-"),var_list[a]]<-Inf
    aa_data[,var_list[a]]<-as.numeric(aa_data[,var_list[a]])
    
  }
  
  min_dist<-matrixStats::rowMins(as.matrix(aa_data[,var_list]))
  
  min_dist[which(min_dist==500)]<-"-"
  min_dist
}

norm_aa_dist<-function(aa_data,norm_var){
  
  aa_data[which(aa_data[,norm_var]=="-"),norm_var]<-Inf
  
  aa_data[,norm_var]<-as.numeric(aa_data[,norm_var])
  
  prot_id<-unique(aa_data$uniprot_id)
  
  norm_value<-rep(NA,nrow(aa_data))
  
  for(p in 1:length(prot_id)){
    
    aaindex<-which(aa_data$uniprot_id==prot_id[p])
    
    sub_df<-aa_data[aaindex,]
    
    norm_value[aaindex]<-round(scale(sub_df[,norm_var]),4)
    
  }
  
  norm_value
  
}

plot_var_LDcut<-function(dist_data,plot_var,cut){
  
  
  plot_df<-dist_data[which(dist_data$LD15RK<=cut),]
  plot_df$distance<-plot_df[,plot_var]

  which(is.na(plot_df$distance))
  
  p<-ggplot(plot_df,aes(x=distance))+
    labs(title=str_wrap(paste0("LD15RK<=",cut," ",nrow(plot_df)," AA ",plot_var),40))+ 
    #xlim(xrange)+
    geom_density(aes(colour=Classification),alpha=0.4,show.legend=FALSE)+
    stat_density(aes(colour=Classification),
                 geom="line",position="identity")+
    theme_minimal()
  
  p
  
}

plot_var_CFcut<-function(dist_data,plot_var,cut){
  
  plot_df<-dist_data[which(dist_data$CF10RK>=cut),]

  
  plot_df$distance<-plot_df[,plot_var]
  #print(which(is.na(plot_df$distance)))
  
  p<-ggplot(plot_df,aes(x=distance))+
    labs(title=str_wrap(paste0("CF10RK>=",cut," ",nrow(plot_df)," AA ",plot_var),40))+ 
    #xlim(xrange)+
    geom_density(aes(colour=Classification),alpha=0.4,show.legend=FALSE)+
    stat_density(aes(colour=Classification),
                 geom="line",position="identity")+
    theme_minimal()
  
}

dist_aa_site<-function(inter_contact,aa_site){

  if(nrow(aa_site)>0){
    
    dist<-inter_contact[,aa_site$id]
    PTMDist<-matrixStats::rowMins(as.matrix(dist))
    if(nrow(aa_site)>1){
      index<-apply(dist,1,which.min)
      PTMSite<-aa_site$label[index]
    }else{
      PTMSite<-aa_site$label
    }
    
  }else{
    PTMDist<-"-"
    PTMSite<-"-"
  }
  
  data.frame(AADist=PTMDist,AASite=PTMSite)
}





