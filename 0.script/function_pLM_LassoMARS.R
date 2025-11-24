## function_pLM_LassoMARS.R
library(earth) 
library(glmnet)
library(pROC)
library(Rfast)
library(dplyr)
library(survival)
library(survminer)
library(stringr)

semi_quantitative<-function(zlfc){
  
  # Converts continuous Z-score LFC (log fold change) values into semi-quantitative categorical scores. 
  # This function applies thresholds to cap extreme values (above 4 and below -8) and then rounds 
  # the normalized values to create discrete phenotypic categories. The transformation enables 
  # categorical analysis of continuous drug response data while preserving the relative magnitude 
  # of effects. Useful for converting continuous screening data into discrete classes for 
  # machine learning model training and categorical statistical analysis.
  
  phenotype<-zlfc
  
  phenotype[which(phenotype>4)]<-4
  phenotype[which(phenotype<(-8))]<-(-8)
  phenotype<-round(phenotype/4)
  
  phenotype
}

screen_site_sep<-function(screen_sites){
  
  # Parses and separates concatenated amino acid mutation strings into structured data frame format.
  # This function processes mutation data that may contain multiple mutations separated by pipes 
  # or semicolons, extracting reference amino acid, alternative amino acid, and position information 
  # into separate columns. Handles complex mutation strings from high-throughput screening data 
  # and returns a tidy data frame with individual mutation records for downstream analysis.
  
  sites_sep <- separate_longer_delim(screen_sites,c(aa_mut), delim = "|")
  sites_sep<- separate_longer_delim(sites_sep,c(aa_mut), delim = ";")
  
  aa_mut<-sites_sep$aa_mut
  sites_sep$aa_ref<-substr(aa_mut,1,1)
  sites_sep$aa_alt<-substr(aa_mut,nchar(aa_mut),nchar(aa_mut))
  sites_sep$aa_pos<-substr(aa_mut,2,nchar(aa_mut)-1)
  
  sites_sep$aa_pos<-as.numeric(sites_sep$aa_pos)
  
  sites_sep[which(sites_sep$aa_pos!=""),]
  
  
}


sample_split<-function(all_info){
  
  # Implements randomized train-validation split for machine learning model development.
  # This function randomly partitions the dataset into training (90%) and validation (10%) sets 
  # while maintaining reproducibility through fixed random seed. Returns a vector indicating 
  # the split assignment for each sample, enabling proper model training and evaluation protocol.
  
  nsample<-nrow(all_info)  
  split<-rep("",nsample)

  set.seed(234)
  random_index<-sample(c(1:nsample),nsample,replace = F)
  ntrain<-ceiling(nsample*0.9)
  train_index<-random_index[1:ntrain]
  vali_index<-random_index[(ntrain+1):nsample]
  
  split[train_index]<-"train"
  split[vali_index]<-"validation"

  
  split
}


convert_feature<-function(prot_eb,info){
  
  # Integrates protein embedding features with structural and biophysical properties for model input.
  # This function combines deep learning-based protein embeddings with manually curated structural 
  # features (CF10RK, LD15RK) and amino acid properties (mass, isoelectric point, hydrophobicity)
  # to create comprehensive feature matrices for machine learning models. Ensures proper feature 
  # alignment and formatting for predictive modeling tasks.

  prot_input<-cbind(prot_eb[,-1],info[,c("CF10RK","CF10","LD15RK","LD15","nAngels.ref","Mass.ref","IP.ref","HM.ref","nAngels.alt","Mass.alt","IP.alt","HM.alt")])
  
  prot_input
}

model_train<-function(prot_eb,sample_info,label){
  
  # Trains ensemble machine learning models using LASSO feature selection and MARS regression.
  # This function implements a two-stage modeling approach: first applies LASSO regularization 
  # to select informative features, then trains Multivariate Adaptive Regression Splines (MARS)
  # on the selected features. Performs cross-validation, computes performance metrics (correlation 
  # and AUC), and returns the trained model with validation results for drug response prediction.
  
  train_index<-which(sample_info$split=="train")
  vali_index<-which(sample_info$split=="validation")
  
  train_info<-sample_info[train_index,]
  train_prot_eb<-prot_eb[train_index,]
  
  vali_info<-sample_info[vali_index,]
  vali_prot_eb<-prot_eb[vali_index,]
  
  print(identical(vali_info$id,vali_prot_eb[,1]))
  
  ####
  lasso<-coef(glmnet(x=train_prot_eb[,-1],y=train_info[,"phenotype"]), s = 0.001)
  select<-lasso@i+1
  print(length(select))
  train_prot_eb<-train_prot_eb[,select]
  vali_prot_eb<-vali_prot_eb[,select]
  ########
  
  prot_eb_input<-convert_feature(train_prot_eb,train_info)
  prot_eb_vali<-convert_feature(vali_prot_eb,vali_info)
  
  
  prot_eb_input$response<-train_info[,"phenotype"]
  
  prot_mars <- earth(response ~ ., data = prot_eb_input,degree=2)
  # summary(prot_mars)
  
  
  prot_model_info<-model_validation(model=prot_mars,vali_input=prot_eb_vali,vali_obs=vali_info)
  prot_model_cor<-prot_model_info[[1]]
  prot_model_roc<-prot_model_info[[2]]
  
  perf<-c(type=label,
          prot_model_cor=round(prot_model_cor,3),
          prot_model_auc=round(prot_model_roc$auc,3))
  
  list(model=prot_mars,select=select,validation=perf)
  
  
}

model_validation<-function(model,vali_input,vali_obs){
  
  # Evaluates trained model performance on validation data using correlation and ROC analysis.
  # This function generates predictions on held-out validation data, computes Pearson correlation 
  # between predicted and observed values, and calculates ROC curves with AUC metrics for 
  # categorical classification performance. Provides comprehensive model assessment for 
  # regression and classification tasks.
  
  fit_vali<-as.numeric(predict(model,vali_input))
  vali_cor<-cor(vali_obs[,"phenotype"],fit_vali)
  
  vali_obs<-vali_site_proc(vali_obs)
  roc <- roc(response = vali_obs$category,predictor = fit_vali,
             direction=">")
  
  list(vali_cor,roc)
  
  
}

vali_site_proc<-function(df){
  
  # Converts continuous phenotypic scores into binary categorical labels for classification tasks.
  # This function transforms quantitative phenotype measurements into discrete "sensitive" 
  # and "non-response/resistant" categories based on predefined thresholds. Enables binary 
  # classification analysis and ROC evaluation of model performance for drug response prediction.
  
  df$category<-"-"
  df$category[which(df[,"phenotype"]<0)]<-"sensitive"
  df$category[which(df[,"phenotype"]>=0)]<-"non-response/resistent"
  
  df
  
}

vali_pred<-function(model,all_eb,all_info){
  
  # Generates predictions using trained models on validation dataset.
  # This function applies pre-trained models to validation data, converts features to appropriate 
  # format, and returns predictions alongside original observation data. Facilitates model 
  # performance evaluation and comparison across different modeling approaches.
  
  index<-which(all_info$split=="validation")
  
  info<-all_info[index,]
  input<-convert_feature(all_eb[index,],info)
  
  vali<-predict(model,input)
  
  info$prediction<-as.numeric(vali)
  
  info
  
}

eb_CF_to_model_input<-function(CF_df,prot,prot_eb,AA_eb,output_prefix){
  
  # Prepares protein and amino acid embedding data for machine learning model input.
  # This function integrates protein-level and amino acid-level embeddings from protT5
  # models, ensures proper alignment between mutation data and embedding features, and exports 
  # formatted input files for model training and prediction. Handles large-scale embedding data 
  # from protein language models and structural predictions.
  
  CF_in<-CF_df[which(CF_df$uniprot_id%in%prot),]
  # print("sample:")
  # print(nrow(CF_in))
  
  write.table(CF_in,
              paste0(output_prefix,"_test_info.txt"),
              append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
  
  ################# c1  ##############
  colnames(prot_eb)<-c("prot",paste0("p_",c(0:1023)))
  use_prot<-prot_eb[which(prot_eb$prot%in%prot),]
  rownames(use_prot)<-use_prot$prot
  eb_prot<-use_prot[CF_in$uniprot_id,]
  print("C1:")
  print(identical(eb_prot[,1],CF_in$uniprot_id))
  print(which(is.na(eb_prot)))
  
  write.table(eb_prot,
              paste0(output_prefix,"_eb_prot_c1.txt"),
              append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
  ################# c2  ##############
  colnames(AA_eb)<-c("prot.aa",paste0("aa_",c(0:1023)))
  use_AA<-AA_eb[which(AA_eb$prot.aa%in%CF_in$id),]
  rownames(use_AA)<-use_AA$prot.aa
  
  feature_AA<-use_AA[CF_in$id,]
  # print(head(feature_AA[,1:20]))
  # print(head(CF_in))
  
  print(identical(feature_AA[,1],CF_in$id))
  print(which(is.na(feature_AA)))
  
  write.table(feature_AA,
              paste0(output_prefix,"_eb_aa_c2.txt"),
              append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
  ###
  print(identical(str_split_fixed(feature_AA[,1],"\\.",n=2)[,1],eb_prot[,1]))
  
  all_eb<-cbind(feature_AA,eb_prot[,-1])
  
  write.table(all_eb,
              paste0(output_prefix,"_eb_all.txt"),
              append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
}


plot_drug_box<-function(drug_phenotype,para,test_folder,drug_scope,drug_list,ofolder){
  
  # Generates multi-panel boxplot visualizations for drug prediction results across multiple parameters.
  # This function creates comprehensive comparison plots showing model predictions for different 
  # drugs and feature sets. Organizes results in grid layout for easy comparison of prediction 
  # performance across various experimental conditions and modeling approaches.
  
  plot_list<-list()
  m<-1
  for(plot_var in para){

    plot_list[[m]]<-drug_prediction(drug_list[1],plot_var,test_folder,drug_scope)
    plot_list[[m+1]]<-drug_prediction(drug_list[2],plot_var,test_folder,drug_scope)
    plot_list[[m+2]]<-drug_prediction(drug_list[3],plot_var,test_folder,drug_scope)
    plot_list[[m+3]]<-drug_prediction(drug_list[4],plot_var,test_folder,drug_scope)
    
    m<-m+4
    
  }
  
  p<-ggarrange(plotlist=plot_list,ncol=4,nrow=length(para),common.legend =T)
  
  ggsave(paste0(ofolder,"Drug_",drug_scope,"_Boxplot.pdf"),p,width=9,height=3.5*length(para))
  

}


boxplot_drug_prediction<-function(plot_data,use_feature,use_label){
  
  # Creates comparative boxplots with statistical testing for drug response predictions.
  # This function visualizes model predictions across different response categories (resistant 
  # vs sensitive) with embedded statistical significance testing using Wilcoxon rank-sum tests. 
  # Includes individual data points as jittered dots to show distribution density and provides 
  # publication-ready visualization of prediction performance.
  
  my_comparisons <- list(c("Resistance","Sensitive/Ambiguous"))
  
  p<-ggplot(plot_data,aes(x = factor(Classification), y = prediction))+
    labs(title=stringr::str_wrap(use_label,width=40))+ 
    xlab("")+
    ylab(use_feature)+
    geom_boxplot(lwd=0.3,fatten=1,outlier.shape = NA)+
    geom_jitter(aes(color=gene),shape=16, position = position_jitter(0.2),alpha=0.5)+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),comparisons = my_comparisons,method ="wilcox.test",method.args = list(alternative = "greater"))+
    theme_minimal() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p
}

drug_prediction<-function(drug,plot_var,test_folder,suffix){
  
  # Executes end-to-end drug response prediction pipeline for specific compounds.
  # This function loads pre-processed data, applies trained models to generate predictions, 
  # classifies responses into resistance categories, and generates evaluation visualizations. 
  # Integrates multiple modeling approaches and provides comprehensive prediction outputs 
  # for drug development applications.
  
  drug_all_eb<-read.table(paste0(test_folder,"Drug_",drug,"_",suffix,"_eb_all.txt"),quote="",header=T,sep="\t",fill=TRUE)
  drug_info<-read.table(paste0(test_folder,"Drug_",drug,"_",suffix,"_test_info.txt"),quote="",header=T,sep="\t",fill=TRUE)
  
  drug_pred<-stress_evaluate(CDKi_model,MEKi_model,PARPi_model,drug_all_eb,drug_info)
  
  drug_pheno_pred<-inner_join(drug_info,drug_pred,by=join_by("id","uniprot_id","aa_pos","aa_ref","aa_alt","gene"))
  
  drug_pheno_pred$Classification<-"Sensitive/Ambiguous"
  drug_pheno_pred$Classification[which(drug_pheno_pred$resistance==1)]<-"Resistance"
  print(plot_var)
  #print(drug_pheno_pred[1,])
  drug_pheno_pred$prediction<-drug_pheno_pred[,plot_var]
  
  write.table(drug_pheno_pred,
              paste0("Drug_",drug,"_",suffix,"_prediction.txt"),
              append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)
  
  boxplot<-boxplot_drug_prediction(drug_pheno_pred,plot_var,drug)
  
  boxplot
}

stress_evaluate<-function(abem_model,bini_model,olap_model,core_prot_eb,core_info){
  
  # Applies multiple trained models to generate ensemble predictions for stress response.
  # This function leverages predictions from different drug response models (Abemaciclib, 
  # Binimetinib) and computes aggregated stress response scores. Enables comprehensive 
  # assessment of cellular stress phenotypes using multi-model consensus approaches.

  
  abem_pred<-test_pred(model=abem_model[["model"]],prot_eb=core_prot_eb[,abem_model[["select"]]],info=core_info)
  
  bini_pred<-test_pred(model=bini_model[["model"]],prot_eb=core_prot_eb[,bini_model[["select"]]],info=core_info)
  
  olap_pred<-test_pred(model=olap_model[["model"]],prot_eb=core_prot_eb[,olap_model[["select"]]],info=core_info)
  
  pred<-cbind(core_info[,c(1:6)],abem_pred=abem_pred,bini_pred=bini_pred,olap_pred=olap_pred)
  
  pred$StressResponse<-rowMeans(as.matrix(pred[,c("abem_pred","bini_pred")]))
  
  pred
  
}


test_pred<-function(model,prot_eb,info){
  
  # Generates predictions using single trained model on test dataset.
  # This function applies individual pre-trained models to new data, handles feature conversion, 
  # and returns numerical predictions. Serves as core prediction engine for model evaluation 
  # and application to new experimental data.
  
  test_input<-convert_feature(prot_eb=prot_eb,info=info)
  
  as.numeric(predict(model,test_input))
  
}

patient_predict_by_site<-function(site_pred,resisent_cut,obs){
  
  # Aggregates site-level predictions to patient-level resistance scores and clinical outcomes.
  # This function converts individual mutation predictions into patient-level metrics by 
  # calculating average stress response and counting resistant sites. Integrates with clinical 
  # data (overall survival) to enable survival analysis and clinical correlation studies.
  
  site_pred$Resistent<-0
  site_pred$Resistent[which(site_pred[,"StressResponse"]>resisent_cut)]<-1
  print(table(site_pred$Resistent))
  
  patient_pred<-site_pred %>%
    group_by(PATIENT_ID) %>%
    summarise(
      StressResponse = mean(StressResponse),
      Resistent_count = sum(Resistent),
      nSite = length(StressResponse)
    )
  
  patient<-as.data.frame(inner_join(patient_pred,obs,by=join_by("PATIENT_ID")))
  patient<-patient[which(patient$OS_MONTHS!=""),]
  patient<-patient[which(patient$OS_STATUS!=""),]
  
  patient
}


plot_OS_surv<-function(data,var,title,pdf_flag){
  
  # Generates Kaplan-Meier survival curves based on continuous prediction variables.
  # This function creates survival plots comparing high vs low prediction groups using median 
  # cutoff, includes confidence intervals and log-rank p-values. Produces publication-ready 
  # survival analysis visualizations for clinical outcome prediction validation.
  
  low_cut<-quantile(data[,var],0.5)
  high_cut<-quantile(data[,var],0.5)
  
  print(low_cut)
  print(high_cut)
  
  data$pred_type<-""
  data$pred_type[which(data[,var]<low_cut)]<-paste0(var,"  low (< Median)")
  data$pred_type[which(data[,var]>=high_cut)]<-paste0(var," high (>= Median)")

  
  print(table(data$pred_type))
  
  plot_data<-data[which(data$pred_type!=""),]
  
  plot_data$OS_STATUS<-as.numeric(str_split_fixed(plot_data$OS_STATUS,":",n=2)[,1])
  
  fit <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ pred_type, data = plot_data)
  
  plot<-ggsurvplot(fit, data = plot_data,conf.int = TRUE, pval = TRUE,title=title) + theme_survminer(base_size = 0.5)
  
  if(pdf_flag){
    
    ggsave(paste0("Surv_quant_",title,".pdf"), print(plot$plot), width = 6, height = 5)
    # 
    # surv_pvalue(fit,plot_data)
  }
  

  plot
  
}

plot_OS_surv_by_cut<-function(data,var,cut,pdf_prefix,pdf_flag){
  
  # Creates survival plots using predefined cutoff values for binary group classification.
  # This function enables flexible survival analysis with user-specified thresholds for 
  # continuous prediction variables. Useful for testing specific clinical decision boundaries 
  # and biomarker cutoff optimization.
  
  data$pred_type<-""
  data$pred_type[which(data[,var]<cut)]<-paste0(var,"<",cut)
  
  data$pred_type[which(data[,var]>=cut)]<-paste0(var,">=",cut)
  
  print(table(data$pred_type))
  
  plot_data<-data[which(data$pred_type!=""),]
  
  plot_data$OS_STATUS<-as.numeric(str_split_fixed(plot_data$OS_STATUS,":",n=2)[,1])
  
  fit <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ pred_type, data = plot_data)
  
  plot<-ggsurvplot(fit, data = plot_data,conf.int = TRUE, pval = TRUE,title=pdf_prefix)
  
  if(pdf_flag){
    
    ggsave(paste0("Surv_cut_",title,".pdf"), print(plot$plot), width = 6, height = 5)
    # 
    # surv_pvalue(fit,plot_data)
  }
  
  plot
  
}


plot_OS_surv_by_cut2<-function(data,var,high_cut,low_cut,pdf_prefix,pdf_flag){
  
  # Generates survival plots comparing extreme groups using dual cutoff approach.
  # This function visualizes survival differences between high and low extremes of prediction 
  # distributions, excluding intermediate values. Provides enhanced contrast for identifying 
  # strong prognostic factors and includes case counts in group labels for clinical interpretation.
  
  data$pred_type<-""
  nlow<-length(which(data[,var]<low_cut))
  data$pred_type[which(data[,var]<low_cut)]<-paste0(var,"<",low_cut,"(",nlow,")")
  
  nhigh<-length(which(data[,var]>=high_cut))
  data$pred_type[which(data[,var]>=high_cut)]<-paste0(var,">=",high_cut,"(",nhigh,")")
  
  print(table(data$pred_type))
  
  plot_data<-data[which(data$pred_type!=""),]
  
  plot_data$OS_STATUS<-as.numeric(str_split_fixed(plot_data$OS_STATUS,":",n=2)[,1])
  
  fit <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ pred_type, data = plot_data)
  
  plot<-ggsurvplot(fit, data = plot_data,conf.int = TRUE, pval = TRUE,palette = c("#00BFC4","#F8766D"),title=pdf_prefix)
  library(scales)
  hue_pal()(2)
  
  if(pdf_flag){
    
    ggsave(paste0("Surv_cut_",title,".pdf"), print(plot$plot), width = 6, height = 5)
    # 
    # surv_pvalue(fit,plot_data)
  }
  
  plot
  
}

plot_OS_surv_by_qcut<-function(data,var,high_cut,low_cut,pdf_prefix,pdf_flag){
  
  # Creates survival plots using quantile-based cutoffs for group definition.
  # This function employs percentile-based thresholds (e.g., top and bottom quartiles) to 
  # define comparison groups, making the analysis robust to distribution skewness. Includes 
  # group sizes in labels and provides flexible quantile-based survival analysis.
  
  print(pdf_prefix)
  data$pred_type<-""
  
  low_value<-quantile(data[,var],low_cut)
  print(low_value)
  
  nlow<-length(which(data[,var]<low_value))
  data$pred_type[which(data[,var]<low_value)]<-paste0(var,"<",low_cut*100,"% (",nlow,")")
  
 
  
  high_value<-quantile(data[,var],high_cut)
  print(high_value)
  
  nhigh<-length(which(data[,var]>high_value))
  data$pred_type[which(data[,var]>high_value)]<-paste0(var,">",high_cut*100,"% (",nhigh,")")
  
  print(table(data$pred_type))
  
  plot_data<-data[which(data$pred_type!=""),]
  
  plot_data$OS_STATUS<-as.numeric(str_split_fixed(plot_data$OS_STATUS,":",n=2)[,1])
  
  fit <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ pred_type, data = plot_data)
  
  plot<-ggsurvplot(fit, data = plot_data,conf.int = TRUE, pval = TRUE,title=pdf_prefix)
  
  if(pdf_flag){
    
    ggsave(paste0("Surv_qcut_",title,".pdf"), print(plot$plot), width = 6, height = 5)
    # 
    # surv_pvalue(fit,plot_data)
  }
  
  plot
  
}

resistence_to_Phos_dist<-function(aa_phenotype,Phos_nb){
  
  # Analyzes spatial relationship between resistance mutations and phosphorylation sites.
  # This function identifies resistant mutations and examines their proximity to known 
  # phosphorylation sites using structural distance metrics. Enables investigation of 
  # potential mechanistic links between drug resistance and post-translational modification 
  # landscapes in protein structures.
  
  
  resp_var<-"StressResponse"
  resisent_cut<-0.2
  
  aa_phenotype$id<-paste0(aa_phenotype$aa_pos,aa_phenotype$aa_ref)
  
  aa_phenotype<-aa_phenotype[,c("uniprot_id","id","aa_pos","aa_ref","aa_alt",resp_var,"LD15RK","CF10RK")]
  aa_phenotype<-unique(aa_phenotype)
  
  aa_phenotype<-aa_phenotype[which(aa_phenotype[,resp_var]>resisent_cut),]
  
  Phos_nb<-Phos_nb[which(Phos_nb$uniprot_id%in%aa_phenotype$uniprot_id),]
  
  aa_phenotype_Phos<-left_join(Phos_nb,aa_phenotype,by=join_by("uniprot_id","id"))
  summary(aa_phenotype_Phos[,resp_var])
  
  aa_phenotype_Phos$Classification<-"Ambiguous"
  aa_phenotype_Phos$Classification[which(aa_phenotype_Phos[,resp_var]>resisent_cut)]<-"Resistance"
  
  print(table(aa_phenotype_Phos$Classification))
  
  aa_phenotype_Phos
}




# resistence_to_Phos_dist<-function(aa_phenotype,Phos_nb){
#   
#   resp_var<-"StressResponse"
#   resisent_cut<-0.2
#   
#   aa_phenotype$id<-paste0(aa_phenotype$aa_pos,aa_phenotype$aa_ref)
#   
#   aa_phenotype<-aa_phenotype[,c("uniprot_id","id","aa_pos","aa_ref","aa_alt",resp_var,"LD15RK","CF10RK")]
#   aa_phenotype<-unique(aa_phenotype)
#   
#   aa_phenotype_Phos<-inner_join(Phos_nb,aa_phenotype,by=join_by("uniprot_id","id"))
#   summary(aa_phenotype_Phos[,resp_var])
#   
#   aa_phenotype_Phos$Classification<-"Ambiguous"
#   aa_phenotype_Phos$Classification[which(aa_phenotype_Phos[,resp_var]>resisent_cut)]<-"Resistance"
#   
#   print(table(aa_phenotype_Phos$Classification))
#   
#   aa_phenotype_Phos
# }

plot_var_dist<-function(dist_data,plot_var,plot_title){
  
  # Creates density distribution plots for distance variables across different classification groups.
  # This function visualizes the distribution of spatial distances or other continuous variables 
  # between resistance-associated sites and functional regions. Uses kernel density estimation 
  # to compare distributions across different phenotypic categories.
  
  plot_df<-dist_data
  
  ngene<-length(unique(plot_df$uniprot_id))
  
  plot_df$distance<-plot_df[,plot_var]
  
  print(which(is.na(plot_df$distance)))
  
  p<-ggplot(plot_df,aes(x=distance))+
    labs(title=str_wrap(plot_title,40))+ 
    #xlim(xrange)+
    geom_density(aes(colour=Classification),alpha=0.4,show.legend=FALSE)+
    stat_density(aes(colour=Classification),
                 geom="line",position="identity")+
    theme_minimal()
  
  p
  
}



boxplot_var_dist<-function(dist_data,plot_var,plot_title){
  
  # Generates comparative boxplots with statistical testing for distance distributions.
  # This function provides statistical comparison of distance metrics between resistance 
  # and ambiguous mutation groups using violin plots and boxplots. Includes Wilcoxon test 
  # results and enables quantitative assessment of spatial distribution differences.
  
  plot_df<-dist_data
  
  plot_df<-plot_df[which(plot_df[,plot_var]!="-"),]
  
  ngene<-length(unique(plot_df$uniprot_id))
  
  plot_df$distance<-as.numeric(plot_df[,plot_var])
  
  print(which(is.na(plot_df$distance)))
  
  p<-ggplot(plot_df,aes(x = factor(Classification), y = distance))+
    labs(title=str_wrap(plot_title,40))+ 
    geom_violin()+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),comparisons = list(c("Resistance","Ambiguous")),method ="wilcox.test",method.args = list(alternative = "less"))+
    geom_boxplot(width = 0.1,outlier.shape = NA)+
    theme_minimal()
  
  p
  
}


min_comb_aa_dist<-function(aa_data,var_list){
  
  # Computes minimum distance from multiple amino acid distance metrics.
  # This function processes multiple distance measurements (e.g., to different structural 
  # features) and calculates the minimum distance for each position. Handles missing values 
  # and returns the most relevant proximity metric for functional analysis of mutation sites.
  
  for(a in 1:length(var_list)){
    print(var_list[a])
    aa_data[which(aa_data[,var_list[a]]=="-"),var_list[a]]<-Inf
    aa_data[,var_list[a]]<-as.numeric(aa_data[,var_list[a]])
    
  }
  
  min_dist<-matrixStats::rowMins(as.matrix(aa_data[,var_list]))
  
  min_dist[which(min_dist==500)]<-"-"
  min_dist
}
