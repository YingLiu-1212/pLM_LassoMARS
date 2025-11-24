# This R script performs comprehensive validation of machine learning-based stress response predictions on MSK clinical cancer data. The code applies pre-trained Lasso-MARS ensemble models to generate cellular stress response scores for MSK patient mutations and analyzes their association with overall survival outcomes. Key analytical components include: generating stress response predictions for individual mutation sites using ensemble modeling of CDK, MEK, and PARP inhibitor responses; aggregating site-level predictions to patient-level metrics including stress response scores and resistance mutation counts; conducting stratified survival analysis by KRAS G12 mutation status to examine mutation-specific effects; performing cancer type-specific subgroup analyses for tumor types with sufficient sample sizes; and generating multi-panel Kaplan-Meier survival plots for comprehensive visualization. The analysis provides independent clinical validation of computational stress response signatures in metastatic cancer patients, enabling assessment of the prognostic value of protein structure-informed predictions across different molecular subtypes and cancer types in a real-world clinical cohort.

library("RColorBrewer")
library(stringr)
options(scipen=999)
options(stringsAsFactors = FALSE)

source("~/Other_project/CRC/reproducibility/script/function_pLM_LassoMARS.R")
setwd("~/Other_project/CRC/reproducibility/e.validation_clinical")

data_folder<-"../a.data/"
model_folder<-"../c.training/"
output_folder<-"../o.output_figures/"

load(paste0(model_folder,"pLM_LassoMARS.RData"))

msk_all_eb<-read.table(paste0("MSK_Core_eb_all.txt"),quote="",header=T,sep="\t",fill=TRUE)
msk_info<-read.table(paste0("MSK_Core_test_info.txt"),quote="",header=T,sep="\t",fill=TRUE)
msk_obs<-read.table(paste0("MSK_Core_patient_info.txt"),quote="",header=T,sep="\t",fill=TRUE)

print(identical(msk_info$id,msk_all_eb[,1]))

msk_obs<-msk_obs[which(!is.na(msk_obs$OS_MONTHS) & !is.na(msk_obs$OS_STATUS)),]

msk_kras_plus<-msk_obs[which(msk_obs$KRAS12G=="Yes"),"PATIENT_ID"]
msk_kras_index<-msk_info$PATIENT_ID%in%msk_kras_plus

####

msk_pred<-stress_evaluate(CDKi_model,MEKi_model,PARPi_model,msk_all_eb,msk_info)

msk_prediction_info<-left_join(msk_pred,msk_info[,1:13],by=join_by("id","uniprot_id","aa_pos","aa_ref","aa_alt","PATIENT_ID"))


write.table(msk_prediction_info,
            paste0("msk_core73_prediction_info.txt"),
            append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)


msk_cut<-0.2

msk_patient<-patient_predict_by_site(msk_pred,msk_cut,msk_obs)
#table(msk_patient$Resistent_count)

write.table(msk_patient,
            paste0("msk_patient_stress_prediction.txt"),
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)


msk_patient_nonKras<-patient_predict_by_site(stress_evaluate(CDKi_model,MEKi_model,PARPi_model,msk_all_eb[which(!msk_kras_index),],msk_info[which(!msk_kras_index),]),msk_cut,msk_obs)


# table(msk_patient_nonKras$Resistent_count)

msk_patient_Kras<-patient_predict_by_site(stress_evaluate(CDKi_model,MEKi_model,PARPi_model,msk_all_eb[which(msk_kras_index),],msk_info[which(msk_kras_index),]),msk_cut,msk_obs)

# table(msk_patient_Kras$Resistent_count)

######


para<-c("StressResponse")

for(i in 1:length(para)){
  
  plot_var<-para[i]
  
  p4<-plot_OS_surv(data=msk_patient,var=plot_var,paste0("MSK (",nrow(msk_patient)," cases)"),pdf_flag=FALSE)
  
  p5<-plot_OS_surv(data=msk_patient_nonKras,var=plot_var,paste0("MSK KRAS12G(WT) (",nrow(msk_patient_nonKras)," cases)"),pdf_flag=FALSE)
  
  p6<-plot_OS_surv(data=msk_patient_Kras,var=plot_var,paste0("MSK KRAS12G(MUT) (",nrow(msk_patient_Kras)," cases)"),pdf_flag=FALSE)
  
  comb<-arrange_ggsurvplots(list(p4,p5,p6),ncol=3,nrow=1,print=FALSE)
  
  ggsave(paste0(output_folder,"Surv_OS_",plot_var,"_msk.pdf"),comb,limitsize = FALSE,width=12,height=4)
}

#####


table(msk_patient$CANCER_TYPE)
type<-names(which(table(msk_patient$CANCER_TYPE)>200))
for(i in 1:length(para)){
  
  plot_var<-para[i]
  
  p<-list()
  for(j in 1:length(type)){
    
    plot_data<-msk_patient[which(msk_patient$CANCER_TYPE==type[j]),]
    
    mut_type<-unique(msk_info[which(msk_info$PATIENT_ID%in%plot_data$PATIENT_ID),"id"])
    
    print(type[j])
    print(length(mut_type))
    
    p[[j]]<-plot_OS_surv(data=plot_data,var=plot_var,paste0("MSK_",type[j]),pdf_flag=FALSE)
  }
  
  comb<-arrange_ggsurvplots(p,ncol=3,nrow=4,print=FALSE)
  
  ggsave(paste0(output_folder,"Surv_OS_",plot_var,"_msk_type.pdf"),comb,limitsize = FALSE,width=18,height=24)
  
  
}



