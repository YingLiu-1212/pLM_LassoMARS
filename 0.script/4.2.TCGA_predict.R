# This R script performs comprehensive survival analysis using machine learning predictions of cellular stress response in TCGA cancer patients. The code applies pre-trained Lasso-MARS ensemble models (CDKi, MEKi, PARPi) to predict stress response scores from TCGA mutation data, then analyzes the relationship between these predictions and clinical outcomes. Key analytical steps include: generating stress response predictions for individual mutation sites using ensemble modeling; aggregating site-level predictions to patient-level scores and resistance counts; stratifying patients by KRAS G12 mutation status to examine mutation-specific effects; performing Kaplan-Meier survival analysis to compare overall survival between high and low stress response groups; conducting subgroup analyses across different cancer types with sufficient sample sizes; and generating multi-panel survival plots for comprehensive visualization. The analysis enables systematic investigation of how predicted cellular stress responses, derived from protein structural and mutational features, correlate with patient survival outcomes in pan-cancer cohorts, providing insights into the prognostic value of computational stress response signatures in clinical oncology.

library(dplyr)
options(scipen=999)
options(stringsAsFactors = FALSE)
source("~/Other_project/CRC/reproducibility/script/function_pLM_LassoMARS.R")

setwd("~/Other_project/CRC/reproducibility/e.validation_clinical")

data_folder<-"../a.data/"
model_folder<-"../c.training/"
output_folder<-"../o.output_figures/"

load(paste0(model_folder,"pLM_LassoMARS.RData"))

tcga_all_eb<-read.table(paste0("TCGA_Core_all_eb_all.txt"),quote="",header=T,sep="\t",fill=TRUE)

tcga_info<-read.table(paste0("TCGA_Core_all_test_info.txt"),quote="",header=T,sep="\t",fill=TRUE)
tcga_obs<-read.table(paste0("TCGA_Core_all_patient_info.txt"),quote="",header=T,sep="\t",fill=TRUE)

tcga_kras_plus<-tcga_obs[which(tcga_obs$KRAS12G=="Yes"),"PATIENT_ID"]
tcga_kras_index<-tcga_info$PATIENT_ID%in%tcga_kras_plus

####

tcga_pred<-stress_evaluate(CDKi_model,MEKi_model,PARPi_model,tcga_all_eb,tcga_info)

tcga_prediction_info<-left_join(tcga_pred,tcga_info[,1:13],by=join_by("id","uniprot_id","aa_pos","aa_ref","aa_alt","PATIENT_ID"))


write.table(tcga_prediction_info,
            paste0("tcga_core73_prediction_info.txt"),
            append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

########

resp_var<-"StressResponse"

tcga_cut<-0.2

tcga_patient<-patient_predict_by_site(tcga_pred,tcga_cut,tcga_obs)
table(tcga_patient$Resistent_count)


write.table(tcga_patient,
            paste0("tcga_patient_stress_prediction.txt"),
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)

tcga_patient_nonKras<-patient_predict_by_site(stress_evaluate(CDKi_model,MEKi_model,PARPi_model,tcga_all_eb[which(!tcga_kras_index),],tcga_info[which(!tcga_kras_index),]),tcga_cut,tcga_obs)

tcga_patient_Kras<-patient_predict_by_site(stress_evaluate(CDKi_model,MEKi_model,PARPi_model,tcga_all_eb[which(tcga_kras_index),],tcga_info[which(tcga_kras_index),]),tcga_cut,tcga_obs)


para<-c("StressResponse")

for(i in 1:length(para)){
  
  plot_var<-para[i]
  
  p1<-plot_OS_surv(data=tcga_patient,var=plot_var,paste0("TCGA (",nrow(tcga_patient)," cases)"),pdf_flag=FALSE)
  
  p2<-plot_OS_surv(data=tcga_patient_nonKras,var=plot_var,paste0("TCGA KRAS12G(WT) (",nrow(tcga_patient_nonKras)," cases)"),pdf_flag=FALSE)
  
  p3<-plot_OS_surv(data=tcga_patient_Kras,var=plot_var,paste0("TCGA KRAS12G(MUT)(",nrow(tcga_patient_Kras)," cases)"),pdf_flag=FALSE)
  
  comb<-arrange_ggsurvplots(list(p1,p2,p3),ncol=3,nrow=1,print=FALSE)
  
  ggsave(paste0(output_folder,"Surv_OS_",plot_var,"_tcga.pdf"),comb,limitsize = FALSE,width=12,height=4)
}



table(tcga_patient$CANCER_TYPE_ACRONYM)
type<-names(which(table(tcga_patient$CANCER_TYPE_ACRONYM)>200))
for(i in 1:length(para)){
  
  plot_var<-para[i]
  
  p<-list()
  for(j in 1:length(type)){
    
    plot_data<-tcga_patient[which(tcga_patient$CANCER_TYPE_ACRONYM==type[j]),]
    
    mut_type<-unique(tcga_info[which(tcga_info$PATIENT_ID%in%plot_data$PATIENT_ID),"id"])
    
    print(type[j])
    print(length(mut_type))
    
    p[[j]]<-plot_OS_surv(plot_data,var=plot_var,paste0("TCGA_",type[j]),pdf_flag=FALSE)
  }
  
  comb<-arrange_ggsurvplots(p,ncol=3,nrow=3,print=FALSE)

  ggsave(paste0(output_folder,"Surv_OS_",plot_var,"_tcga_type.pdf"),comb,limitsize = FALSE,width=18,height=18)
  
  
}
#####


