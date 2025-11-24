# This R script performs comprehensive survival analysis based on the spatial proximity between cancer mutations and phosphorylation sites in COSMIC-annotated resistance genes across TCGA and MSK clinical cohorts. The code implements a specialized distance metric that calculates the average minimum distance from the three closest phosphorylation sites for each mutation, providing a robust measure of mutation-phosphosite spatial relationships. Key analytical steps include: processing clinical mutation data from TCGA and MSK datasets; calculating patient-level phosphorylation proximity metrics using the average minimum distance function; stratifying patients by KRAS G12 mutation status for subgroup analysis; generating density distribution plots to compare phosphorylation distance patterns between cohorts; and performing survival analysis using dual cutoff approaches to compare overall survival between patients with mutations close to versus distant from phosphorylation sites. This analysis enables systematic investigation of how spatial relationships between mutations and regulatory phosphorylation sites influence clinical outcomes, providing insights into the structural determinants of cancer progression and potential biomarkers for patient stratification.

avg_min_n<-function(nums,n){
  
  min_n <- sort(nums)[1:min(n,length(nums))]
  mean_min_n <- mean(min_n)
  mean_min_n
}

library(pROC)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggpubr)
library("RColorBrewer")
library(stringr)
options(scipen=999)
options(stringsAsFactors = FALSE)
setwd("~/Other_project/CRC/reproducibility/f.phos_dist_effect")
source("~/Other_project/CRC/reproducibility/script/function_pLM_LassoMARS.R")

#######

data_folder<-"../a.data/"
clinical_folder<-"../e.validation_clinical/"
output_folder<-"../o.output_figures/"


load(paste0(data_folder,"Unires_Phos_AA_PDB_neighbour.RData"))

tcga_all_eb<-read.table(paste0(clinical_folder,"TCGA_Core_all_eb_all.txt"),quote="",header=T,sep="\t",fill=TRUE)

tcga_info<-read.table(paste0(clinical_folder,"TCGA_Core_all_test_info.txt"),quote="",header=T,sep="\t",fill=TRUE)
tcga_obs<-read.table(paste0(clinical_folder,"TCGA_Core_all_patient_info.txt"),quote="",header=T,sep="\t",fill=TRUE)

msk_all_eb<-read.table(paste0(clinical_folder,"MSK_Core_eb_all.txt"),quote="",header=T,sep="\t",fill=TRUE)
msk_info<-read.table(paste0(clinical_folder,"MSK_Core_test_info.txt"),quote="",header=T,sep="\t",fill=TRUE)
msk_obs<-read.table(paste0(clinical_folder,"MSK_Core_patient_info.txt"),quote="",header=T,sep="\t",fill=TRUE)

#####

cosmic_resistence<-read.table(paste0(data_folder,"COSMIC_Resistance_Mutations.txt"),quote="",header=T,sep="\t",fill=TRUE)
cosmic_resistence_prot<-unique(cosmic_resistence$uniprot_id)
Phos_nb<-Phos_nb[which(Phos_nb$uniprot_id%in%cosmic_resistence_prot),]
########

tcga_info$id<-paste0(tcga_info$aa_pos,tcga_info$aa_ref)

kras_patient<-unique(tcga_info[which(tcga_info$uniprot_id=="P01116"& tcga_info$id=="12G"),"PATIENT_ID"])  

tcga_obs$KRAS12G<-"No"
tcga_obs$KRAS12G[which(tcga_obs$PATIENT_ID%in%kras_patient)]<-"Yes"
table(tcga_obs$KRAS12G)

tcga_aa_nb<-inner_join(unique(tcga_info[,c("PATIENT_ID","uniprot_id","id")]),Phos_nb,by=join_by("uniprot_id","id"))
print(length(unique(tcga_aa_nb$uniprot_id)))

tcga_patient_nb<-tcga_aa_nb %>%
  group_by(PATIENT_ID) %>%
  summarise(
    # PhosDist1 = min(PhosDist),
    # PhosDist2 = avg_min_n(PhosDist,2),
    PhosDist3 = avg_min_n(PhosDist,3),
    Phos_around=length(which(PhosDist<20))
  )

tcga_patient<-as.data.frame(inner_join(tcga_patient_nb,tcga_obs,by=join_by("PATIENT_ID")))

tcga_patient_nonKras<-tcga_patient[which(tcga_patient$KRAS12G=="No"),]
tcga_patient_Kras<-tcga_patient[which(tcga_patient$KRAS12G=="Yes"),]

######

msk_info$id<-paste0(msk_info$aa_pos,msk_info$aa_ref)

kras_patient<-unique(msk_info[which(msk_info$uniprot_id=="P01116"& msk_info$id=="12G"),"PATIENT_ID"])  ###4340 

msk_obs$KRAS12G<-"No"
msk_obs$KRAS12G[which(msk_obs$PATIENT_ID%in%kras_patient)]<-"Yes"
table(msk_obs$KRAS12G)

msk_aa_nb<-inner_join(unique(msk_info[,c("PATIENT_ID","uniprot_id","id")]),Phos_nb,by=join_by("uniprot_id","id"))

print(length(unique(msk_aa_nb$uniprot_id))) ##5

msk_patient_nb<-msk_aa_nb %>%
  group_by(PATIENT_ID) %>%
  summarise(
    # PhosDist1 = min(PhosDist),
    # PhosDist2 = avg_min_n(PhosDist,2),
    PhosDist3 = avg_min_n(PhosDist,3),
    Phos_around=length(which(PhosDist<20))
  )

msk_patient<-as.data.frame(inner_join(msk_patient_nb,msk_obs,by=join_by("PATIENT_ID")))

msk_patient_nonKras<-msk_patient[which(msk_patient$KRAS12G=="No"),]
msk_patient_Kras<-msk_patient[which(msk_patient$KRAS12G=="Yes"),]
#########

clinical_all<-rbind(cbind(orgin="TCGA",tcga_patient[,1:2]),cbind(orgin="MSK",msk_patient[,1:2]))

p<-ggplot(clinical_all,aes(x=PhosDist3))+
  labs(title=str_wrap(paste0("PhosDist3 distribution"),40))+ 
  xlim(c(-1,75))+
  geom_density(aes(colour=orgin),alpha=0.4,show.legend=FALSE)+
  stat_density(aes(colour=orgin),
               geom="line",position="identity")+
  theme_minimal()+
  theme(legend.position = "top") 


plot_var<-"PhosDist3"

dist_low<-20
dist_high<-20

p1<-plot_OS_surv_by_cut2(data=tcga_patient,var=plot_var,high_cut=dist_high,low_cut=dist_low,paste0("TCGA (",nrow(tcga_patient)," cases)"),pdf_flag=FALSE)

p2<-plot_OS_surv_by_cut2(data=msk_patient,var=plot_var,high_cut=dist_high,low_cut=dist_low,paste0("MSK (",nrow(msk_patient)," cases)"),pdf_flag=FALSE)

dist_low<-15
dist_high<-30

p3<-plot_OS_surv_by_cut2(data=tcga_patient,var=plot_var,high_cut=dist_high,low_cut=dist_low,paste0("TCGA (",nrow(tcga_patient)," cases)"),pdf_flag=FALSE)

p4<-plot_OS_surv_by_cut2(data=msk_patient,var=plot_var,high_cut=dist_high,low_cut=dist_low,paste0("MSK (",nrow(msk_patient)," cases)"),pdf_flag=FALSE)


comb<-arrange_ggsurvplots(list(p1,p2,p3,p4),ncol=2,nrow=2,print=FALSE)

ggsave(paste0(output_folder,"Surv_OS_",plot_var,"_clinical_COSMIC_core.pdf"),comb,limitsize = FALSE,width=10,height=8)


