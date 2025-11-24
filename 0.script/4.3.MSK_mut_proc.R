# This R script processes MSK cancer clinical data and integrates it with protein structural features for clinical outcome prediction. The code systematically prepares MSK metastatic cancer mutation data for machine learning analysis by combining clinical information with AlphaFold2-predicted structural topology metrics. Key processing steps include: mapping MSK mutation data to UniProt identifiers using gene correspondence tables; filtering and integrating structural features (CF10, CF10RK, LD15, LD15RK) for mutation positions; extracting and standardizing clinical variables including overall survival, microsatellite instability status, cancer types, and patient demographics; creating specialized annotations for KRAS G12 mutations to enable subgroup analysis; preparing machine learning-ready datasets by combining structural topology with normalized amino acid biophysical properties; and generating formatted input files for predictive modeling. This processed dataset enables the application of pre-trained stress response models to MSK patient data, facilitating the validation of computational predictions in an independent clinical cohort and supporting the development of structure-informed prognostic tools for metastatic cancer patients.

options(scipen=999)
options(stringsAsFactors = FALSE)

source("~/Other_project/CRC/reproducibility/script/function_pLM_LassoMARS.R")
setwd("~/Other_project/CRC/reproducibility/e.validation_clinical")


data_folder<-"../a.data/"
model_folder<-"../c.training/"
output_folder<-"../o.output_figures/"


prot_gene<-read.table(paste0(data_folder,"UP000005640_9606_prot_gene.txt"),quote="",header=T,sep="\t",fill=TRUE)

msk_comb<-read.table(paste0(data_folder,"msk_met_2021_mut_patient.csv"),quote="\"",header=T,sep=",",fill=TRUE)

msk_comb<-left_join(msk_comb,prot_gene,by=join_by("Hugo_Symbol"=="gene_name"))

msk_comb<-msk_comb[!is.na(msk_comb$uniprot_id),]
msk_comb<-unique(msk_comb)

table(msk_comb$STAGE_AT_DIAGNOSIS)

write.table(msk_comb,
            paste0("MSK_Pan_Can.txt"),
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)


load(file=paste0(data_folder,"AF2_WG_CF.RData"))

msk_comb$Protein_position<-as.character(msk_comb$Protein_position)
msk_cf<-inner_join(msk_comb,WG_CF,by=join_by("uniprot_id"=="Uniprot","Protein_position"=="Pos","aa_from"=="AA"))

# write.table(msk_cf,
#             paste0("MSK_Pan_Can_CF_full.txt"),
#             append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)


msk_cf_reduce<-msk_cf[,c("uniprot_id","Hugo_Symbol","pdb_chain","chain_len","Protein_position","HGVSp_Short","aa_from","aa_to","CF10","CF10RK","LD15","LD15RK","Tumor_Sample_Barcode","PATIENT_ID","AGE_AT_SEQUENCING","MSI_SCORE","MSI_TYPE","OS_MONTHS","OS_STATUS","SAMPLE_TYPE","CANCER_TYPE")]

colnames(msk_cf_reduce)[1:6]<-c("UniprotID","gene","pdb_chain","chain_len","Pos","aa_variant")


msk_cf_reduce<-unique(msk_cf_reduce)

write.table(msk_cf_reduce,
            paste0("MSK_Pan_Can_CF_info.txt"),
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)

mut_pos<-unique(msk_cf_reduce[,c("UniprotID","Pos")])


# write.table(mut_pos,
#             paste0("MSK_mut_pos_for_eb.txt"),
#             append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)



############


load(file = paste0(data_folder,"core73_eb.RData"))

AA_property<-read.table(paste0(data_folder,"aa_property.txt"),quote="",header=T,sep="\t",fill=TRUE)
AA_property_scale<-AA_property[,c("code","nAngels","Mass","IP","HM")]
AA_property_scale$nAngels<-round(scale(AA_property$nAngels)[,1],2)
AA_property_scale$Mass<-round(scale(AA_property$Mass)[,1],2)
AA_property_scale$IP<-round(scale(AA_property$IP)[,1],2)
AA_property_scale$HM<-round(scale(AA_property$HM)[,1],2)

CF_all<-read.table(paste0("MSK_Pan_Can_CF_info.txt"),quote="",header=T,sep="\t",fill=TRUE)


prot_scope<-read.table(paste0(data_folder,"Prot_scope_common.csv"),quote="",header=F,sep=",",fill=TRUE)
colnames(prot_scope)<-c("UNIPROT","overlap")

CF_all$id<-paste0(CF_all$Pos,CF_all$aa_from)

CF_use<-CF_all[,c("PATIENT_ID","id","UniprotID","Pos","aa_from","aa_to","CF10","CF10RK","LD15","LD15RK")]
colnames(CF_use)<-c("PATIENT_ID","id","uniprot_id","aa_pos","aa_ref","aa_alt","CF10","CF10RK","LD15","LD15RK")
CF_use<-unique(CF_use)

patient_info<-CF_all[,c("PATIENT_ID","AGE_AT_SEQUENCING","MSI_SCORE","MSI_TYPE","OS_MONTHS","OS_STATUS","SAMPLE_TYPE","CANCER_TYPE")]
patient_info<-unique(patient_info)
print(nrow(patient_info))
print(length(unique(patient_info$PATIENT_ID)))

kras_patient<-CF_all[which(CF_all$UniprotID=="P01116"& CF_all$id=="12G"),"PATIENT_ID"]  ### 

patient_info$KRAS12G<-"No"
patient_info$KRAS12G[which(patient_info$PATIENT_ID%in%kras_patient)]<-"Yes"
table(patient_info$KRAS12G)

write.table(patient_info,
            paste0("MSK_Core_patient_info.txt"),
            append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# table(patient_info$DFS_EVENT)
# table(patient_info$OS_STATUS)

CF_use<-left_join(CF_use,AA_property_scale,by=join_by("aa_ref"=="code"))
CF_use<-left_join(CF_use,AA_property_scale,by=join_by("aa_alt"=="code"),suffix = c(".ref", ".alt"))
CF_use$id<-paste0(CF_use$uniprot_id,".",CF_use$aa_pos)


CF_use<-unique(CF_use[which(!is.na(CF_use$CF10) & !is.na(CF_use$LD15)),])
print(which(is.na(CF_use)))


prot_in<-prot_scope$UNIPROT
prefix<-paste0("MSK_Core")

eb_CF_to_model_input(CF_use,prot_in,prot_eb,AA_eb,prefix)

