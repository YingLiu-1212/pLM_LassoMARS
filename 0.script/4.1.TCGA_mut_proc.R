# This R script processes TCGA (The Cancer Genome Atlas) pan-cancer mutation data and integrates it with protein structural features from AlphaFold2 predictions to create comprehensive datasets for clinical outcome analysis. The code performs systematic data integration by combining clinical information from TCGA patients with structural topology metrics (CF10, CF10RK, LD15, LD15RK) and amino acid biophysical properties. Key processing steps include: filtering and merging TCGA mutation data with protein-gene correspondence tables to ensure accurate UniProt identifier mapping; integrating structural features from AlphaFold2 predictions for mutation positions; extracting and standardizing clinical variables including overall survival (OS), progression-free survival (PFS), cancer types, and treatment information; creating specialized patient annotations such as KRAS G12 mutation status for subgroup analysis; preparing machine learning-ready datasets by combining structural features with normalized amino acid properties; and generating formatted input files for predictive modeling. 

library(dplyr)
options(scipen=999)
options(stringsAsFactors = FALSE)
source("~/Other_project/CRC/reproducibility/script/function_pLM_LassoMARS.R")
setwd("~/Other_project/CRC/reproducibility/e.validation_clinical")

data_folder<-"../a.data/"

prot_gene<-read.table(paste0(data_folder,"UP000005640_9606_prot_gene.txt"),quote="",header=T,sep="\t",fill=TRUE)

tcga_comb<-read.table(paste0(data_folder,"all_tcga_pan_can_atlas_2018.csv"),quote="",header=T,sep=",",fill=TRUE)

tcga_comb<-left_join(tcga_comb,prot_gene,by=join_by("Hugo_Symbol"=="gene_name"))

tcga_comb<-tcga_comb[!is.na(tcga_comb$uniprot_id),]
tcga_comb<-unique(tcga_comb)

# table(tcga_comb$PATH_T_STAGE)
# table(tcga_comb$TREATMENT_TYPE)

write.table(tcga_comb,
            paste0("TCGA_Pan_Can_atlas_2018.txt"),
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)


load(file=paste0(data_folder,"AF2_WG_CF.RData"))

WG_CF$Pos<-as.character(WG_CF$Pos)
tcga_cf<-inner_join(tcga_comb,WG_CF,by=join_by("uniprot_id"=="Uniprot","Protein_position"=="Pos","aa_from"=="AA"))

# write.table(tcga_cf,
#             paste0("TCGA_Pan_Can_CF_full.txt"),
#             append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)


tcga_cf_reduce<-tcga_cf[,c("uniprot_id","Hugo_Symbol","pdb_chain","chain_len","Protein_position","HGVSp_Short","aa_from","aa_to","CF10","CF10RK","LD15","LD15RK","Tumor_Sample_Barcode","PATIENT_ID","CANCER_TYPE_ACRONYM","AGE","AJCC_PATHOLOGIC_TUMOR_STAGE","OS_STATUS","OS_MONTHS","PFS_STATUS","PFS_MONTHS","START_DATE","STOP_DATE","TREATMENT_TYPE","AGENT")]

colnames(tcga_cf_reduce)[1:6]<-c("UniprotID","gene","pdb_chain","chain_len","Pos","aa_variant")


tcga_cf_reduce<-unique(tcga_cf_reduce)

write.table(tcga_cf_reduce,
            paste0("TCGA_Pan_Can_CF_info.txt"),
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)

load(file = paste0(data_folder,"core73_eb.RData"))

######

AA_property<-read.table(paste0(data_folder,"aa_property.txt"),quote="",header=T,sep="\t",fill=TRUE)
AA_property_scale<-AA_property[,c("code","nAngels","Mass","IP","HM")]
AA_property_scale$nAngels<-round(scale(AA_property$nAngels)[,1],2)
AA_property_scale$Mass<-round(scale(AA_property$Mass)[,1],2)
AA_property_scale$IP<-round(scale(AA_property$IP)[,1],2)
AA_property_scale$HM<-round(scale(AA_property$HM)[,1],2)

CF_all<-read.table(paste0("TCGA_Pan_Can_CF_info.txt"),quote="",header=T,sep="\t",fill=TRUE)

# ajcc_stage<-CF_all$AJCC_PATHOLOGIC_TUMOR_STAGE
# CF_all$Stage<-"0"
# CF_all$Stage[which(ajcc_stage%in%c("STAGE 0","STAGE I","STAGE I/II (NOS)","STAGE IA","STAGE IB"))]<-"I"
# CF_all$Stage[which(ajcc_stage%in%c("STAGE II","STAGE IIA","STAGE IIB","STAGE IIC"))]<-"II"
# CF_all$Stage[which(ajcc_stage%in%c("STAGE III","STAGE IIIA","STAGE IIIB","STAGE IIIC"))]<-"III"
# CF_all$Stage[which(ajcc_stage%in%c("STAGE IV","STAGE IVA","STAGE IVB","STAGE IVC"))]<-"IV"
# seperate_stage<-unique(CF_all$Stage)


#######
prot_scope<-read.table(paste0(data_folder,"Prot_scope_common.csv"),quote="",header=F,sep=",",fill=TRUE)
colnames(prot_scope)<-c("UNIPROT","overlap")

CF_all$id<-paste0(CF_all$Pos,CF_all$aa_from)

CF_use<-CF_all[,c("PATIENT_ID","id","UniprotID","Pos","aa_from","aa_to","CF10","CF10RK","LD15","LD15RK")]
colnames(CF_use)<-c("PATIENT_ID","id","uniprot_id","aa_pos","aa_ref","aa_alt","CF10","CF10RK","LD15","LD15RK")
CF_use<-unique(CF_use)

patient_info<-CF_all[,c("PATIENT_ID","Stage","PFS_STATUS","PFS_MONTHS","OS_MONTHS","OS_STATUS","CANCER_TYPE_ACRONYM","AGE")]
patient_info<-unique(patient_info)
# print(nrow(patient_info))
# print(length(unique(patient_info$PATIENT_ID)))


kras_patient<-CF_all[which(CF_all$UniprotID=="P01116"& CF_all$id=="12G"),"PATIENT_ID"]  ### 814

patient_info$KRAS12G<-"No"
patient_info$KRAS12G[which(patient_info$PATIENT_ID%in%kras_patient)]<-"Yes"
table(patient_info$KRAS12G)

write.table(patient_info,
            paste0("TCGA_Core_all_patient_info.txt"),
            append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)


CF_use<-left_join(CF_use,AA_property_scale,by=join_by("aa_ref"=="code"))
CF_use<-left_join(CF_use,AA_property_scale,by=join_by("aa_alt"=="code"),suffix = c(".ref", ".alt"))
CF_use$id<-paste0(CF_use$uniprot_id,".",CF_use$aa_pos)


CF_use<-unique(CF_use[which(!is.na(CF_use$CF10) & !is.na(CF_use$LD15)),])
print(which(is.na(CF_use)))


prot_in<-prot_scope$UNIPROT
prefix<-paste0("TCGA_Core_all")

eb_CF_to_model_input(CF_use,prot_in,prot_eb,AA_eb,prefix)



