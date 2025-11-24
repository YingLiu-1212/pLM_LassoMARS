# This R script performs comprehensive processing and analysis of cross-drug resistance screening data from published CRISPR-based functional genomics studies. The code systematically processes drug resistance phenotypes across nine cancer therapeutics (DabCet, trametinib, pictilisib, sotorasib, adagrasib, gefitinib, osimertinib, olaparib, niraparib) to prepare datasets for machine learning prediction of drug resistance mechanisms. Key processing steps include: parsing and filtering screening data to identify missense mutations with significant resistance or sensitivity phenotypes; implementing rigorous quality controls using FDR thresholds and log2 fold change variability filters; integrating structural topology features from AlphaFold2 predictions with amino acid biophysical properties; preparing protein and amino acid embeddings for model input; and generating formatted datasets for both full protein sets and core protein subsets. The script handles complex mutation annotation data by separating concatenated amino acid variant strings and ensures data integrity through multiple validation checks. Output includes comprehensive phenotype-structure datasets suitable for cross-drug resistance prediction using pre-trained machine learning models, enabling systematic investigation of shared and drug-specific resistance mechanisms across multiple therapeutic classes.

library(tidyr)
# library(Biostrings)
library(dplyr)
options(scipen=999)
options(stringsAsFactors = FALSE)
source("~/Other_project/CRC/reproducibility/script/function_pLM_LassoMARS.R")
setwd("~/Other_project/CRC/reproducibility/d.validation_cross_drug")


data_folder<-"../a.data/"
model_folder<-"../c.training/"
output_folder<-"../o.output_figures/"

data_label<-"cross_drug"

#### published screen data process #######

Druglib<-read.table(paste0(data_folder,data_label,"_screen_data.txt"),header=TRUE,sep="\t",fill=TRUE)
Druglib<-Druglib[which(!is.na(Druglib$CRISPR.PAM.Sequence)),]

Druglib<-Druglib[which(Druglib$most_severe_consequence%in%c("missense")),]

Druglib$id<-paste0(Druglib$CRISPR.PAM.Sequence,"_",Druglib$editor)

#table(Druglib$most_severe_consequence)
#length(unique(Druglib$id))

info<-colnames(Druglib)
which(grepl("FDR_",info))
which(grepl("_resistance",info))
which(grepl("variant_class_",info))

Drug_screens<-Druglib[,c(121,which(grepl("FDR_",info)),which(grepl("_resistance",info)),which(grepl("variant_class_",info)))]

Drug_L2FC<-as.matrix(Druglib[,c(which(grepl("L2FC_",info)))])
rownames(Drug_L2FC)<-Druglib[,121]

range_L2FC<-rowMaxs(Drug_L2FC,value = TRUE)-rowMins(Drug_L2FC,value = TRUE)

flag<-range_L2FC>3 & range_L2FC<4
fdr_negative<-0.5

length(which(flag))
#######

drug1<-Drug_screens[,c("id","FDR_MAGeCK_HT29_DebCet","debcet_resistance","variant_class_DabCet")]
# tmp<-drug1[which(drug1$FDR_MAGeCK_HT29_DebCet>0.2 & flag),]

drug1_phen<-data.frame(id=drug1[,"id"],drug="DabCet")
drug1_phen$resistance<-NA
drug1_phen$resistance[which(drug1$debcet_resistance==TRUE)]<-1
drug1_phen$resistance[which(drug1$FDR_MAGeCK_HT29_DebCet>fdr_negative & flag)]<-0
drug1_phen$resistance[which(drug1$variant_class_DabCet=="drug-sensitising")]<-0
drug1_phen<-drug1_phen[which(!is.na(drug1_phen$resistance)),]
table(drug1_phen$resistance)

####

drug2<-Drug_screens[,c("id","FDR_MAGeCK_HT29_Tram","trametinib_resistance","variant_class_Tram")]
drug2_phen<-data.frame(id=drug2[,"id"],drug="trametinib")
drug2_phen$resistance<-NA
drug2_phen$resistance[which(drug2$trametinib_resistance==TRUE)]<-1
drug2_phen$resistance[which(drug2$FDR_MAGeCK_HT29_Tram>fdr_negative & flag)]<-0
drug2_phen$resistance[which(drug2$variant_class_Tram=="drug-sensitising")]<-0
drug2_phen<-drug2_phen[which(!is.na(drug2_phen$resistance)),]
table(drug2_phen$resistance)

# ####
drug3<-Drug_screens[,c("id","FDR_MAGeCK_HT29_Pict","pictilisib_resistance","variant_class_Pict")]
drug3_phen<-data.frame(id=drug3[,"id"],drug="pictilisib")
drug3_phen$resistance<-NA

drug3_phen$resistance[which(drug3$pictilisib_resistance==TRUE)]<-1
drug3_phen$resistance[which(drug3$FDR_MAGeCK_HT29_Pict>fdr_negative & flag)]<-0
drug3_phen$resistance[which(drug3$variant_class_Pict=="drug-sensitising")]<-0
drug3_phen<-drug3_phen[which(!is.na(drug3_phen$resistance)),]
table(drug3_phen$resistance)

####
drug4<-Drug_screens[,c("id","FDR_MAGeCK_H23_Sotor","sotorasib_resistance","variant_class_Sotor")]
drug4_phen<-data.frame(id=drug4[,"id"],drug="sotorasib")
drug4_phen$resistance<-NA

drug4_phen$resistance[which(drug4$sotorasib_resistance==TRUE)]<-1
drug4_phen$resistance[which(drug4$FDR_MAGeCK_H23_Sotor>fdr_negative & flag)]<-0
drug4_phen$resistance[which(drug4$variant_class_Sotor=="drug-sensitising")]<-0
drug4_phen<-drug4_phen[which(!is.na(drug4_phen$resistance)),]
table(drug4_phen$resistance)
####

drug5<-Drug_screens[,c("id","FDR_MAGeCK_H23_Adag","adagrasib_resistance","variant_class_Adag")]
drug5_phen<-data.frame(id=drug5[,"id"],drug="adagrasib")
drug5_phen$resistance<-NA
drug5_phen$resistance[which(drug5$adagrasib_resistance==TRUE)]<-1
drug5_phen$resistance[which(drug5$FDR_MAGeCK_H23_Adag>fdr_negative & flag)]<-0
drug5_phen$resistance[which(drug5$variant_class_Adag=="drug-sensitising")]<-0
drug5_phen<-drug5_phen[which(!is.na(drug5_phen$resistance)),]
table(drug5_phen$resistance)

####
drug6<-Drug_screens[,c("id","FDR_MAGeCK_PC9_Gefit","gefitinib_resistance","variant_class_Gefit")]
drug6_phen<-data.frame(id=drug6[,"id"],drug="gefitinib")
drug6_phen$resistance<-NA

drug6_phen$resistance[which(drug6$gefitinib_resistance==TRUE)]<-1
drug6_phen$resistance[which(drug6$FDR_MAGeCK_PC9_Gefit>fdr_negative & flag)]<-0
drug6_phen$resistance[which(drug6$variant_class_Gefit=="drug-sensitising")]<-0
drug6_phen<-drug6_phen[which(!is.na(drug6_phen$resistance)),]
table(drug6_phen$resistance)
####

drug7<-Drug_screens[,c("id","FDR_MAGeCK_PC9_Osim","osimertinib_resistance","variant_class_Osim")]
drug7_phen<-data.frame(id=drug7[,"id"],drug="osimertinib")
drug7_phen$resistance<-NA

drug7_phen$resistance[which(drug7$osimertinib_resistance==TRUE)]<-1
drug7_phen$resistance[which(drug7$FDR_MAGeCK_PC9_Osim>fdr_negative & flag)]<-0
drug7_phen$resistance[which(drug7$variant_class_Osim=="drug-sensitising")]<-0
drug7_phen<-drug7_phen[which(!is.na(drug7_phen$resistance)),]
table(drug7_phen$resistance)

####
drug8<-Drug_screens[,c("id","FDR_MAGeCK_MHHES1_Olap","olaparib_resistance","variant_class_Olap")]
drug8_phen<-data.frame(id=drug8[,"id"],drug="olaparib")
drug8_phen$resistance<-NA

drug8_phen$resistance[which(drug8$olaparib_resistance==TRUE)]<-1
drug8_phen$resistance[which(drug8$FDR_MAGeCK_MHHES1_Ola>fdr_negative & flag)]<-0
drug8_phen$resistance[which(drug8$variant_class_Olap=="drug-sensitising")]<-0
drug8_phen<-drug8_phen[which(!is.na(drug8_phen$resistance)),]
table(drug8_phen$resistance)

######
drug9<-Drug_screens[,c("id","FDR_MAGeCK_MHHES1_Nirap","niraparib_resistance","variant_class_Nirap")]
drug9_phen<-data.frame(id=drug9[,"id"],drug="niraparib")
drug9_phen$resistance<-NA

drug9_phen$resistance[which(drug9$niraparib_resistance==TRUE)]<-1
drug9_phen$resistance[which(drug9$FDR_MAGeCK_MHHES1_Nirap>fdr_negative & flag)]<-0
drug9_phen$resistance[which(drug9$variant_class_Nirap=="drug-sensitising")]<-0
drug9_phen<-drug9_phen[which(!is.na(drug9_phen$resistance)),]
table(drug9_phen$resistance)

######
drug_phen<-rbind(drug1_phen,drug2_phen,drug3_phen,drug4_phen,drug5_phen,drug6_phen,drug7_phen,drug8_phen,drug9_phen)


write.table(drug_phen,
            paste0(data_label,"_screen_phenotype.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = T)


#### lib editing info process #######

load(file=paste0(data_folder,"AF2_WG_CF.RData"))

lib_info<-read.table(paste0(data_folder,data_label,"_sglib.csv"),header=T,sep=",",fill=TRUE)

edit_AA<-read.table(paste0(data_folder,data_label,"_all_editing_at_aa_confident.csv"),header=T,sep=",",fill=TRUE)
edit_AA$editor<-toupper(edit_AA$editor)
edit_AA$id<-paste0(edit_AA$sg_pam,"_",edit_AA$editor)

edit_AA_info<-inner_join(edit_AA,lib_info,by=join_by("id","editor"))


edit_AA_reduced<-edit_AA_info[,c("sg_pam","editor","uniprot_id","aa_variant.single_nt.multiple_nt.","mut_type.single_nt.multiple_nt.","CDS_gene","tx_class","Pos","aa_ref","Change")]
edit_AA_reduced<-edit_AA_reduced[which(!is.na(edit_AA_reduced$sg_pam)),]

edit_AA_reduced<-edit_AA_reduced[which(edit_AA_reduced$tx_class=="canonical"),]

edit_AA_reduced<-unique(edit_AA_reduced)

colnames(edit_AA_reduced)[8:10]<-c("Pos_org","aa_ref_org","Change_org")

write.table(edit_AA_reduced,
            paste0(data_label,"_edit_AA_canonical.txt"),
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)

edit_AA_reduced$id<-paste0(edit_AA_reduced$sg_pam,"_",edit_AA_reduced$editor)


drug_phenotype<-read.table(paste0(data_label,"_screen_phenotype.csv"),quote="",header=T,sep=",",fill=TRUE)

edit_AA_phenotype<-full_join(edit_AA_reduced,drug_phenotype,by="id")

edit_AA_phenotype<-edit_AA_phenotype[which(!is.na(edit_AA_phenotype$resistance)),]
edit_AA_phenotype<-edit_AA_phenotype[which(!is.na(edit_AA_phenotype$uniprot_id)),]

edit_AA_phenotype<-edit_AA_phenotype[,c("uniprot_id","CDS_gene","aa_variant.single_nt.multiple_nt.","drug","resistance")]
colnames(edit_AA_phenotype)[3]<-"aa_variant"
edit_AA_phenotype<-edit_AA_phenotype[order(edit_AA_phenotype$resistance,decreasing = T),]

edit_AA_phenotype<-unique(edit_AA_phenotype)

edit_AA_sep<-separate_longer_delim(edit_AA_phenotype,c(aa_variant), delim = "|")
edit_AA_sep<-separate_longer_delim(edit_AA_sep,c(aa_variant), delim = ";")
edit_AA_sep<-edit_AA_sep[which(edit_AA_sep$aa_variant!=""),]

aa_mut<-edit_AA_sep$aa_variant
edit_AA_sep$aa_from<-substr(aa_mut,1,1)
edit_AA_sep$aa_to<-substr(aa_mut,nchar(aa_mut),nchar(aa_mut))
edit_AA_sep$aa_pos<-substr(aa_mut,2,nchar(aa_mut)-1)

edit_AA_sep$id<-paste0(edit_AA_sep$aa_pos,edit_AA_sep$aa_from)

edit_AA_sep<-edit_AA_sep[which(edit_AA_sep$aa_from!=edit_AA_sep$aa_to),]

edit_AA_sep<-edit_AA_sep[order(edit_AA_sep$resistance,decreasing = T),]
edit_AA_sep<-edit_AA_sep[which(!duplicated(edit_AA_sep[,-5])),]

AA_phenotype_CF<-inner_join(edit_AA_sep,WG_CF,by=join_by("uniprot_id"=="Uniprot","id"=="id"))

write.table(AA_phenotype_CF,
            paste0(data_label,"_edit_AA_phenotype_AF2CF.txt"),
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)


#### model prediction prepare #######

prot_eb<-read.table(paste0(data_folder,"embeddings_per_protein.txt"),header=F,quote="",sep="\t",fill=TRUE)

AA_eb<-read.table(paste0(data_folder,"cross_drug_eb_per_residue_all.txt"),header=F,quote="",sep="\t",fill=TRUE)


AA_property<-read.table(paste0(data_folder,"aa_property.txt"),quote="",header=T,sep="\t",fill=TRUE)
AA_property_scale<-AA_property[,c("code","nAngels","Mass","IP","HM")]
AA_property_scale$nAngels<-round(scale(AA_property$nAngels)[,1],2)
AA_property_scale$Mass<-round(scale(AA_property$Mass)[,1],2)
AA_property_scale$IP<-round(scale(AA_property$IP)[,1],2)
AA_property_scale$HM<-round(scale(AA_property$HM)[,1],2)


CF_info<-read.table(paste0(data_label,"_edit_AA_phenotype_AF2CF.txt"),quote="",header=T,sep="\t",fill=TRUE)

prot_in<-unique(CF_info$uniprot_id)

CF_info<-CF_info[,c("id","uniprot_id","CDS_gene","aa_pos","aa_from","aa_to","CF10","CF10RK","LD15","LD15RK","drug","resistance")]

colnames(CF_info)<-c("id","uniprot_id","gene","aa_pos","aa_ref","aa_alt","CF10","CF10RK","LD15","LD15RK","drug","resistance")

CF_use<-CF_info
length(unique(CF_use$uniprot_id))
table(CF_use$drug)

CF_use<-left_join(CF_use,AA_property_scale,by=join_by("aa_ref"=="code"))
CF_use<-left_join(CF_use,AA_property_scale,by=join_by("aa_alt"=="code"),suffix = c(".ref", ".alt"))
CF_use$id<-paste0(CF_use$uniprot_id,".",CF_use$aa_pos)
CF_use<-unique(CF_use[which(!is.na(CF_use$CF10)),])
print(which(is.na(CF_use$CF10)))

write.table(CF_use,
            paste0(data_label,"_full_phenotype_info.txt"),
            append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

####
prot_scope<-read.table(paste0(data_folder,"Prot_scope_common.csv"),quote="",header=F,sep=",",fill=TRUE)
colnames(prot_scope)<-c("UNIPROT","overlap")

CF_core<-CF_use[which(CF_use$uniprot_id%in%prot_scope[,1]),]
condition<-intersect(names(which(table(CF_core$drug)>20)),names(which(table(CF_core[CF_core$resistance==1,"drug"])>10)))


# condition<-names(which(table(CF_use$drug)>20))
for(i in 1:length(condition)){
  CF_use_phe<-CF_use[which(CF_use$drug==condition[i]),]
  CF_use_phe<-unique(CF_use_phe)
  print(condition[i])
  
  print(table(CF_use_phe$resistance))
  
  if(length(which(CF_use_phe$resistance==1))>20){
    
    prefix<-paste0("Drug_",condition[i],"_full")
    
    eb_CF_to_model_input(CF_use_phe,prot_in,prot_eb,AA_eb,prefix)
  }
}




CF_use<-CF_use[which(CF_use$uniprot_id%in%prot_scope[,1]),]
length(unique(CF_use$uniprot_id))
table(CF_use$drug)

write.table(CF_use,
            paste0(data_label,"_Core_phenotype_info.txt"),
            append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

for(i in 1:length(condition)){
  
  CF_use_phe<-CF_use[which(CF_use$drug==condition[i]),]
  
  CF_use_phe<-unique(CF_use_phe)
  
  prot_in<-prot_scope$UNIPROT
  
  if(length(which(CF_use_phe$resistance==1))>10){
    
    prefix<-paste0("Drug_",condition[i],"_Core")
    
    eb_CF_to_model_input(CF_use_phe,prot_in,prot_eb,AA_eb,prefix)
    
    
  }
  
}

#### model prediction #######


load(paste0(model_folder,"pLM_LassoMARS.RData"))


drug_phenotype<-read.table(paste0(data_label,"_Core_phenotype_info.txt"),quote="",header=T,sep="\t",fill=TRUE)
unique(drug_phenotype$gene)

drug_full_phenotype<-read.table(paste0(data_label,"_full_phenotype_info.txt"),quote="",header=T,sep="\t",fill=TRUE)
unique(drug_full_phenotype$gene)

para<-c(
  "abem_pred",
  "bini_pred",
  "olap_pred",
  "StressResponse")

plot_drug_box(drug_phenotype,para,"./","Core",condition,output_folder)
plot_drug_box(drug_full_phenotype,para,"./","full",condition,output_folder)


