# This R script performs comprehensive data preprocessing and feature engineering for training machine learning models to predict drug response based on screening data. The code systematically processes experimental data from three cancer drugs (Abemaciclib, Binimetinib, Olaparib) to prepare training datasets for predictive modeling. Key preprocessing steps include: filtering and processing screening data to identify positive and negative training sites; handling amino acid mutations by separating concatenated mutation strings into structured format; generating semi-quantitative phenotypic scores from continuous Z-score LFC values; integrating structural topology features from AlphaFold2 predictions with amino acid biophysical properties; preparing protein and amino acid embeddings from pLM; and creating train-validation splits for model development. The script implements rigorous data quality controls including removal of nonsense mutations, deduplication of training examples, and integration of knockout phenotype data. Output includes formatted training datasets with comprehensive feature sets suitable for machine learning algorithms, particularly Lasso-MARS ensemble models for drug response prediction.

library(tidyr)
library(dplyr)
# options(scipen=999)
options(stringsAsFactors = FALSE)
source("~/Other_project/CRC/reproducibility/script/function_pLM_LassoMARS.R")
setwd("~/Other_project/CRC/reproducibility/b.training_prepare")

data_folder<-"../a.data/"

drugs<-c("abem","bini","olap")
for(drug in drugs){
  
  screens<-read.table(paste0("screen_",drug,"_lib_all_potential_Anno.csv"),quote="",header=T,sep=",",fill=TRUE)
  
  screens<-screens[which(screens$tx_class=="canonical"),]
  
  neg_sgRNA<-screens[which(screens$rra==1 & nchar(screens$aa_variant.single_nt.multiple_nt.)>20 & abs(screens$zlfc)<1),"id"]
  # print(length(neg_sgRNA))
  
  negative_candidates<-screens[which(screens$id%in%neg_sgRNA),]
  
  
  negative_candidates<-negative_candidates[!grepl("Nonsense",negative_candidates$mut_type),]
  
  
  positive_site<-read.table(paste0("screen_",drug,"_lib_rra01_potential_Anno_canonical.csv"),quote="",header=T,sep=",",fill=TRUE)
  
  negative_candidates<-negative_candidates[which(negative_candidates$uniprot_id%in%positive_site$uniprot_id),]
  
  # print(nrow(negative_candidates))
  
  write.table(negative_candidates,
              paste0("screen_sgRNA_negative_candidates_",drug,".txt"),
              append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
  
  
}


########

negative_site_all<-read.table(paste0("screen_lib_editing_select_for_training_negative.txt"),quote="",header=T,sep="\t",fill=TRUE)
negative_site_all$id<-paste0(negative_site_all$sg_seq,"_",negative_site_all$sg_type)


drugs<-c("abem","bini","olap")

for(drug in drugs){
  
  drug_site<-read.table(paste0("screen_",drug,"_lib_rra01_potential_Anno_canonical.csv"),quote="",header=T,sep=",",fill=TRUE)
  print(nrow(drug_site))
  drug_site<-drug_site[,c("id","gene","sg_seq","zlfc","rra","uniprot_id","aa_variant.single_nt.multiple_nt.","mut_type.single_nt.multiple_nt.")]
  
  colnames(drug_site)<-c("id","gene","sg_seq","zlfc","rra","uniprot_id","aa_mut","mut_type")
  
  
  drug_site<-drug_site[!grepl("Nonsense",drug_site$mut_type),]
  
  drug_site_sep<-screen_site_sep(drug_site)
  
  positive_sites<-drug_site_sep[which(drug_site_sep$aa_ref!=drug_site_sep$aa_alt & drug_site_sep$aa_alt!="*"),]
  
  positive_sites$phenotype<-semi_quantitative(positive_sites$zlfc)
  
  negative_candidates<-read.table(paste0("screen_sgRNA_negative_candidates_",drug,".txt"),quote="",header=T,sep="\t",fill=TRUE)
  
  negative_candidates<-negative_candidates[,c("id","gene","sg_seq","zlfc","rra","uniprot_id","aa_variant.single_nt.multiple_nt.","mut_type.single_nt.multiple_nt.")]
  
  colnames(negative_candidates)<-c("id","gene","sg_seq","zlfc","rra","uniprot_id","aa_mut","mut_type")
  
  negative_candidates_sep<-screen_site_sep(negative_candidates)
  
  negative_candidates_sep<-negative_candidates_sep[which(negative_candidates_sep$aa_ref!=negative_candidates_sep$aa_alt &  negative_candidates_sep$aa_alt!="*"),]
  
  
  negative_candidates_sep<-negative_candidates_sep[order(abs(negative_candidates_sep$zlfc)),]
  
  tmp<-paste0(negative_candidates_sep$uniprot_id,"_",round(negative_candidates_sep$aa_pos/10))
  
  
  negative_sites<-negative_candidates_sep[!duplicated(tmp),]
  negative_sites$phenotype<-0
  
  negative_sites<-inner_join(negative_sites,negative_site_all[,c("id","aa_variant")],by=join_by("id"=="id","aa_mut"=="aa_variant"))
  
  
  training_sites<-rbind(positive_sites,negative_sites)
  
  training_sites<-training_sites[order(training_sites$gene),]
  
  training_sites<-training_sites[,-c(3,8)]
  
  training_sites<-training_sites[order(abs(training_sites$zlfc),decreasing = T),]
  
  training_sites<-training_sites[!duplicated(training_sites[,c("uniprot_id","aa_mut","phenotype")]),]
  
  write.table(training_sites,
              paste0("training_sites_",drug,".txt"),
              append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
  
}

#########
file_list<-c("abem","bini","olap")

for(i in 1:length(file_list)){
  
  ko_phenotype<-read.table(paste0("screen_BBK_sgRNA_",file_list[i],".txt"),quote="",header=T,sep="\t",fill=TRUE)
  
  ko_phenotype$KO<-semi_quantitative(ko_phenotype$zlfc)
  
  ko_phenotype<-ko_phenotype[order(abs(ko_phenotype$zlfc),decreasing = T),]
  ko_phenotype<-ko_phenotype[!duplicated(ko_phenotype$gene),]
  
  ko_phenotype<-ko_phenotype[,c("gene","zlfc","KO")]
  colnames(ko_phenotype)<-c("gene","KO_zlfc","KO")
  
  write.table(ko_phenotype,
              paste0("KO_phenotype_",file_list[i],".txt"),
              append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
  
  
}




######

load(paste0(data_folder,"AF2_WG_CF.RData"))

for(drug in c("abem","bini","olap")){
  
  training_site<-read.table(paste0("training_sites_",drug,".txt"),quote="",header=T,sep="\t",fill=TRUE)
  
  training_ko<-read.table(paste0("KO_phenotype_",drug,".txt"),quote="",header=T,sep="\t",fill=TRUE)
  
  #colnames(training_ko)<-c("gene","KO","uniprot")
  
  training_site<-left_join(training_site,training_ko,by=join_by("gene"=="gene"))
  
  training_site$KO[is.na(training_site$KO)]<-0
  
  training_site$id<-paste0(training_site$aa_pos,training_site$aa_ref)
  
  sub_CF<-WG_CF[which(WG_CF$Uniprot%in%training_site$uniprot_id),]
  
  training_site_CF<-inner_join(training_site,sub_CF,by=join_by("uniprot_id"=="Uniprot","id"=="id"))
  
  write.table(training_site_CF,
              paste0("training_site_CF_",drug,".txt"),
              append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
  
}

########


AA_property<-read.table(paste0(data_folder,"aa_property.txt"),quote="",header=T,sep="\t",fill=TRUE)

AA_property_scale<-AA_property[,c("code","nAngels","Mass","IP","HM")]

AA_property_scale$nAngels<-round(scale(AA_property$nAngels)[,1],2)

AA_property_scale$Mass<-round(scale(AA_property$Mass)[,1],2)
AA_property_scale$IP<-round(scale(AA_property$IP)[,1],2)
AA_property_scale$HM<-round(scale(AA_property$HM)[,1],2)


condition<-c("abem","bini","olap")
output_label<-"train"

for(i in 1:length(condition)){
  
  print(condition[i])
  CF_info<-read.table(paste0("training_site_CF_",condition[i],".txt"),quote="",header=T,sep="\t",fill=TRUE)
  
  CF_use<-CF_info[,c("id","uniprot_id","gene","aa_pos","aa_ref","aa_alt","zlfc","phenotype","KO","CF10","CF10RK","LD15","LD15RK")]
  
  CF_use$AF<-1
  
  CF_use<-left_join(CF_use,AA_property_scale,by=join_by("aa_ref"=="code"))
  
  CF_use<-left_join(CF_use,AA_property_scale,by=join_by("aa_alt"=="code"),suffix = c(".ref", ".alt"))
  
  
  CF_use$id<-paste0(CF_use$uniprot_id,".",CF_use$aa_pos)
  
  CF_use<-unique(CF_use[which(!is.na(CF_use$CF10)),])
  
  CF_use$split<-sample_split(CF_use)
  
  write.table(CF_use,
              paste0("mars_",output_label,"_",condition[i],"_info.txt"),
              append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
  ################# c1  ##############
  train_prot<-read.table(paste0("eb_per_protein_",condition[i],".txt"),quote="",header=F,sep="\t",fill=TRUE)
  colnames(train_prot)<-c("prot",paste0("p_",c(0:1023)))
  train_prot<-unique(train_prot)
  rownames(train_prot)<-train_prot$prot
  eb_prot<-train_prot[CF_use$uniprot_id,]
  
  print(identical(eb_prot[,1],CF_use$uniprot_id))
  print(which(is.na(eb_prot)))
  
  write.table(eb_prot,
              paste0("mars_",output_label,"_",condition[i],"_eb_prot_c1.txt"),
              append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
  ################# c2  ##############
  train_AA<-read.table(paste0("eb_per_residue_",condition[i],".txt"),quote="",header=F,sep="\t",fill=TRUE)
  colnames(train_AA)<-c("prot.aa",paste0("aa_",c(0:1023)))
  
  train_AA<-unique(train_AA)
  rownames(train_AA)<-train_AA$prot.aa
  
  eb_AA<-train_AA[CF_use$id,]
  print(which(is.na(eb_AA)))
  
  write.table(eb_AA,
              paste0("mars_",output_label,"_",condition[i],"_eb_aa_c2.txt"),
              append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
}

