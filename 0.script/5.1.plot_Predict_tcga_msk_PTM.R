# This R script performs comparative analysis of the spatial relationship between drug resistance-associated mutations and phosphorylation sites in kinase proteins across TCGA and MSK clinical cohorts. The code analyzes how mutations predicted to confer drug resistance are distributed relative to known phosphorylation sites in a curated set of kinases. Using pre-computed structural distance data from UniProt and PDB resources, the script calculates minimum distances between resistance mutations and phosphorylation sites, then generates comparative boxplot visualizations with statistical testing. The analysis specifically examines whether resistance mutations in kinases show distinct spatial patterns relative to regulatory phosphorylation sites, potentially revealing mechanisms by which mutations disrupt kinase signaling and confer drug resistance. 

library(dplyr)
library(ggplot2)
library(ggpubr)
library("RColorBrewer")
library(stringr)
options(scipen=999)
options(stringsAsFactors = FALSE)
setwd("~/Other_project/CRC/reproducibility/f.phos_dist_effect")
source("~/Other_project/CRC/reproducibility/script/function_pLM_LassoMARS.R")


data_folder<-"../a.data/"
clinical_folder<-"../e.validation_clinical/"
output_folder<-"../o.output_figures/"

load(paste0(data_folder,"Unires_Phos_AA_PDB_neighbour.RData"))

## kinase_list from a.data/GO_for_Core73.txt ######
kinase_list<-c("Q15831", "P36507", "Q15759", "P28482", "P10398", "P42345", "P15056", "Q00535", "P04049", "O75914", "Q16539", "P27361", "P51812", "Q02750", "P19784", "P53778", "P24941", "P06493", "P45983", "Q00526", "Q13177", "P45985", "O14733", "Q13153", "P46734")

tcga_aa_phenotype<-read.table(paste0(clinical_folder,"tcga_core73_prediction_info.txt"),quote="",header=T,sep="\t",fill=TRUE)

# tcga_aa_phenotype<-tcga_aa_phenotype[which(tcga_aa_phenotype$uniprot_id%in%cosmic_resistence_prot),]
tcga_aa_phenotype<-tcga_aa_phenotype[which(tcga_aa_phenotype$uniprot_id%in%kinase_list),]

tcga_aa_phenotype_Phos<-resistence_to_Phos_dist(tcga_aa_phenotype,Phos_nb)

msk_aa_phenotype<-read.table(paste0(clinical_folder,"msk_core73_prediction_info.txt"),quote="",header=T,sep="\t",fill=TRUE)

# msk_aa_phenotype<-msk_aa_phenotype[which(msk_aa_phenotype$uniprot_id%in%cosmic_resistence_prot),]

msk_aa_phenotype<-msk_aa_phenotype[which(msk_aa_phenotype$uniprot_id%in%kinase_list),]

msk_aa_phenotype_Phos<-resistence_to_Phos_dist(msk_aa_phenotype,Phos_nb)

length(unique(tcga_aa_phenotype_Phos$uniprot_id)) 
length(unique(msk_aa_phenotype_Phos$uniprot_id)) 




tcga_aa_phenotype_Phos$Classification<-factor(tcga_aa_phenotype_Phos$Classification,levels=c("Resistance","Ambiguous"))

msk_aa_phenotype_Phos$Classification<-factor(msk_aa_phenotype_Phos$Classification,levels=c("Resistance","Ambiguous"))


plot_list<-list()
plot_list[[1]]<-boxplot_var_dist(tcga_aa_phenotype_Phos,"PhosDist",plot_title="Distance of mutated AAs from phosphorylated AAs in TCGA patients")

plot_list[[2]]<-boxplot_var_dist(msk_aa_phenotype_Phos,"PhosDist",plot_title="Distance of mutated AAs from phosphorylated AAs in MSK patients")

comb<-ggarrange(plotlist=plot_list,ncol=2,nrow=1,common.legend = T)

ggsave(paste0(output_folder,"PDB_AA_Unires_clinical_prediction_boxplot.pdf"),comb,limitsize = FALSE,width=8,height=4)





