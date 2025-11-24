# This R script performs comparative spatial analysis of drug resistance mutations relative to phosphorylation sites using data from COSMIC database and experimental screening studies. The code analyzes the structural proximity between known resistance-associated mutations and regulatory phosphorylation sites to investigate potential mechanisms of drug resistance. Key analytical steps include: processing COSMIC database annotations to identify clinically documented resistance mutations; integrating experimental screening data from Abemaciclib and Binimetinib studies to identify functionally validated resistance sites; calculating minimum distances between resistance mutations and phosphorylation sites using pre-computed structural neighbor data; generating comparative boxplot visualizations with statistical testing to examine spatial distribution patterns; and conducting cross-dataset comparison between curated clinical mutations (COSMIC) and experimentally derived resistance sites (screening data). 

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
output_folder<-"../o.output_figures/"
screen_folder<-"../c.training/"

load(paste0(data_folder,"Unires_Phos_AA_PDB_neighbour.RData"))

cosmic_aa_phenotype<-read.table(paste0(data_folder,"COSMIC_Resistance_Mutations.txt"),quote="",header=T,sep="\t",fill=TRUE)



cosmic_aa_phenotype$id<-paste0(cosmic_aa_phenotype$aa_pos,cosmic_aa_phenotype$aa_ref)
cosmic_aa_phenotype$resistance<-1

cosmic_aa_phenotype<-unique(cosmic_aa_phenotype[,c("uniprot_id","id","GENE_SYMBOL","resistance")])

cosmic_Phos_nb<-Phos_nb[which(Phos_nb$uniprot_id%in%unique(cosmic_aa_phenotype$uniprot_id)),]

cosmic_aa_phenotype_Phos<-left_join(cosmic_Phos_nb,cosmic_aa_phenotype,by=join_by("uniprot_id","id"))


cosmic_aa_phenotype_Phos$Classification<-"Ambiguous"
cosmic_aa_phenotype_Phos$Classification[which(cosmic_aa_phenotype_Phos$resistance==1)]<-"Resistance"

table(cosmic_aa_phenotype_Phos$Classification)

length(unique(cosmic_aa_phenotype_Phos$uniprot_id)) ##26
########
abem_phenotype<-read.table(paste0(screen_folder,"mars_train_abem_info.txt"),quote="",header=T,sep="\t",fill=TRUE)
bini_phenotype<-read.table(paste0(screen_folder,"mars_train_bini_info.txt"),quote="",header=T,sep="\t",fill=TRUE)

screen_AB_aa_phenotype<-rbind(abem_phenotype,bini_phenotype)

screen_AB_aa_phenotype$id<-paste0(screen_AB_aa_phenotype$aa_pos,screen_AB_aa_phenotype$aa_ref)

screen_AB_aa_phenotype<-screen_AB_aa_phenotype[which(screen_AB_aa_phenotype$KO<1),]

screen_AB_aa_phenotype<-unique(screen_AB_aa_phenotype)

AB_aa_phenotype<-screen_AB_aa_phenotype %>%
  group_by(uniprot_id,id) %>%
  summarise(
    phenotype = sum(phenotype),
  )

AB_aa_phenotype<-AB_aa_phenotype[which(AB_aa_phenotype$phenotype>=1),]

length(unique(AB_aa_phenotype$uniprot_id)) ## 56

AB_Phos_nb<-Phos_nb[which(Phos_nb$uniprot_id%in%unique(AB_aa_phenotype$uniprot_id)),]

AB_aa_phenotype_Phos<-left_join(AB_Phos_nb,AB_aa_phenotype,by=join_by("uniprot_id","id"))

AB_aa_phenotype_Phos$Classification<-"Ambiguous"
AB_aa_phenotype_Phos$Classification[which(AB_aa_phenotype_Phos$phenotype>=1)]<-"Resistance"

table(AB_aa_phenotype_Phos$Classification)
length(unique(AB_aa_phenotype_Phos$uniprot_id)) ## 27


######



plot_list<-list()

cosmic_aa_phenotype_Phos$Classification<-factor(cosmic_aa_phenotype_Phos$Classification,levels=c("Resistance","Ambiguous"))

AB_aa_phenotype_Phos$Classification<-factor(AB_aa_phenotype_Phos$Classification,levels=c("Resistance","Ambiguous"))

plot_list[[1]]<-boxplot_var_dist(cosmic_aa_phenotype_Phos,"PhosDist",plot_title="Distance of mutated AAs from phosphorylated AAs in COSMIC annotated resistent mutations")

plot_list[[2]]<-boxplot_var_dist(AB_aa_phenotype_Phos,"PhosDist",plot_title="Distance of resistent mutations from phosphorylated AAs in Abem&Bini screens")

comb<-ggarrange(plotlist=plot_list,ncol=2,nrow=1,common.legend = T)

ggsave(paste0(output_folder,"PDB_AA_Unires_Cosmic_AB_box.pdf"),comb,limitsize = FALSE,width=8,height=4)





