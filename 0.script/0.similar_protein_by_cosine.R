# This R script performs protein similarity analysis using cosine distance metrics to expand drug screening targets based on structural and functional similarities. The code processes pre-computed protein cosine distance matrices to identify similar proteins for three cancer drugs (Abemaciclib, Binimetinib, Olaparib) using a defined cosine similarity threshold. For each drug's known target proteins, the script identifies structurally similar proteins by filtering the cosine distance matrix, then aggregates and counts the expanded protein sets. The analysis generates comprehensive output files containing the original targets along with their similar protein counterparts and corresponding gene mappings. Additionally, the script performs comparative analysis across all three drug target sets by identifying common proteins shared between them and creates a Venn diagram visualization to illustrate the overlaps and unique proteins in each drug's expanded target scope. 

library(tidyr)
library(dplyr)
library(stringr)
library(VennDiagram) 
options(scipen=9)
options(stringsAsFactors = FALSE)


setwd("~/Other_project/CRC/reproducibility")

data_folder<-"./a.data/"
figure_folder<-"./o.output_figures/"

load(paste0(data_folder,"prot_cosine_dist.RData"))
prot_gene<-read.table(paste0(data_folder,"UP000005640_9606_prot_gene.txt"),header=T,quote="",sep="\t",fill=TRUE)

cosine_cut<-0.06
file_list<-c("abem","bini","olap")

for(p in 1:length(file_list)){
  
  training_prot<-read.table(paste0(data_folder,"screen_prot_",file_list[p],".txt"),header=T,quote="",sep="\t",fill=TRUE)
  
  screen_prot<-training_prot$uniprot_id
  prot_eb_cosine<-cosine_dist_all[,which(colnames(cosine_dist_all)%in%training_prot$uniprot_id)]
  
  screen_gene_extend<-NULL
  prot_scope_comb<-NULL
  similar_prot<-rep("",length(screen_prot))
  similar_gene<-rep("",length(screen_prot))
  extend_n<-rep(0,length(screen_prot))
  
  for(i in 1:length(screen_prot)){
    
    cos_dist<-prot_eb_cosine[,screen_prot[i]]
    
    sim_prot<-rownames(prot_eb_cosine)[which(cos_dist<cosine_cut)]
    
    extend_n[i]<-length(sim_prot)
    
    sim_gene<-prot_gene[which(prot_gene$uniprot_id%in%sim_prot),]
    
    screen_gene_extend<-rbind(screen_gene_extend,sim_gene)
    
    similar_gene[i]<-paste0(sim_gene$gene_name,collapse = ";")
    similar_prot[i]<-paste0(sim_prot,collapse = ";")
    
    
    
  }
  print(sum(extend_n))
  
  training_prot<-data.frame(training_prot)
  training_prot$extend_n<-extend_n
  training_prot$similar_prot<-similar_prot
  training_prot$similar_gene<-similar_gene
  
  prot_scope_comb<-rbind(prot_scope_comb,unique(screen_gene_extend))
  
  write.table(training_prot,
              paste0(data_folder,"prot_similar_",file_list[p],".csv"),append = FALSE, quote = FALSE, sep = ",",row.names = FALSE, col.names = TRUE)
  
  write.table(prot_scope_comb,
              paste0(data_folder,"Prot_",file_list[p],"_scope.txt"),append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
  
}
##############

abem_prot<-read.table(paste0(data_folder,"Prot_abem_scope.txt"),header=T,quote="",sep="\t",fill=TRUE)
bini_prot<-read.table(paste0(data_folder,"Prot_bini_scope.txt"),header=T,quote="",sep="\t",fill=TRUE)
olap_prot<-read.table(paste0(data_folder,"Prot_olap_scope.txt"),header=T,quote="",sep="\t",fill=TRUE)

abem_prot<-unique(abem_prot)
bini_prot<-unique(bini_prot)
olap_prot<-unique(olap_prot)

tmp<-table(c(abem_prot$uniprot_id,bini_prot$uniprot_id,olap_prot$uniprot_id))

common_prot<-tmp[which(tmp==3)]

length(common_prot)

write.table(common_prot,
            paste0(data_folder,"Prot_scope_common.csv"),
            append = FALSE, quote = F, sep = ",",row.names = F, col.names = F)

venn.plot <- venn.diagram(  
  x = list(  
    Abem = abem_prot$uniprot_id,  
    Bini = bini_prot$uniprot_id,  
    Olap = olap_prot$uniprot_id  
  ),  
  category.names = c("Abem", "Bini", "Olap"),  
  filename = NULL,  
  output = TRUE,  
  imagetype = "pdf",  
  height = 480,  
  width = 480,  
  resolution = 300,  
  compression = "lzw",  
  units = "px",  
  lwd = 2,  
  cat.cex = 2,  
  cat.fontface = "bold",  
  cat.dist = c(0.09, 0.09, 0.09),  
  label.col = "black",  
  cex = 1.5,  
  imaginary.rect = FALSE  
)  

pdf(file=paste0(figure_folder,"Prot_scope_venn_diagram.pdf"), width=6, height=6)
grid.draw(venn.plot)
dev.off()



