## 5.0.Unires_nb_PDB_dist_proc.R
library(matrixStats)
library(dplyr)
library(bio3d)
library(stringr)
options(scipen=999)
options(stringsAsFactors = FALSE)
setwd("~/Other_project/CRC/reproducibility/f.phos_dist_effect")
source("~/Other_project/CRC/reproducibility/script/_forAnalysis.R")


data_folder<-"../a.data/"
model_folder<-"../c.training/"


prot_gene<-read.table(paste0(data_folder,"UP000005640_9606_prot_gene.txt"),quote="",header=T,sep="\t",fill=TRUE)

PTM<-read.table(paste0(data_folder,"uniprotkb_PDB_UniRes.txt"),quote="\"",header=T,sep="\t",fill=TRUE)
PTM<-PTM[,c("Uniprot","pdb_chain","chain_len","id","Pos","AA","Type","PTM")]

PTM<-PTM[which(PTM$Uniprot%in%prot_gene$uniprot_id),]
PTM$label<-paste0(PTM$Type,"-",PTM$id)
PTM<-unique(PTM)
PTM<-PTM[which(PTM$Type!="N-acetylalanine"),]

##
PDB_folder<-"~/Published_data/PDB/"
file_path<-read.table(paste0(PDB_folder,"PDB_file_path_all.txt"),quote="",header=F,sep=",",fill=TRUE)
file_path<-file_path[,c(1,2,6)]
colnames(file_path)<-c("folder","Uniprot","pdb_id")
file_path<-unique(file_path)

protein_pdb<-PTM[,c("Uniprot","pdb_chain","chain_len")]
protein_pdb<-protein_pdb[which(protein_pdb$chain_len>50),]

protein_pdb<-unique(protein_pdb)
protein_pdb$pdb_id<-str_split_fixed(protein_pdb$pdb_chain,"_",n=2)[,1]
protein_pdb_path<-left_join(protein_pdb,file_path,by=c("Uniprot","pdb_id"))
########
PTM_phos<-PTM[which(PTM$Type%in%c("Phosphoserine","Phosphothreonine","Phosphotyrosine")),]
PTM_lysine<-PTM[which(PTM$Type%in%c("N6-acetyllysine","N6-succinyllysine")),]

Phosphoserine<-PTM[which(PTM$Type=="Phosphoserine"),]
Phosphothreonine<-PTM[which(PTM$Type=="Phosphothreonine"),]
Phosphotyrosine<-PTM[which(PTM$Type=="Phosphotyrosine"),]

succinyllysine<-PTM[which(PTM$Type=="N6-succinyllysine"),]
acetyllysine<-PTM[which(PTM$Type=="N6-acetyllysine"),]
######

prot_list<-unique(PTM[,c("Uniprot","pdb_chain")])
#which(prot_list=="O14497")

pdb_folder<-"~/Published_data/PDB/pdb_file/"

neighbour_PTM<-NULL

## nrow(prot_list)2496
# split<-5

for(split in c(1:5)){
  print(split)
  subsets<-c((500*(split-1)+1):min(nrow(prot_list),(500*split)))
  
  for(i in subsets){
    if(i%%100==0)print(i)
    
    chain<-prot_list[i,2]
    pid<-prot_list[i,1]
    
    pdb<-str_split_fixed(chain,"_",2)[1]
    chain_id<-str_split_fixed(chain,"_",2)[2]
    
    file_index<-which(protein_pdb_path$pdb_chain==chain & protein_pdb_path$Uniprot==pid)
    file_path<-paste0(protein_pdb_path$folder[file_index],"/")
    
    pdb_name<-paste0(pdb_folder,file_path,pdb,".pdb")
    
    if(file.exists(pdb_name)){
      
      pdb <- read.pdb(pdb_name)
      
      ca.inds <- atom.select(pdb, "calpha",chain=chain_id)
      if(length(ca.inds$atom)<=10){next}
      
      distance<-dm(pdb, inds = ca.inds,mask.lower=FALSE)
      diag(distance)<-0
      select_atom<-pdb$atom[ca.inds$atom, ]
      
      resid<-select_atom$resid
      resno<-select_atom$resno
      resid_no<-paste0(resno,aa321(resid))
      
      colnames(distance)<-rownames(distance)<-resid_no
      
      inter_contact<-round(distance,4)
      
      ###
      # inter_contact<-read.table(pdb_name,quote="",header=T,sep="\t",fill=TRUE)
      # id<-inter_contact[,1]
      # rownames(inter_contact)<-id
      # inter_contact<-inter_contact[,-1]
      # colnames(inter_contact)<-id
      
      id<-resid_no
      
      neighbour<-data.frame(uniprot_id=pid,pdb_chain=chain,id=id)
      
      ptm_site<-PTM[which(PTM$pdb_chain==chain & PTM$id%in%id),]
      
      ptm_around<-dist_aa_site(inter_contact,ptm_site)
      neighbour$PTMDist<-ptm_around$AADist
      neighbour$PTMSite<-ptm_around$AASite
      
      phos_site<-PTM_phos[which(PTM_phos$pdb_chain==chain& PTM_phos$id%in%id),]
      
      phos_around<-dist_aa_site(inter_contact,phos_site)
      neighbour$PhosDist<-phos_around$AADist
      neighbour$PhosSite<-phos_around$AASite
      
      
      #######
      lysine_site<-PTM_lysine[which(PTM_lysine$pdb_chain==chain& PTM_lysine$id%in%id),]
      
      lysine_around<-dist_aa_site(inter_contact,lysine_site)
      
      neighbour$lysineDist<-lysine_around$AADist
      neighbour$lysineSite<-lysine_around$AASite
      
      #####
      
      Phosphoserine_site<-Phosphoserine[which(Phosphoserine$pdb_chain==chain& Phosphoserine$id%in%id),]
      
      Phosphoserine_around<-dist_aa_site(inter_contact,Phosphoserine_site)
      
      neighbour$PhosphoserineDist<-Phosphoserine_around$AADist
      neighbour$PhosphoserineSite<-Phosphoserine_around$AASite
      
      
      ######
      Phosphothreonine_site<-Phosphothreonine[which(Phosphothreonine$pdb_chain==chain& Phosphothreonine$id%in%id),]
      
      Phosphothreonine_around<-dist_aa_site(inter_contact,Phosphothreonine_site)
      
      neighbour$PhosphothreonineDist<-Phosphothreonine_around$AADist
      neighbour$PhosphothreonineSite<-Phosphothreonine_around$AASite
      
      
      #########
      Phosphotyrosine_site<-Phosphotyrosine[which(Phosphotyrosine$pdb_chain==chain& Phosphotyrosine$id%in%id),]
      
      Phosphotyrosine_around<-dist_aa_site(inter_contact,Phosphotyrosine_site)
      
      neighbour$PhosphotyrosineDist<-Phosphotyrosine_around$AADist
      neighbour$PhosphotyrosineSite<-Phosphotyrosine_around$AASite
      
      
      ####
      succinyllysine_site<-succinyllysine[which(succinyllysine$pdb_chain==chain& succinyllysine$id%in%id),]
      
      succinyllysine_around<-dist_aa_site(inter_contact,succinyllysine_site)
      
      neighbour$succinyllysineDist<-succinyllysine_around$AADist
      neighbour$succinyllysineSite<-succinyllysine_around$AASite
      
      
      #######
      acetyllysine_site<-acetyllysine[which(acetyllysine$pdb_chain==chain& acetyllysine$id%in%id),]
      
      acetyllysine_around<-dist_aa_site(inter_contact,acetyllysine_site)
      
      neighbour$acetyllysineDist<-acetyllysine_around$AADist
      neighbour$acetyllysineSite<-acetyllysine_around$AASite
      
      
      neighbour_PTM<-rbind(neighbour_PTM,neighbour)
      
    } ## if(file.exists(pdb_name))
    
  } ## for(i in 1:length(prot_list))
  
  write.table(neighbour_PTM,
              paste0("Unires_AA_PDB_neighbour_",split,".txt"),
              append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)
  
  
}

######

pdb_aa_all<-NULL
for(split in 1:5){
  pdb_aa<-read.table(paste0("Unires_AA_PDB_neighbour_",split,".txt"),quote="",header=T,sep="\t",fill=TRUE)
  pdb_aa_all<-rbind(pdb_aa_all,pdb_aa)
}

write.table(pdb_aa_all,
            paste0("Unires_AA_PDB_neighbour.txt"),
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)

#####
# PTM_nb<-pdb_aa_all[which(pdb_aa_all$PTMSite!="-"),1:5]
# 
# PTM_nb$PTMDist<-as.numeric(PTM_nb$PTMDist)
# 
# PTM_nb<-PTM_nb[order(PTM_nb$PTMDist),]
# 
# PTM_nb<-PTM_nb[!duplicated(PTM_nb[,c("uniprot_id","id")]),]
# 
# write.table(PTM_nb,
#             paste0("UniresPTM_AA_PDB_neighbour_dedup.txt"),
#             append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)


#####
Phos_nb<-pdb_aa_all[which(pdb_aa_all$PhosSite!="-"),c("uniprot_id","pdb_chain","id","PhosDist","PhosSite")]

Phos_nb$PhosDist<-as.numeric(Phos_nb$PhosDist)

Phos_nb<-Phos_nb[order(Phos_nb$PhosDist),]

Phos_nb<-Phos_nb[!duplicated(Phos_nb[,c("uniprot_id","id")]),]

write.table(Phos_nb,
            paste0("UniresPhos_AA_PDB_neighbour_dedup.txt"),
            append = FALSE, quote = F, sep = "\t",row.names = F, col.names = T)

save(Phos_nb,file=paste0("Unires_Phos_AA_PDB_neighbour.RData"))


