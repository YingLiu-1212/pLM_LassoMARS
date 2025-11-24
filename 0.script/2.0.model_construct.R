# This R script implements the construction and evaluation of machine learning models for predicting drug response to three cancer therapeutics (CDK inhibitors, MEK inhibitors, and PARP inhibitors). The code employs an ensemble modeling approach combining LASSO feature selection with Multivariate Adaptive Regression Splines (MARS) to build predictive models based on protein and amino acid embeddings integrated with structural and biophysical features. Key steps include: loading pre-processed training data containing protein embeddings, structural topology features, and drug response phenotypes; training separate models for each drug class using the Lasso-MARS pipeline; generating model visualization plots and exporting coefficient tables for interpretability analysis; performing validation on held-out test sets to assess prediction accuracy; and conducting comprehensive performance evaluation using ROC analysis with AUC calculations. The script produces publication-ready visualizations including model component plots and combined ROC curves, enabling systematic comparison of model performance across different drug classes and providing validated predictive tools for drug response prediction in cancer research.

library(RColorBrewer)
library(stringr)
# options(scipen=999)
options(stringsAsFactors = FALSE)

source("~/Other_project/CRC/reproducibility/script/function_pLM_LassoMARS.R")
setwd("~/Other_project/CRC/reproducibility/c.training")

output_folder<-"../o.output_figures/"

abem_info<-read.table(paste0("mars_train_abem_info.txt"),quote="",header=T,sep="\t",fill=TRUE)
abem_eb<-read.table(paste0("mars_train_abem_eb_all.txt"),quote="",header=T,sep="\t",fill=TRUE)

bini_info<-read.table(paste0("mars_train_bini_info.txt"),quote="",header=T,sep="\t",fill=TRUE)
bini_eb<-read.table(paste0("mars_train_bini_eb_all.txt"),quote="",header=T,sep="\t",fill=TRUE)

olap_info<-read.table(paste0("mars_train_olap_info.txt"),quote="",header=T,sep="\t",fill=TRUE)
olap_eb<-read.table(paste0("mars_train_olap_eb_all.txt"),quote="",header=T,sep="\t",fill=TRUE)


save.image(file = "./model_train_data.RData")

#########

load("./model_train_data.RData")

model_name<-"pLM_LassoMARS"

identical(abem_eb[,1],abem_info$id)

identical(bini_eb[,1],bini_info$id)

CDKi_model<-model_train(abem_eb,abem_info,label="CDKi")

MEKi_model<-model_train(bini_eb,bini_info,label="MEKi")

PARPi_model<-model_train(olap_eb,olap_info,label="PARPi")


######

pdf(paste0(output_folder,"plot_model_CDK.pdf"))
plot(CDKi_model[["model"]])
dev.off()
write.table(CDKi_model[["model"]]$coefficients,
            paste0("model_CDK_coefficients.csv"),
            append = FALSE, quote = F, sep = ",",row.names = T, col.names = NA)


pdf(paste0(output_folder,"plot_model_MEK.pdf"))
plot(MEKi_model[["model"]])
dev.off()
write.table(MEKi_model[["model"]]$coefficients,
            paste0("model_MEK_coefficients.csv"),
            append = FALSE, quote = F, sep = ",",row.names = T, col.names = NA)

pdf(paste0(output_folder,"plot_model_PARP.pdf"))
plot(PARPi_model[["model"]])
dev.off()
write.table(PARPi_model[["model"]]$coefficients,
            paste0("model_PARP_coefficients.csv"),
            append = FALSE, quote = F, sep = ",",row.names = T, col.names = NA)
# 
#########
save(CDKi_model,MEKi_model,PARPi_model,file = paste0(model_name,".RData"))

#######

CDKi_validation<-vali_pred(CDKi_model[["model"]],abem_eb[,CDKi_model[["select"]]],abem_info)

MEKi_validation<-vali_pred(MEKi_model[["model"]],bini_eb[,MEKi_model[["select"]]],bini_info)

PARPi_validation<-vali_pred(PARPi_model[["model"]],olap_eb[,PARPi_model[["select"]]],olap_info)

CDKi_validation<-vali_site_proc(CDKi_validation)
MEKi_validation<-vali_site_proc(MEKi_validation)
PARPi_validation<-vali_site_proc(PARPi_validation)


roc_PARPi <- roc(response = PARPi_validation$category,predictor = PARPi_validation$prediction)
roc_CDKi <- roc(response = CDKi_validation$category,predictor = CDKi_validation$prediction)
roc_MEKi <- roc(response = MEKi_validation$category,predictor = MEKi_validation$prediction)

auc_PARPi <- round(roc_PARPi$auc,2)
auc_CDKi <- round(roc_CDKi$auc,2)
auc_MEKi <- round(roc_MEKi$auc,2)

color_set<-brewer.pal(3, 'Dark2')

pdf(file=paste0(output_folder,"plot_ROC_vali_comb.pdf"), width=5, height=5)

plot_roc_PARPi<-plot(roc_PARPi, col=color_set[1], main="ROC of validate datasets",legacy.axes=T) 
plot_roc_CDKi<-plot(roc_CDKi, add=TRUE, col=color_set[2])
plot_roc_MEKi<-plot(roc_MEKi, add=TRUE, col=color_set[3])

legend("bottomright",box.lwd = 0,
       legend=c(paste0("PARPi AUC=",round(plot_roc_PARPi$auc,2)),paste0("CDKi AUC=",round(plot_roc_CDKi$auc,2)),paste0("MEKi AUC=",round(plot_roc_MEKi$auc,2))),cex = 0.8,
       col=color_set,lty=1)

dev.off()


