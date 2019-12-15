#GNU GENERAL PUBLIC LICENSE - Version 3
#
#Affy_gene_expression_Preprocess_and_analysis.r Copyright (C) 2019
#
# This program is free software: you can redistribute it and/or modify
#
# it under the terms of the GNU General Public License as published by
#
# the Free Software Foundation, either version 3 of the License, or
#
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
#
# but WITHOUT ANY WARRANTY; without even the implied warranty of
#
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#
# GNU General Public License for more details.
#
# 
#
# You should have received a copy of the GNU General Public License
#
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact information: augustoanguitaruiz@gmail.com
#
#Anguita-Ruiz A, Segura-Delgado A, Alcala R & Alcala-Fdez J



















# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #
# PART 1
#These lines are instructions for running the "./GeneSeqRules/required_codes/Part_1_Affy_gene_expression_Preprocessing_and_analysis.r" R script, which is intended for the pre-processing and analysis of raw fluorescence (.cel) files. This script performs:
#- The merging of individual .cel files.
#- The visual inspection of raw data and the generation of quality control plots.
#- The RMA normalization process.
#- The annotation process.
#- The differentially expressed gene analysis.
# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #






# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #  Loading and preparing the dataset
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #




# Loading required packages. If you miss some of them, please, install it using the install.packages("BiocManager") and then the BiocManager::install() function:

require(Matrix)
require(lattice)
require(fdrtool)
require(rpart)
require(oligo)
require(ArrayExpress)
require(limma)
require(Biobase)
require(Biostrings)
require(genefilter)
require(ggplot2)
library(magrittr)
library(dplyr)
require(pd.hugene.1.1.st.v1)
library(pca3d)



# Set working directory

setwd("./GeneSeqRules")



# Load input fluorescence raw .cel files.

celpath = "./Data/raw_cel_files/"
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list)



# Add annotations to loaded cel data. Name specifications and substitutions should be modified according to the study design.

ph = data@phenoData
ph
ph@data
rownames(ph@data) <- sub(".*?_(.+)", "\\1", rownames(ph@data))
rownames(ph@data) <- gsub("\\..*","",rownames(ph@data))
rownames(ph@data) <- gsub("Losers","Losers_",rownames(ph@data))
rownames(ph@data) <- gsub("Regainers","Regainers_",rownames(ph@data))
rownames(ph@data) 
Experimental_group <- gsub("\\_.*","",rownames(ph@data))
Timerecord <- gsub(".*_","",rownames(ph@data)) 
Experimental_group_and_timerecord  <- paste(Experimental_group,Timerecord)
Experimental_group_and_timerecord  <- gsub(" ","_",Experimental_group_and_timerecord)
rownames(ph@data)
Experimental_group
Timerecord
Experimental_group_and_timerecord
ph@data[ ,1] = rownames(ph@data) 
nindinviduals <- length(rownames(ph@data))



# Retrieve intensity info

head(oligo::intensity(data))



# Retrieve annotation info

feat = data@featureData
feat@data



# Retrieve probe set IDs

featureNames(data)




# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #  Quality control plots and Normalization
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #




# Plot all sample histograms together

pmexp = pm(data)
sampleNames = vector()
logs = vector()
for (i in 1:57)
{
sampleNames = c(sampleNames,rep(ph@data[i,1],dim(pmexp)[1]))
logs = c(logs,log2(pmexp[,i]))
}
logData = data.frame(logInt=logs,sampleName=sampleNames)
dataHist2 = ggplot(logData, aes(logInt, colour = sampleName)) 
dataHist2 + geom_density()
ggsave("./Output_Plots/Pre_processing/histograms.pdf") 



# Plot each sample histogram separately

for (i in 1:nindinviduals)
{
name = paste("histogram",i,".pdf",sep="")
pdf(paste("./Output_Plots/Pre_processing/",name,sep=""))
hist(data[,i],lwd=2,which='all',ylab='Density',xlab='Log2 intensities',main=ph@data$sample[i])
dev.off()
}



# Plot each sample boxplot separately

name = "boxplot.pdf"
pdf(paste("./Output_Plots/Pre_processing/",name,sep=""))
boxplot(data,which='all',col='red',names=ph@data$sample) 
dev.off()



# Plot each sample MA plot separately

for (i in 1:nindinviduals)
{
name = paste("MAplot",i,".jpg",sep="")
jpeg(paste("./Output_Plots/Pre_processing/MA",name,sep=""))
MAplot(data,which=i)
dev.off()
}



# Perform RMA Data Normalization

data.rma = rma(data)
data.matrix = exprs(data.rma)



# Plot each sample normalized MA plot separately

for (i in 1:nindinviduals)
{
name = paste("MAplotnorm",ph@data[i,1],".jpg",sep="")
jpeg(paste("./Output_Plots/Pre_processing/",name,sep=''))
MAplot(data.rma,which=i)
dev.off()
}



# Plot sample normalized boxplot

name = "./Output_Plots/Pre_processing/boxplotnorm.pdf"
pdf(name)
boxplot(data.matrix,col='red')
dev.off()



# Plot PCA

a <- t(data.matrix)
a <- a[ , apply(a, 2, var) != 0]
data.PC = prcomp(a,scale.=TRUE)
factor_group <- gsub("WeightRegainers_T0","T0",Experimental_group_and_timerecord)
factor_group <- as.factor(gsub("WeightLosers_T0","T0",factor_group))
pca2d(data.PC,group=factor_group)
pca3d(data.PC,group=factor_group,legend="top")




# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #  Array annotation
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #




#Load the specific array package.

require("hgu133plus2.db") 
require("pd.hg.u133.plus.2")



#Check if studied probes are available within the loaded package.

data.rma = rma(data)
data.matrix = exprs(data.rma)
(rownames(data.rma)) %in% keys(hgu133plus2.db) %>% summary
#   Mode    TRUE 
#logical   54675 



columns(hgu133plus2.db)
ae.annots <- AnnotationDbi::select(x= hgu133plus2.db,
  keys    = rownames(data.rma),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL","ONTOLOGY"),
  keytype = "PROBEID"  )
dim(ae.annots)
all(rownames(data.rma) %in% ae.annots$PROBEID)
dup.ids <- ae.annots$PROBEID[duplicated(ae.annots$PROBEID)] %>% 
  unique %>%
  sort
ae.annots[ae.annots$PROBEID == dup.ids[1], ]
collapser <- function(x){
  x %>% unique %>% sort %>% paste(collapse = "|")
  }



# Redefinition of ae.annots

ae.annots <- AnnotationDbi::select(
  x       = hgu133plus2.db,
  keys    = rownames(data.rma),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL","ONTOLOGY"),
  keytype = "PROBEID"
  ) %>%
  group_by(PROBEID) %>%
  summarise_all(funs(collapser)) %>%
  ungroup
ae.annots %>% tail
all(ae.annots$PROBEID == rownames(data.rma))
 fd <- new("AnnotatedDataFrame",
  data = data.frame(ae.annots[, -1], stringsAsFactors = FALSE)
  )
rownames(fd) <- ae.annots$PROBEID
featureData(data.rma) <- fd
anotacion <- featureData(data.rma)
anotacion <- anotacion@data
write.table(x=anotacion, file="./Output_files/Probe_annotations/PROBE_ANNOT_hugeneu133ttranscriptcluster",sep="\t")
data.matrix = exprs(data.rma)



# Save output file

data.matrix2 <-  data.matrix
colnames(data.matrix2) <- rownames(ph@data)
rownames(data.matrix2) <- paste(rownames(data.matrix2),anotacion[,3],sep="/")
data.matrix2 <- as.data.frame(data.matrix2)
data.matrix2 <- data.matrix2[,order(colnames(data.matrix2))]
colnames(data.matrix2)
write.csv2(x=data.matrix2, file="./Output_files/Normalized_Array/NORM_DATA_ALL_ARRAYS_HUMAN_GSE103766.csv")




# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #  Differentially Expressed (DE) Analysis
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #




# Each experimental group will be assesed separately

ph@data[ ,2] = Experimental_group_and_timerecord
colnames(ph@data)[2]="Experimental_GROUP"
ph@data
groupsG = ph@data$Experimental_GROUP
fG = as.factor(groupsG)
patient <- gsub("WeightRegainers_","WeightRegainers.",rownames(ph@data))
patient <- gsub("WeightLosers_","WeightLosers.",patient)
patient  <- gsub("\\_.*","",patient )
patient  <- gsub("ers.","ers_",patient)
patient
ph@data[ ,3] = patient
colnames(ph@data)[3]="Patient"
ph@data
groupsP = as.factor(ph@data$Patient)
fP = groupsP



###################
############ Experimental Group 1: "WeightLosers"
##################

WeightLosers__NORM <- data.rma[,which((ph@data[,2] %in% c("WeightLosers_T0","WeightLosers_T1","WeightLosers_T2")))]
rownames(WeightLosers__NORM@phenoData@data) <- sub(".*?_(.+)", "\\1", rownames(WeightLosers__NORM@phenoData@data))
rownames(WeightLosers__NORM@phenoData@data) <- gsub("\\..*","",rownames(WeightLosers__NORM@phenoData@data))
rownames(WeightLosers__NORM@phenoData@data) <- gsub("Losers","Losers_",rownames(WeightLosers__NORM@phenoData@data))
rownames(WeightLosers__NORM@phenoData@data) <- gsub("Regainers","Regainers_",rownames(WeightLosers__NORM@phenoData@data))
rownames(WeightLosers__NORM@phenoData@data) 
Experimental_group_wr <- gsub("\\_.*","",rownames(WeightLosers__NORM@phenoData@data))
Timerecord_wr <- gsub(".*_","",rownames(WeightLosers__NORM@phenoData@data)) 
Experimental_group_and_timerecord_wr  <- paste(Experimental_group_wr,Timerecord_wr)
Experimental_group_and_timerecord_wr  <- gsub(" ","_",Experimental_group_and_timerecord_wr)
patient_wr <- gsub("WeightRegainers_","WeightRegainers.",
rownames(WeightLosers__NORM@phenoData@data) )
patient_wr <- gsub("WeightLosers_","WeightLosers.",patient_wr)
patient_wr  <- gsub("\\_.*","",patient_wr)
patient_wr  <- gsub("ers.","ers_",patient_wr)
patient_wr
Experimental_GROUP = Timerecord_wr
groupsG = as.factor(Experimental_GROUP)
fG = groupsG
inds = patient_wr
groupsP = as.factor(inds)
fP = groupsP
WeightLosers__NORM@phenoData@data[,1] <- fG
WeightLosers__NORM@phenoData@data[,2] <- fP
colnames(WeightLosers__NORM@phenoData@data) <- c("Experimental_Group","Patient")
WeightLosers__NORM@phenoData@data
table(WeightLosers__NORM@phenoData@data$Experimental_Group)



# Perform all time comparisons (Baseline vs T1, T1 vs T2 and T2 vs Baseline.

paired.design = model.matrix(~ fP + fG)
data.fit = lmFit(WeightLosers__NORM,paired.design)
TODOS_VS_VLCDt0 = eBayes(data.fit)
head(data.fit$coefficients)
head(TODOS_VS_VLCDt0$p.value)



# Generate volcano plots for DE genes in each time interval. 



#WS-VS-START

name = "./Output_Plots/Volcano_Plots/Group_1/Volcano_end_vs_START.pdf"
pdf(name)
volcanoplot(TODOS_VS_VLCDt0,coef=8,highlight=10)
dev.off()



#WL-VS-START

name = "./Output_Plots/Volcano_Plots/Group_1/Volcano_WL_vs_START.pdf"
pdf(name)
volcanoplot(TODOS_VS_VLCDt0,coef=7,highlight=10)
dev.off()



#WS-vs-WL 

cont.matrix = makeContrasts(fGT2-fGT1,levels=paired.design)
colnames(data.fit$coefficients)[1] <- "Intercept"
data.contr = contrasts.fit(data.fit,cont.matrix)
WS_vs_WL_VLCD = eBayes(data.contr)
head(data.contr$coefficients)
head(WS_vs_WL_VLCD$p.value)
name = "./Output_Plots/Volcano_Plots/Group_1/Volcano_end_vs_WL.pdf"
pdf(name)
volcanoplot(WS_vs_WL_VLCD,coef=1,highlight=10)
dev.off()



# Filter genes by Pvalue and show the name of loci and probes.



#Names WL-VS-START

WL_VS_START_VLCD <- rownames(TODOS_VS_VLCDt0[which(TODOS_VS_VLCDt0$p.value[,7] < 0.01),])
WL_VS_START_VLCD  <- anotacion[rownames(anotacion) %in% WL_VS_START_VLCD,]
WL_VS_START_VLCD 



#Names WS-VS-START

WS_VS_START_VLCD <- rownames(TODOS_VS_VLCDt0[which(TODOS_VS_VLCDt0$p.value[,8] < 0.01),])
WS_VS_START_VLCD  <- anotacion[rownames(anotacion) %in% WS_VS_START_VLCD,]
WS_VS_START_VLCD 



#Names WS-VS-WL

WS_VS_WL_VLCD_ <- names(WS_vs_WL_VLCD$p.value[which(WS_vs_WL_VLCD$p.value < 0.01),])
WS_VS_WL_VLCD_ <- anotacion[rownames(anotacion) %in% WS_VS_WL_VLCD_,]
WS_VS_WL_VLCD_



# SELECT D.E. MULTIPLE TESTING ADJUSTED GENES.



#WL-VS-START

options(digits=2)
tab_WL_VS_START_VLCD = topTable(TODOS_VS_VLCDt0,coef=7,number=500,adjust.method="bonferroni")
topgenes = tab_WL_VS_START_VLCD[tab_WL_VS_START_VLCD[, "adj.P.Val"] < 0.05, ]
dim(topgenes)
topups_WL_VS_START_VLCD = topgenes[topgenes[, "logFC"] >= 0.5, ]
dim(topups_WL_VS_START_VLCD)
topdowns_WL_VS_START_VLCD = topgenes[topgenes[, "logFC"] <= -0.5, ]
dim(topdowns_WL_VS_START_VLCD)



#WS-VS-START

options(digits=2)
tab_WS_VS_START_VLCD = topTable(TODOS_VS_VLCDt0,coef=8,number=500,adjust.method="bonferroni")
topgenes = tab_WS_VS_START_VLCD[tab_WS_VS_START_VLCD[, "adj.P.Val"] < 0.05, ]
dim(topgenes)
topups_WS_VS_START_VLCD = topgenes[topgenes[, "logFC"] >= 0.5, ]
dim(topups_WS_VS_START_VLCD)
topdowns_WS_VS_START_VLCD = topgenes[topgenes[, "logFC"] <= -0.5, ]
dim(topdowns_WS_VS_START_VLCD)



#WS-VS-WL

options(digits=2)
tab_WS_VS_WL_VLCD = topTable(WS_vs_WL_VLCD,coef=1,number=500,adjust.method="bonferroni")
topgenes = tab_WS_VS_WL_VLCD[tab_WS_VS_WL_VLCD[, "adj.P.Val"] < 0.05, ]
dim(topgenes)
topups_WS_VS_WL_VLCD = topgenes[topgenes[, "logFC"] >= 0.5, ]
dim(topups_WS_VS_WL_VLCD)
topdowns_WS_VS_WL_VLCD = topgenes[topgenes[, "logFC"] <= -0.5, ]
dim(topdowns_WS_VS_WL_VLCD)




# SAVE LISTS FOR D.E. MULTIPLE TESTING ADJUSTED GENES.



#WL-VS-START

options(digits=2)
tab_WL_VS_START_VLCD = topTable(TODOS_VS_VLCDt0,coef=7,number=1000,adjust.method="bonferroni")
topgenes = tab_WL_VS_START_VLCD[tab_WL_VS_START_VLCD[, "P.Value"] < 0.0001, ]
dim(topgenes)
write.csv2(topgenes,"./Output_files/DE_genes/Group_1/WL_VS_START_wl_PVALUE_0_0001.csv")
rownames(topgenes)
WL_VS_START_VLCD <- anotacion[rownames(anotacion) %in% rownames(topgenes),]
names_genes_WL_VS_START_VLCD <- paste(rownames(WL_VS_START_VLCD),WL_VS_START_VLCD[,3],sep="/")
write.csv2(names_genes_WL_VS_START_VLCD,"./Output_files/DE_genes/Group_1/WL_VS_START_wl_PVALUE_0_0001_NOMBRES.csv")



#WS-VS-START

options(digits=2)
tab_WS_VS_START_VLCD = topTable(TODOS_VS_VLCDt0,coef=8,number=1000,adjust.method="bonferroni")
topgenes = tab_WS_VS_START_VLCD[tab_WS_VS_START_VLCD[, "P.Value"] < 0.00001, ]
dim(topgenes)
write.csv2(topgenes,"./Output_files/DE_genes/Group_1/END_VS_START_wl_PVALUE_0_00001.csv")
rownames(topgenes)
WS_VS_START_VLCD <- anotacion[rownames(anotacion) %in% rownames(topgenes),]
names_genes_WS_VS_START_VLCD <- paste(rownames(WS_VS_START_VLCD),WS_VS_START_VLCD[,3],sep="/")
write.csv2(names_genes_WS_VS_START_VLCD,"./Output_files/DE_genes/Group_1/END_VS_START_wl_PVALUE_0_00001_NOMBRES.csv")



#WS-VS-WL

options(digits=2)
tab_WS_VS_WL_VLCD = topTable(WS_vs_WL_VLCD,coef=1,number=1000,adjust.method="bonferroni")
topgenes = tab_WS_VS_WL_VLCD[tab_WS_VS_WL_VLCD[, "P.Value"] < 0.001, ]
dim(topgenes)
write.csv2(topgenes,"./Output_files/DE_genes/Group_1/end_VS_WL_wl_PVALUE_0_001.csv")
rownames(topgenes)
WS_VS_WL_VLCD_ <- anotacion[rownames(anotacion) %in% rownames(topgenes),]
names_genes_WS_VS_WL_VLCD_ <- paste(rownames(WS_VS_WL_VLCD_),WS_VS_WL_VLCD_[,3],sep="/")
write.csv2(names_genes_WS_VS_WL_VLCD_,"./Output_files/DE_genes/Group_1/end_VS_WL_wl_PVALUE_0_001_nombres.csv")



###################
############ Experimental Group 2: "WeightRegainers"
##################

WeightLosers__NORM <- data.rma[,which((ph@data[,2] %in% c("WeightRegainers_T0","WeightRegainers_T1","WeightRegainers_T2")))]
rownames(WeightLosers__NORM@phenoData@data) <- sub(".*?_(.+)", "\\1", rownames(WeightLosers__NORM@phenoData@data))
rownames(WeightLosers__NORM@phenoData@data) <- gsub("\\..*","",rownames(WeightLosers__NORM@phenoData@data))
rownames(WeightLosers__NORM@phenoData@data) <- gsub("Losers","Losers_",rownames(WeightLosers__NORM@phenoData@data))
rownames(WeightLosers__NORM@phenoData@data) <- gsub("Regainers","Regainers_",rownames(WeightLosers__NORM@phenoData@data))
rownames(WeightLosers__NORM@phenoData@data) 
Experimental_group_wr <- gsub("\\_.*","",rownames(WeightLosers__NORM@phenoData@data))
Timerecord_wr <- gsub(".*_","",rownames(WeightLosers__NORM@phenoData@data)) 
Experimental_group_and_timerecord_wr  <- paste(Experimental_group_wr,Timerecord_wr)
Experimental_group_and_timerecord_wr  <- gsub(" ","_",Experimental_group_and_timerecord_wr)
patient_wr <- gsub("WeightRegainers_","WeightRegainers.",
rownames(WeightLosers__NORM@phenoData@data) )
patient_wr <- gsub("WeightLosers_","WeightLosers.",patient_wr)
patient_wr  <- gsub("\\_.*","",patient_wr)
patient_wr  <- gsub("ers.","ers_",patient_wr)
patient_wr
Experimental_GROUP = Timerecord_wr
groupsG = as.factor(Experimental_GROUP)
fG = groupsG
inds = patient_wr
groupsP = as.factor(inds)
fP = groupsP
WeightLosers__NORM@phenoData@data[,1] <- fG
WeightLosers__NORM@phenoData@data[,2] <- fP
colnames(WeightLosers__NORM@phenoData@data) <- c("Experimental_Group","Patient")
WeightLosers__NORM@phenoData@data
table(WeightLosers__NORM@phenoData@data$Experimental_Group)



# Perform all time comparisons (Baseline vs T1, T1 vs T2 and T2 vs Baseline.

paired.design = model.matrix(~ fP + fG)
data.fit = lmFit(WeightLosers__NORM,paired.design)
TODOS_VS_VLCDt0 = eBayes(data.fit)
head(data.fit$coefficients)
head(TODOS_VS_VLCDt0$p.value)



# Generate volcano plots for DE genes in each time interval. 



#WS-VS-START

name = "./Output_Plots/Volcano_Plots/Group_2/Volcano_end_vs_START.pdf"
pdf(name)
volcanoplot(TODOS_VS_VLCDt0,coef=8,highlight=10)
dev.off()



#WL-VS-START

name = "./Output_Plots/Volcano_Plots/Group_2/Volcano_WL_vs_START.pdf"
pdf(name)
volcanoplot(TODOS_VS_VLCDt0,coef=7,highlight=10)
dev.off()



#WS-vs-WL 

cont.matrix = makeContrasts(fGT2-fGT1,levels=paired.design)
colnames(data.fit$coefficients)[1] <- "Intercept"
data.contr = contrasts.fit(data.fit,cont.matrix)
WS_vs_WL_VLCD = eBayes(data.contr)
head(data.contr$coefficients)
head(WS_vs_WL_VLCD$p.value)
name = "./Output_Plots/Volcano_Plots/Group_2/Volcano_end_vs_WL.pdf"
pdf(name)
volcanoplot(WS_vs_WL_VLCD,coef=1,highlight=10)
dev.off()



# Filter genes by Pvalue and show the name of loci and probes.



#Names WL-VS-START

WL_VS_START_VLCD <- rownames(TODOS_VS_VLCDt0[which(TODOS_VS_VLCDt0$p.value[,7] < 0.01),])
WL_VS_START_VLCD  <- anotacion[rownames(anotacion) %in% WL_VS_START_VLCD,]
WL_VS_START_VLCD 



#Names WS-VS-START


WS_VS_START_VLCD <- rownames(TODOS_VS_VLCDt0[which(TODOS_VS_VLCDt0$p.value[,8] < 0.01),])
WS_VS_START_VLCD  <- anotacion[rownames(anotacion) %in% WS_VS_START_VLCD,]
WS_VS_START_VLCD 



#Names WS-VS-WL

WS_VS_WL_VLCD_ <- names(WS_vs_WL_VLCD$p.value[which(WS_vs_WL_VLCD$p.value < 0.01),])
WS_VS_WL_VLCD_ <- anotacion[rownames(anotacion) %in% WS_VS_WL_VLCD_,]
WS_VS_WL_VLCD_



# SELECT D.E. MULTIPLE TESTING ADJUSTED GENES.



#WL-VS-START

options(digits=2)
tab_WL_VS_START_VLCD = topTable(TODOS_VS_VLCDt0,coef=7,number=500,adjust.method="bonferroni")
topgenes = tab_WL_VS_START_VLCD[tab_WL_VS_START_VLCD[, "adj.P.Val"] < 0.05, ]
dim(topgenes)
topups_WL_VS_START_VLCD = topgenes[topgenes[, "logFC"] >= 0.5, ]
dim(topups_WL_VS_START_VLCD)
topdowns_WL_VS_START_VLCD = topgenes[topgenes[, "logFC"] <= -0.5, ]
dim(topdowns_WL_VS_START_VLCD)



#WS-VS-START

options(digits=2)
tab_WS_VS_START_VLCD = topTable(TODOS_VS_VLCDt0,coef=8,number=500,adjust.method="bonferroni")
topgenes = tab_WS_VS_START_VLCD[tab_WS_VS_START_VLCD[, "adj.P.Val"] < 0.05, ]
dim(topgenes)
topups_WS_VS_START_VLCD = topgenes[topgenes[, "logFC"] >= 0.5, ]
dim(topups_WS_VS_START_VLCD)
topdowns_WS_VS_START_VLCD = topgenes[topgenes[, "logFC"] <= -0.5, ]
dim(topdowns_WS_VS_START_VLCD)



#WS-VS-WL

options(digits=2)
tab_WS_VS_WL_VLCD = topTable(WS_vs_WL_VLCD,coef=1,number=500,adjust.method="bonferroni")
topgenes = tab_WS_VS_WL_VLCD[tab_WS_VS_WL_VLCD[, "adj.P.Val"] < 0.05, ]
dim(topgenes)
topups_WS_VS_WL_VLCD = topgenes[topgenes[, "logFC"] >= 0.5, ]
dim(topups_WS_VS_WL_VLCD)
topdowns_WS_VS_WL_VLCD = topgenes[topgenes[, "logFC"] <= -0.5, ]
dim(topdowns_WS_VS_WL_VLCD)




# SAVE LISTS FOR D.E. MULTIPLE TESTING ADJUSTED GENES.



#WL-VS-START

options(digits=2)
tab_WL_VS_START_VLCD = topTable(TODOS_VS_VLCDt0,coef=7,number=1000,adjust.method="bonferroni")
topgenes = tab_WL_VS_START_VLCD[tab_WL_VS_START_VLCD[, "P.Value"] < 0.0001, ]
dim(topgenes)
write.csv2(topgenes,"./Output_files/DE_genes/Group_2/WL_VS_START_wr_PVALUE_0_0001.csv")
rownames(topgenes)
WL_VS_START_VLCD <- anotacion[rownames(anotacion) %in% rownames(topgenes),]
names_genes_WL_VS_START_VLCD <- paste(rownames(WL_VS_START_VLCD),WL_VS_START_VLCD[,3],sep="/")
write.csv2(names_genes_WL_VS_START_VLCD,"./Output_files/DE_genes/Group_2/WL_VS_START_wR_PVALUE_0_0001_NOMBRES.csv")



#WS-VS-START

options(digits=2)
tab_WS_VS_START_VLCD = topTable(TODOS_VS_VLCDt0,coef=8,number=1000,adjust.method="bonferroni")
topgenes = tab_WS_VS_START_VLCD[tab_WS_VS_START_VLCD[, "P.Value"] < 0.00001, ]
dim(topgenes)
write.csv2(topgenes,"./Output_files/DE_genes/Group_2/END_VS_START_wR_PVALUE_0_00001.csv")
rownames(topgenes)
WS_VS_START_VLCD <- anotacion[rownames(anotacion) %in% rownames(topgenes),]
names_genes_WS_VS_START_VLCD <- paste(rownames(WS_VS_START_VLCD),WS_VS_START_VLCD[,3],sep="/")
write.csv2(names_genes_WS_VS_START_VLCD,"./Output_files/DE_genes/Group_2/END_VS_START_wR_PVALUE_0_00001_NOMBRES.csv")



#WS-VS-WL

options(digits=2)
tab_WS_VS_WL_VLCD = topTable(WS_vs_WL_VLCD,coef=1,number=1000,adjust.method="bonferroni")
topgenes = tab_WS_VS_WL_VLCD[tab_WS_VS_WL_VLCD[, "P.Value"] < 0.001, ]
dim(topgenes)
write.csv2(topgenes,"./Output_files/DE_genes/Group_2/end_VS_WL_wR_PVALUE_0_001.csv")
rownames(topgenes)
WS_VS_WL_VLCD_ <- anotacion[rownames(anotacion) %in% rownames(topgenes),]
names_genes_WS_VS_WL_VLCD_ <- paste(rownames(WS_VS_WL_VLCD_),WS_VS_WL_VLCD_[,3],sep="/")
write.csv2(names_genes_WS_VS_WL_VLCD_,"./Output_files/DE_genes/Group_2/end_VS_WL_wR_PVALUE_0_001_nombres.csv")




#................................. End of code





































# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #
# PART 2
#These lines are intended for the discretization of absolute gene expression values and the generation of the sequence database.
# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #





# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #  Loading input files and preparing working dataset
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #




# Loading required packages. If you miss some of them, please, install it using the install.packages() function:

library(stringr)
library(schoolmath)
library(plyr)
library(stringr)




# Load Normalized array data

BBDD_input <- read.csv2("./Output_files/Normalized_Array/NORM_DATA_ALL_ARRAYS_HUMAN_GSE103766.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1)



# Load lists and names of DE genes by experimental group



# For group_2

all_sig_WR_WL_START <- read.csv2("./Output_files/DE_genes/Group_2/WL_VS_START_wR_PVALUE_0_0001_NOMBRES.csv")
all_sig_WR_WL_START <- as.character(unique(all_sig_WR_WL_START[,2]))
all_sig_WR_WL_START
all_sig_WR_END_WL <- read.csv2("./Output_files/DE_genes/Group_2/end_VS_WL_wR_PVALUE_0_001_nombres.csv")
all_sig_WR_END_WL <- as.character(unique(all_sig_WR_END_WL[,2]))
all_sig_WR_END_WL



# For group_1

all_sig_Wl_WL_START <- read.csv2("./Output_files/DE_genes/Group_1/WL_VS_START_wl_PVALUE_0_0001_NOMBRES.csv")
all_sig_Wl_WL_START <- as.character(unique(all_sig_Wl_WL_START[,2]))
all_sig_Wl_WL_START
all_sig_Wl_END_WL <- read.csv2("./Output_files/DE_genes/Group_1/end_VS_WL_wl_PVALUE_0_001_nombres.csv")
all_sig_Wl_END_WL <- as.character(unique(all_sig_Wl_END_WL[,2]))
all_sig_Wl_END_WL



# Identify all significantly DE probes from loaded lists

all_sig <- c(all_sig_WR_WL_START,all_sig_WR_END_WL,all_sig_Wl_WL_START,all_sig_Wl_END_WL)
all_sig 
all_sig <- unique(all_sig)
all_sig 
all_sig_ <- rownames(BBDD_input)[which(rownames(BBDD_input) %in% all_sig)]
all_sig_ <- as.data.frame(all_sig_)
colnames(all_sig_)[1] <- "X"
all_sig_



# Load statistics of selected DE genes by experimental group. Note: these input files have been generated by the code Part_1.



# For group_2

all_sig_WR_WL_START_conFC <- read.csv2("./Output_files/DE_genes/Group_2/WL_VS_START_wr_PVALUE_0_0001.csv")
all_sig_WR_WL_START_conFC[,1] <- paste(all_sig_WR_WL_START_conFC$X,"/",all_sig_WR_WL_START_conFC$SYMBOL,sep="")
SIGTMP <- all_sig_
SIGTMP_ <-merge(SIGTMP,all_sig_WR_WL_START_conFC[,c(1,6)],all.x=T,by="X")
SIGTMP_[is.na(SIGTMP_[,2]),2] <- 999
SIGTMP_
write.csv2(SIGTMP_,"./Output_files/DE_genes/Group_2/WL_VS_START_wr_FC.csv")
all_sig_WR_END_WL_conFC <- read.csv2("./Output_files/DE_genes/Group_2/end_VS_WL_wR_PVALUE_0_001.csv")
all_sig_WR_END_WL_conFC[,1] <- paste(all_sig_WR_END_WL_conFC$X,"/",all_sig_WR_END_WL_conFC$SYMBOL,sep="")
SIGTMP <- all_sig_
SIGTMP_ <-merge(SIGTMP,all_sig_WR_END_WL_conFC[,c(1,6)],all.x=T,by="X")
SIGTMP_[is.na(SIGTMP_[,2]),2] <- 999
SIGTMP_
write.csv2(SIGTMP_,"./Output_files/DE_genes/Group_2/END_WL_wr_FC.csv")



# For group_1

all_sig_Wl_WL_START_conFC <- read.csv2("./Output_files/DE_genes/Group_1/WL_VS_START_wl_PVALUE_0_0001.csv")
all_sig_Wl_WL_START_conFC[,1] <- paste(all_sig_Wl_WL_START_conFC$X,"/",all_sig_Wl_WL_START_conFC$SYMBOL,sep="")
SIGTMP <- all_sig_
SIGTMP_ <-merge(SIGTMP,all_sig_Wl_WL_START_conFC[,c(1,6)],all.x=T,by="X")
SIGTMP_[is.na(SIGTMP_[,2]),2] <- 999
SIGTMP_
write.csv2(SIGTMP_,"./Output_files/DE_genes/Group_1/WL_START_wL_FC.csv")
all_sig_Wl_END_WL_conFC <- read.csv2("./Output_files/DE_genes/Group_1/end_VS_WL_wl_PVALUE_0_001.csv")
all_sig_Wl_END_WL_conFC[,1] <- paste(all_sig_Wl_END_WL_conFC$X,"/",all_sig_Wl_END_WL_conFC$SYMBOL,sep="")
SIGTMP <- all_sig_
SIGTMP_ <-merge(SIGTMP,all_sig_Wl_END_WL_conFC[,c(1,6)],all.x=T,by="X")
SIGTMP_[is.na(SIGTMP_[,2]),2] <- 999
SIGTMP_
write.csv2(SIGTMP_,"./Output_files/DE_genes/Group_1/END_WL_wL_FC.csv")



# Load generated Foldchange files by experimental group. Note: in these files, a FC of 999 refers to a non-significant change for the particular time interval.



#Group 1 T1-vs-Baseline

LCD_WL_START__ <- read.csv2("./Output_files/DE_genes/Group_1/WL_START_wL_FC.csv",row.names=1)



#Relax FoldChanges at 5 %

LCD_WL_START__[which(is.negative(LCD_WL_START__$logFC)),2] <- LCD_WL_START__[which(is.negative(LCD_WL_START__$logFC)),2] - (LCD_WL_START__[which(is.negative(LCD_WL_START__$logFC)),2]*0.5)
LCD_WL_START__[which(!is.negative(LCD_WL_START__$logFC)),2] <- LCD_WL_START__[which(!is.negative(LCD_WL_START__$logFC)),2] - (LCD_WL_START__[which(!is.negative(LCD_WL_START__$logFC)),2]*0.5)



#Group 1 T2-vs-t1

LCD_WS_WL__ <- read.csv2("./Output_files/DE_genes/Group_1/END_WL_wL_FC.csv",row.names=1)



#Relax FoldChanges at 5 %

LCD_WS_WL__[which(is.negative(LCD_WS_WL__$logFC)),2] <- LCD_WS_WL__[which(is.negative(LCD_WS_WL__$logFC)),2] - (LCD_WS_WL__[which(is.negative(LCD_WS_WL__$logFC)),2]*0.5)
LCD_WS_WL__[which(!is.negative(LCD_WS_WL__$logFC)),2] <- LCD_WS_WL__[which(!is.negative(LCD_WS_WL__$logFC)),2] - (LCD_WS_WL__[which(!is.negative(LCD_WS_WL__$logFC)),2]*0.5)



#Group 2 T1-vs-Baseline

VLCD_WL_START__ <- read.csv2("./Output_files/DE_genes/Group_2/WL_VS_START_wr_FC.csv",row.names=1)



#Relax FoldChanges at 5 %

VLCD_WL_START__[which(is.negative(VLCD_WL_START__$logFC)),2] <- VLCD_WL_START__[which(is.negative(VLCD_WL_START__$logFC)),2] - (VLCD_WL_START__[which(is.negative(VLCD_WL_START__$logFC)),2]*0.5)
VLCD_WL_START__[which(!is.negative(VLCD_WL_START__$logFC)),2] <- VLCD_WL_START__[which(!is.negative(VLCD_WL_START__$logFC)),2] - (VLCD_WL_START__[which(!is.negative(VLCD_WL_START__$logFC)),2]*0.5)



#Group 2 T2-vs-t1


VLCD_WS_WL__ <- read.csv2("./Output_files/DE_genes/Group_2/END_WL_wr_FC.csv",row.names=1)



#Relax FoldChanges at 5 %

VLCD_WS_WL__[which(is.negative(VLCD_WS_WL__$logFC)),2] <- VLCD_WS_WL__[which(is.negative(VLCD_WS_WL__$logFC)),2] - (VLCD_WS_WL__[which(is.negative(VLCD_WS_WL__$logFC)),2]*0.5)
VLCD_WS_WL__[which(!is.negative(VLCD_WS_WL__$logFC)),2] <- VLCD_WS_WL__[which(!is.negative(VLCD_WS_WL__$logFC)),2] - (VLCD_WS_WL__[which(!is.negative(VLCD_WS_WL__$logFC)),2]*0.5)



# Subset the normalized microarray input matrix according to selected DE probes.

BBDD_input <- BBDD_input[rownames(BBDD_input) %in% LCD_WL_START__$X,]




# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #  Discretization of gene expression values by experimental condition and generation of the sequence database
#Note: The discretization process is conducted separately in each experimental condition (data scope)
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #




##################################################################################################################
# Group 1
##################################################################################################################



# Select subjects (from the normalized expression matrix) which belong to the experimental group 1.

LDL_BD <- BBDD_input[, which(gsub("\\_.*","",colnames(BBDD_input)) == "WeightLosers")] 
dim(LDL_BD)
colnames(LDL_BD)



# Preparation and inspection of the resulting dataset

LDL_BD <- data.frame(lapply(LDL_BD, as.character), stringsAsFactors=FALSE)
LDL_BD <- data.frame(lapply(LDL_BD, as.numeric), stringsAsFactors=FALSE)
str(LDL_BD)
dim(LDL_BD) 
colnames(LDL_BD) 
sapply(LDL_BD, function(y) sum(length(which(is.na(y))))) 
LDL_BD <- data.frame(lapply(LDL_BD, as.character), stringsAsFactors=FALSE)
LDL_BD <- data.frame(lapply(LDL_BD, as.numeric), stringsAsFactors=FALSE)



#Study of NAs by individual/time record

sapply(LDL_BD, function(y) sum(length(which(is.na(y))))) 
str(LDL_BD) 



#Generation of new working DB

NUEVA_LDL_BD <- LDL_BD 



# Discretization process 

eliminar <- seq(1, length(LDL_BD), by=3)
vector <- 1:length(LDL_BD)
vector <- vector[-eliminar]
for( i in vector ) {
r <- i-1
NUEVA_LDL_BD[,i] <- log2(LDL_BD[,i]/LDL_BD[,r])
}
NUEVA_LDL_BD <- NUEVA_LDL_BD[,-eliminar] # Remove all absolute baseline data and let only interval changes.
SELECCCION_T2WL_T1START <- seq(1, length(NUEVA_LDL_BD), by=2)
SELECCCION_T3WS_T2WL <- seq(0, length(NUEVA_LDL_BD), by=2)[-1]



# Comparing each individual/time/probe interval FoldChange vs the mean FoldChange of the same proble and time interval in all individuals. Note= 499.5 refers to non-significan DE genes.

average_T2WL_T1START <- apply(NUEVA_LDL_BD[,SELECCCION_T2WL_T1START],1,mean,na.rm=TRUE)
average_T2WL_T1START
average_T2WL_T1START[which(LCD_WL_START__$logFC == 499.50)] <- 499.50
average_T2WL_T1START
average_T3WS_T2WL <- apply(NUEVA_LDL_BD[,SELECCCION_T3WS_T2WL],1,mean,na.rm=TRUE)
average_T3WS_T2WL
average_T3WS_T2WL[which(LCD_WS_WL__$logFC == 499.50)] <- 499.50
average_T3WS_T2WL




# Below, the values of change in gene expression are substituted by 1,2 o 3 (according to the next coding scheme) agreed on if they overcome or not the mean foldchanges previously estimated.
#Coding scheme
#1 = Downregulated
#2 = Up-regulated
#3 = No change

for( i in SELECCCION_T2WL_T1START) {

NUEVA_LDL_BD[which(is.negative(average_T2WL_T1START)),i][ average_T2WL_T1START[which(is.negative(average_T2WL_T1START))] >= NUEVA_LDL_BD[which(is.negative(average_T2WL_T1START)),i ]] <- -3

NUEVA_LDL_BD[which(is.negative(average_T2WL_T1START)),i][ average_T2WL_T1START[which(is.negative(average_T2WL_T1START))] < NUEVA_LDL_BD[which(is.negative(average_T2WL_T1START)),i ]] <- 0


NUEVA_LDL_BD[which(!is.negative(average_T2WL_T1START)),i][ average_T2WL_T1START[which(!is.negative(average_T2WL_T1START))] <= NUEVA_LDL_BD[which(!is.negative(average_T2WL_T1START)),i ]] <- 3

NUEVA_LDL_BD[which(!is.negative(average_T2WL_T1START)),i][ average_T2WL_T1START[which(!is.negative(average_T2WL_T1START))] > NUEVA_LDL_BD[which(!is.negative(average_T2WL_T1START)),i ]] <- 0

NUEVA_LDL_BD[,i][ NUEVA_LDL_BD[,i] == 3 ] <- 2 # Up-regulated
NUEVA_LDL_BD[,i][ NUEVA_LDL_BD[,i] == -3 ] <- 1 # Downregulated
NUEVA_LDL_BD[,i][ NUEVA_LDL_BD[,i] == 0 ] <- 3 #No change


}

for( i in SELECCCION_T3WS_T2WL) {

NUEVA_LDL_BD[which(is.negative(average_T3WS_T2WL)),i][ average_T3WS_T2WL[which(is.negative(average_T3WS_T2WL))] >= NUEVA_LDL_BD[which(is.negative(average_T3WS_T2WL)),i ]] <- -3

NUEVA_LDL_BD[which(is.negative(average_T3WS_T2WL)),i][ average_T3WS_T2WL[which(is.negative(average_T3WS_T2WL))] < NUEVA_LDL_BD[which(is.negative(average_T3WS_T2WL)),i ]] <- 0


NUEVA_LDL_BD[which(!is.negative(average_T3WS_T2WL)),i][ average_T3WS_T2WL[which(!is.negative(average_T3WS_T2WL))] <= NUEVA_LDL_BD[which(!is.negative(average_T3WS_T2WL)),i ]] <- 3

NUEVA_LDL_BD[which(!is.negative(average_T3WS_T2WL)),i][ average_T3WS_T2WL[which(!is.negative(average_T3WS_T2WL))] > NUEVA_LDL_BD[which(!is.negative(average_T3WS_T2WL)),i ]] <- 0

NUEVA_LDL_BD[,i][ NUEVA_LDL_BD[,i] == 3 ] <- 2 # Up-regulated
NUEVA_LDL_BD[,i][ NUEVA_LDL_BD[,i] == -3 ] <- 1 # Downregulated
NUEVA_LDL_BD[,i][ NUEVA_LDL_BD[,i] == 0 ] <- 3 #No change


}

rownames(NUEVA_LDL_BD) <- rownames(BBDD_input)



# The next command will show how many genes/probes change their expression by time interval. 

lapply(NUEVA_LDL_BD, table)
str(NUEVA_LDL_BD) #study of the dataframe generated so far. It contains the integer values ​​of change by time interval and individual (each column corresponds to an individual and time interval:  _2_VLCD (refers interval T2-T1) or _3_VLCD (refers interval T3-T2 per individual)



# Generate the file with correspondence between gene/probe names and order identifiers in the output sequence database

nuevos_nombres_genes <- sprintf("%03d", 1:nrow(NUEVA_LDL_BD)) 
CORRESPONDENCIA_NAMES_GENES <- data.frame(rownames(BBDD_input),nuevos_nombres_genes)
write.table(CORRESPONDENCIA_NAMES_GENES, "./Output_files/Correspondence_file/Correspondence_gene_names",sep="\t",col.names=F)  # genero un archivo que contiene en una columna los nombres de los genes iniciales "anotados" y los numeros generados por gen. Guardo este archivo en la carpeta.



# Assign new probe identifiers (numbers)

rownames(NUEVA_VLDL_BD) <- nuevos_nombres_genes 



# Preparation of final sequence database for the Group 1. This database will present as many rows as individuals within the group. The columns will represent the change in the expression of each probe/gene by time interval. The columns with -1 refers the end of a time interval while the columns with 2 refers to the end of the database.

filas <- ncol(NUEVA_LDL_BD)/2
columnas <- (nrow(NUEVA_LDL_BD)*2)+3
FINAL__LDL_BD <- data.frame(matrix(NA, nrow = filas , ncol = columnas )) 
eliminar_2 <- seq(1, ncol(NUEVA_LDL_BD), by=2)
vector_2 <- 1:ncol(NUEVA_LDL_BD)
vector_2 <- vector_2[-eliminar_2]
filas_alberto <- rep(1:filas, each=2) 
for( i in vector_2 ) {
r <- i-1
genes_individuo_t1 <- paste(as.character(NUEVA_LDL_BD[,r]), nuevos_nombres_genes, sep='')
sepa <- "-1"
genes_individuo_t2 <- paste(as.character(NUEVA_LDL_BD[,i]), nuevos_nombres_genes, sep='')
sepa <- "-1"
fin <- "-2"
transaction_1 <- c(genes_individuo_t1,sepa,genes_individuo_t2,sepa,fin)
FINAL__LDL_BD[filas_alberto[i],] <- transaction_1
}
dim(FINAL__LDL_BD)
class(FINAL__LDL_BD)
str(FINAL__LDL_BD)
FINAL__LDL_BD
NOMBRES_INDIVIDUOS <- colnames(NUEVA_LDL_BD)[seq(1, ncol(NUEVA_LDL_BD), by=2)]
NOMBRES_INDIVIDUOS <- gsub("ers_","ers",NOMBRES_INDIVIDUOS)
NOMBRES_INDIVIDUOS <- str_extract(NOMBRES_INDIVIDUOS, "[^_]+")
NOMBRES_INDIVIDUOS
FINAL__LDL_BD <- as.data.frame(lapply(FINAL__LDL_BD, as.numeric))
rownames(FINAL__LDL_BD) <- NOMBRES_INDIVIDUOS
colnames(FINAL__LDL_BD) <- as.character(CORRESPONDENCIA_NAMES_GENES[,1])
colnames(FINAL__LDL_BD)[c((length(as.character(CORRESPONDENCIA_NAMES_GENES[,1]))+2):(length(colnames(FINAL__LDL_BD))-2))] <- as.character(CORRESPONDENCIA_NAMES_GENES[,1])
head(FINAL__LDL_BD)



# Substitute by NA all genes/probes which do not change their expression.
GENES_ANODINOS_A_SUSTITUIR_POR_na <- as.numeric(paste("3",nuevos_nombres_genes,sep=""))
for( i in 1:length(GENES_ANODINOS_A_SUSTITUIR_POR_na) ) {
FINAL__LDL_BD[FINAL__LDL_BD == GENES_ANODINOS_A_SUSTITUIR_POR_na[i]] <- NA
}
head(FINAL__LDL_BD)



# Save final SEQUENCE DATA BASE OF GROUP 1.

write.csv2(x=FINAL__LDL_BD, file="./Output_files/Sequence_databases/Group_1/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_1.csv")

list <- list()
for (i in 1:nrow(FINAL__LDL_BD)) {

list_i <- as.character(FINAL__LDL_BD[i,])

list_i[list_i == "NA"] <- NA

list_i <- list_i[-which(is.na(list_i))]

list[[i]] <- list_i

}

lapply(list, write, "./Output_files/Sequence_databases/Group_1/Formated_Sequence_DB/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_1_formated", append=TRUE, ncolumns=1000)


##################################################################################################################
# Group 2
##################################################################################################################



# Select subjects (from the normalized expression matrix) which belong to the experimental group 2.

VLDL_BD <- BBDD_input[, which(gsub("\\_.*","",colnames(BBDD_input)) == "WeightRegainers")] 
dim(VLDL_BD) 
colnames(VLDL_BD) 



# Preparation and inspection of the resulting dataset

VLDL_BD <- data.frame(lapply(VLDL_BD, as.character), stringsAsFactors=FALSE)
VLDL_BD <- data.frame(lapply(VLDL_BD, as.numeric), stringsAsFactors=FALSE)



#Study of NAs by individual/time record


sapply(VLDL_BD, function(y) sum(length(which(is.na(y))))) 
str(VLDL_BD) 



#Generation of new working DB

NUEVA_VLDL_BD <- VLDL_BD



# Discretization process 

eliminar <- seq(1, length(VLDL_BD), by=3) 
vector <- 1:length(VLDL_BD)
vector <- vector[-eliminar]
for( i in vector ) {
r <- i-1
NUEVA_VLDL_BD[,i] <- log2(VLDL_BD[,i]/VLDL_BD[,r]) 
}
NUEVA_VLDL_BD <- NUEVA_VLDL_BD[,-eliminar] # Remove all absolute baseline data and let only interval changes.
SELECCCION_T2WL_T1START <- seq(1, length(NUEVA_VLDL_BD), by=2)
SELECCCION_T3WS_T2WL <- seq(0, length(NUEVA_VLDL_BD), by=2)[-1]



# Comparing each individual/time/probe interval FoldChange vs the mean FoldChange of the same proble and time interval in all individuals. Note= 499.5 refers to non-significan DE genes.

average_T2WL_T1START <- apply(NUEVA_VLDL_BD[,SELECCCION_T2WL_T1START],1,mean,na.rm=TRUE)
average_T2WL_T1START
average_T2WL_T1START[which(VLCD_WL_START__$logFC == 499.50)] <- 499.50
average_T2WL_T1START
average_T3WS_T2WL <- apply(NUEVA_VLDL_BD[,SELECCCION_T3WS_T2WL],1,mean,na.rm=TRUE)
average_T3WS_T2WL
average_T3WS_T2WL[which(VLCD_WS_WL__$logFC == 499.50)] <- 499.50
average_T3WS_T2WL




# Below, the values of change in gene expression are substituted by 1,2 o 3 (according to the next coding scheme) agreed on if they overcome or not the mean foldchanges previously estimated.
#Coding scheme
#1 = Downregulated
#2 = Up-regulated
#3 = No change

for( i in SELECCCION_T2WL_T1START) {

NUEVA_VLDL_BD[which(is.negative(average_T2WL_T1START)),i][ average_T2WL_T1START[which(is.negative(average_T2WL_T1START))] >= NUEVA_VLDL_BD[which(is.negative(average_T2WL_T1START)),i ]] <- -3

NUEVA_VLDL_BD[which(is.negative(average_T2WL_T1START)),i][ average_T2WL_T1START[which(is.negative(average_T2WL_T1START))] < NUEVA_VLDL_BD[which(is.negative(average_T2WL_T1START)),i ]] <- 0


NUEVA_VLDL_BD[which(!is.negative(average_T2WL_T1START)),i][ average_T2WL_T1START[which(!is.negative(average_T2WL_T1START))] <= NUEVA_VLDL_BD[which(!is.negative(average_T2WL_T1START)),i ]] <- 3

NUEVA_VLDL_BD[which(!is.negative(average_T2WL_T1START)),i][ average_T2WL_T1START[which(!is.negative(average_T2WL_T1START))] > NUEVA_VLDL_BD[which(!is.negative(average_T2WL_T1START)),i ]] <- 0

NUEVA_VLDL_BD[,i][ NUEVA_VLDL_BD[,i] == 3 ] <- 2 # Up-regulated
NUEVA_VLDL_BD[,i][ NUEVA_VLDL_BD[,i] == -3 ] <- 1 # Downregulated
NUEVA_VLDL_BD[,i][ NUEVA_VLDL_BD[,i] == 0 ] <- 3 #No change


}

for( i in SELECCCION_T3WS_T2WL) {

NUEVA_VLDL_BD[which(is.negative(average_T3WS_T2WL)),i][ average_T3WS_T2WL[which(is.negative(average_T3WS_T2WL))] >= NUEVA_VLDL_BD[which(is.negative(average_T3WS_T2WL)),i ]] <- -3

NUEVA_VLDL_BD[which(is.negative(average_T3WS_T2WL)),i][ average_T3WS_T2WL[which(is.negative(average_T3WS_T2WL))] < NUEVA_VLDL_BD[which(is.negative(average_T3WS_T2WL)),i ]] <- 0


NUEVA_VLDL_BD[which(!is.negative(average_T3WS_T2WL)),i][ average_T3WS_T2WL[which(!is.negative(average_T3WS_T2WL))] <= NUEVA_VLDL_BD[which(!is.negative(average_T3WS_T2WL)),i ]] <- 3

NUEVA_VLDL_BD[which(!is.negative(average_T3WS_T2WL)),i][ average_T3WS_T2WL[which(!is.negative(average_T3WS_T2WL))] > NUEVA_VLDL_BD[which(!is.negative(average_T3WS_T2WL)),i ]] <- 0

NUEVA_VLDL_BD[,i][ NUEVA_VLDL_BD[,i] == 3 ] <- 2 # Up-regulated
NUEVA_VLDL_BD[,i][ NUEVA_VLDL_BD[,i] == -3 ] <- 1 # Downregulated
NUEVA_VLDL_BD[,i][ NUEVA_VLDL_BD[,i] == 0 ] <- 3 #No change

}
rownames(NUEVA_VLDL_BD) <- rownames(BBDD_input)



# The next command will show how many genes/probes change their expression by time interval. 

lapply(NUEVA_VLDL_BD, table)

str(NUEVA_VLDL_BD) #study of the dataframe generated so far. It contains the integer values ​​of change by time interval and individual (each column corresponds to an individual and time interval:  _2_VLCD (refers interval T2-T1) or _3_VLCD (refers interval T3-T2 per individual)



# Generate the file with correspondence between gene/probe names and order identifiers in the output sequence database

nuevos_nombres_genes <- sprintf("%03d", 1:nrow(NUEVA_VLDL_BD)) 
CORRESPONDENCIA_NAMES_GENES <- data.frame(rownames(BBDD_input),nuevos_nombres_genes)
write.table(CORRESPONDENCIA_NAMES_GENES, "./Output_files/Correspondence_file/Correspondence_gene_names",sep="\t",col.names=F)  # genero un archivo que contiene en una columna los nombres de los genes iniciales "anotados" y los numeros generados por gen. Guardo este archivo en la carpeta.



# Assign new probe identifiers (numbers)

rownames(NUEVA_VLDL_BD) <- nuevos_nombres_genes  



# Preparation of final sequence database for the Group 1. This database will present as many rows as individuals within the group. The columns will represent the change in the expression of each probe/gene by time interval. The columns with -1 refers the end of a time interval while the columns with 2 refers to the end of the database.

filas <- ncol(NUEVA_VLDL_BD)/2
columnas <- (nrow(NUEVA_VLDL_BD)*2)+3
FINAL__VLDL_BD <- data.frame(matrix(NA, nrow = filas , ncol = columnas )) 
eliminar_2 <- seq(1, ncol(NUEVA_VLDL_BD), by=2)
vector_2 <- 1:ncol(NUEVA_VLDL_BD)
vector_2 <- vector_2[-eliminar_2]
filas_alberto <- rep(1:filas, each=2) 
for( i in vector_2 ) {
r <- i-1
genes_individuo_t1 <- paste(as.character(NUEVA_VLDL_BD[,r]), nuevos_nombres_genes, sep='')
sepa <- "-1"
genes_individuo_t2 <- paste(as.character(NUEVA_VLDL_BD[,i]), nuevos_nombres_genes, sep='')
sepa <- "-1"
fin <- "-2"
transaction_1 <- c(genes_individuo_t1,sepa,genes_individuo_t2,sepa,fin)
FINAL__VLDL_BD[filas_alberto[i],] <- transaction_1
}
dim(FINAL__VLDL_BD)
class(FINAL__VLDL_BD)
str(FINAL__VLDL_BD)
FINAL__VLDL_BD[1:5,3457:3461]
NOMBRES_INDIVIDUOS <- colnames(NUEVA_VLDL_BD)[seq(1, ncol(NUEVA_VLDL_BD), by=2)]
NOMBRES_INDIVIDUOS <- gsub("ers_","ers",NOMBRES_INDIVIDUOS)
NOMBRES_INDIVIDUOS <- str_extract(NOMBRES_INDIVIDUOS, "[^_]+")
NOMBRES_INDIVIDUOS
FINAL__VLDL_BD <- as.data.frame(lapply(FINAL__VLDL_BD, as.numeric))
rownames(FINAL__VLDL_BD) <- NOMBRES_INDIVIDUOS
colnames(FINAL__VLDL_BD) <- as.character(CORRESPONDENCIA_NAMES_GENES[,1])
colnames(FINAL__VLDL_BD)[c((length(as.character(CORRESPONDENCIA_NAMES_GENES[,1]))+2):(length(colnames(FINAL__VLDL_BD))-2))] <- as.character(CORRESPONDENCIA_NAMES_GENES[,1])
colnames(FINAL__VLDL_BD)
head(FINAL__VLDL_BD)



# Substitute by NA all genes/probes which do not change their expression.

GENES_ANODINOS_A_SUSTITUIR_POR_na <- as.numeric(paste("3",nuevos_nombres_genes,sep=""))
for( i in 1:length(GENES_ANODINOS_A_SUSTITUIR_POR_na) ) {
FINAL__VLDL_BD[FINAL__VLDL_BD == GENES_ANODINOS_A_SUSTITUIR_POR_na[i]] <- NA
}
head(FINAL__VLDL_BD)



# Save final SEQUENCE DATA BASE OF GROUP 1.

write.csv2(x=FINAL__VLDL_BD, file="./Output_files/Sequence_databases/Group_2/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_2.csv")

list <- list()
for (i in 1:nrow(FINAL__VLDL_BD)) {

list_i <- as.character(FINAL__VLDL_BD[i,])

list_i[list_i == "NA"] <- NA

list_i <- list_i[-which(is.na(list_i))]

list[[i]] <- list_i

}

lapply(list, write, "./Output_files/Sequence_databases/Group_2/Formated_Sequence_DB/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_2_formated", append=TRUE, ncolumns=1000)

#................................. End of code
