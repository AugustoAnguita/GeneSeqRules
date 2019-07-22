#GNU GENERAL PUBLIC LICENSE - Version 3
#
#Code_1_Affy_gene_expression_Preprocess_and_analysis.r Copyright (C) 2019
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
