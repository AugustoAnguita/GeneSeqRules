#GNU GENERAL PUBLIC LICENSE - Version 3
#
#Part_6_generating_annotation_file_with_go_and_kegg_terms_by_probe.r Copyright (C) 2019
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

require("hgu133plus2.db")
require("pd.hg.u133.plus.2") # Note: These packages should be changed according to the affymetrix platform under study.
library("plyr")
library("readr")

# Set working directory

setwd("./GeneSeqRules")



# Load input fluorescence raw .cel files.

BBDD_input <- read.table("./Output_files/Correspondence_file/Correspondence_gene_names", header=FALSE)
BBDD_input[,2] <- as.character(BBDD_input[,2])



# Prepare probe names

lista_affy_human_probeids <- gsub("\\/.*","",BBDD_input[,2])
length(lista_affy_human_probeids)
rm(BBDD_input)



# Check the number of probes availables in the package.

summary(lista_affy_human_probeids %in% keys(hgu133plus2.db))  # Note: These packages should be changed according to the affymetrix platform under study.
columns(hgu133plus2.db)  # Note: These packages should be changed according to the affymetrix platform under study.
tlength <- length(lista_affy_human_probeids)
interval <- tlength/20
secuencias <- seq(from = interval, to = tlength, by = interval)
secuencias
j = 1
for (i in secuencias) {
ae.annots <- AnnotationDbi::select(x= hgu133plus2.db,
  keys    = lista_affy_human_probeids[c(j:floor(i))],  columns = c("SYMBOL","GO","GOALL","PATH"),  keytype = "PROBEID"  ) # Note: These packages should be changed according to the affymetrix platform under study.
dim(ae.annots)
head(ae.annots)
print(j)
print(floor(i))
j = ceiling(i)
str(ae.annots)
ae.annots <- ae.annots[,-c(4,7)] 
ae.annots <- unique(ae.annots)
str(ae.annots)
write.csv2(ae.annots, file = paste("./Data/ae_annot/building/ae.annots",i,".csv"))
}
mydir = "./Data/ae_annot/building"
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
myfiles
dat_csv = ldply(myfiles, read_csv2)
dat_csv = dat_csv[,-1]
dat_csv = dat_csv[,-8]
dat_csv = unique(dat_csv)
write.csv2(dat_csv, file = paste("./Data/ae_annot/final/ae.annots.csv"))

