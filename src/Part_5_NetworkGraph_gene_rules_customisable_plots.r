#GNU GENERAL PUBLIC LICENSE - Version 3
#
#Part_7_NetworkGraph_gene_rules_customisable_plots.r Copyright (C) 2019
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

library(arules)
library(arulesViz)
library(visNetwork) 
library(Rgraphviz)
library(stringr)
library(gtools)
library(tidyr)
library(Hmisc)
library(dplyr)



##############################

# DATA LOAD AND PREPROCESSING

##############################

# Set working directory

setwd("./GeneSeqRules")


# Note: change the next file path for repeating the process with each group output rule.

Group_X_rules <- read.PMML("./Output_files/Output_rules/Output_rules_more_measures/pmml/Output_rules_raw_group_1_measures_biologicmeasures.pmml")
Group_X_rules



# Convert to dataframe

Group_X_rules <- DATAFRAME(Group_X_rules)
str(Group_X_rules)
rules <- Group_X_rules
dim(rules)



# Select Strong ARs

rules <- subset(rules, CF >= 0 & LIFT >= 1 & CONF >= 0.7 )

dim(rules)

rules[,1] <- sapply(rules[,1],function(x) gsub('^.|.$','',as.character(x)))
rules[,2] <- sapply(rules[,2],function(x) gsub('^.|.$',"",as.character(x))) 

NUM_lhs <- sapply(rules[,1],str_count, pattern = ",")+1 
NUM_rhs <- sapply(rules[,2],str_count, pattern = ",")+1

to_funct <- function(x){ 1:x }

NUM_lhs_ <- sapply(NUM_lhs,to_funct)
NUM_rhs_ <- sapply(NUM_rhs,to_funct)

a <- 0

for (i in 1:nrow(rules)) {

	b <- nrow(crossing(as.data.frame(NUM_lhs_[i]),as.data.frame(NUM_rhs_[i])))
	a <- a + b

				}



NUM_lhs <- sapply(rules[,1],str_count, pattern = ",")+1 
NUM_rhs <- sapply(rules[,2],str_count, pattern = ",")+1


RES_M <- as.data.frame(matrix(NA,a,ncol(rules)))
names(RES_M) <- names(rules)
str(RES_M)


a <- 1
for ( i in 1:nrow(rules) ) {


	if ( !grepl(",", rules[i,1]) & !grepl(",", rules[i,2]) ) { #No and No


		RES_M[a,1] <- rules[i,1]
		RES_M[a,2] <- rules[i,2]
		RES_M[a,3:ncol(rules)] <- rules[i,3:ncol(rules)]
		a <- a+1


		} else if ( grepl(",", rules[i,1]) & !grepl(",", rules[i,2]) ) { #Yes and No
	

			lhss <- str_split_fixed(rules[i,1], ",", n=NUM_lhs[i])

			n <- nrow(expand.grid(lhss, rules[i,2]))
			
			RES_M[a:(a-1+n),1:2] <- as.matrix(expand.grid(lhss, rules[i,2]))
			
			RES_M[a:(a-1+n),3:ncol(rules)] <- rules[i,3:ncol(rules)]

			a <- a+n


			}  else if ( grepl(",", rules[i,1]) & grepl(",", rules[i,2]) ) { #Yes Yes


				lhss <- str_split_fixed(rules[i,1], ",", n=NUM_lhs[i])
				rhss <- str_split_fixed(rules[i,2], ",", n=NUM_rhs[i])

				n <- nrow(expand.grid(lhss, rhss))

				RES_M[a:(a-1+n),1:2] <- as.matrix(expand.grid(lhss, rhss))

				RES_M[a:(a-1+n),3:ncol(rules)] <- rules[i,3:ncol(rules)]

				a <- a+n


				}  else if ( !grepl(",", rules[i,1]) & grepl(",", rules[i,2]) ) { #No Yes


					rhss <- str_split_fixed(rules[i,2], ",", n=NUM_rhs[i])

					n <- nrow(expand.grid(rules[i,1], rhss))
			
					RES_M[a:(a-1+n),1:2] <- as.matrix(expand.grid(rules[i,1], rhss))
			
					RES_M[a:(a-1+n),3:ncol(rules)] <- rules[i,3:ncol(rules)]

					a <- a+n


					} else { }


	}









RES_M$LHS <- gsub("\\|[^\\]]*\\=", "=", RES_M$LHS, perl=TRUE)
RES_M$RHS <- gsub("\\|[^\\]]*\\=", "=", RES_M$RHS, perl=TRUE)



# CREATION OF DF1 (EDGES)


RES_M


dflhs <- RES_M[,c(1,3:ncol(RES_M))]
dfrhs <- RES_M[,c(2:ncol(RES_M))]

names(dflhs)[1] <- "UNIQUE_PROBES"
names(dfrhs)[1] <- "UNIQUE_PROBES"

dfmerge <- rbind(dflhs,dfrhs)

dfmerge[,1] <- as.factor(dfmerge[,1])



dfmerge_means <- dfmerge 
dfmerge_means <- group_by(dfmerge_means,UNIQUE_PROBES)
dfmerge_means <- summarise_all(dfmerge_means,"mean")
names(dfmerge_means)[2:ncol(dfmerge_means)] <- paste(names(dfmerge_means)[2:ncol(dfmerge_means)],"_mean",sep="")
print(dfmerge_means, n=nrow(dfmerge_means))
dfmerge_means <- as.data.frame(dfmerge_means)






PROBES <- gsub("\\/.*","",dfmerge_means$UNIQUE_PROBES)
PROBES
UNIQUE_GENES <- str_match(dfmerge_means$UNIQUE_PROBES, "/(.*?)=")[,2]
UNIQUE_GENES
REGULATION_TYPE <- sub(".*=", "", dfmerge_means$UNIQUE_PROBES)
UNIQUE_GENES <- gsub("\\|.*","",UNIQUE_GENES)
UNIQUE_GENES


DF1 <- data.frame(dfmerge_means[,1],PROBES,UNIQUE_GENES,REGULATION_TYPE,dfmerge_means[,c(2:ncol(dfmerge_means))])
names(DF1)[1] <- "UNIQUE_PROBES"


ngrp <- 2 # Customisable: This argument will change according to the selected grouping variable. Here it's has been set as 2, since we have chosen the type of gene expression change (upregulation/downregulation) as grouping variable.
matrix_of_groups <- as.data.frame(matrix(NA,ngrp+1,dim(DF1)[2]))
names(matrix_of_groups) <- names(DF1)

DF1_ <- rbind(matrix_of_groups, DF1)
DF1_[1:(ngrp+1),1] <- c(0,1:ngrp)
DF1_[2:(ngrp+1),4] <- c(0,1) # Customisable: The number four refers to the selected column. 
DF1_$UNIQUE_PROBES <- as.factor(DF1_$UNIQUE_PROBES)
DF1_$REGULATION_TYPE <- as.factor(DF1_$REGULATION_TYPE)

str(DF1_)

GROUPs <- seq(0, 5, by = 0.10)  #INTERVALS FOR GROUPS BASED ON BIOLOGICAL MEASURES.

MEASURES_CAT <- as.data.frame(lapply(DF1_[,c("BP_mean","MF_mean","CC_mean","SP_mean","TF_mean")],function(x) cut2(x,GROUPs)))
MEASURES_CAT <- lapply(MEASURES_CAT,droplevels)
 
names(MEASURES_CAT) <- paste(names(MEASURES_CAT),"_cat",sep="")

DF1_def <- data.frame(DF1_,MEASURES_CAT)

DF1_def <- DF1_def[c((1+1):nrow(DF1_def)),c("REGULATION_TYPE","UNIQUE_PROBES")] # Customisable: In the case of changing the grouping variable, the colname "REGULATION_TYPE" will change according to the new selected grouping variable.

str(DF1_def)

DF1_def 

edges <- DF1_def

edges <- edges[order(edges$REGULATION_TYPE,edges$UNIQUE_PROBES),] # Customisable:  In the case of changing the grouping variable, the colname "REGULATION_TYPE" will change according to the new selected grouping variable.

edges #DF1 (VERTICES) ***************** FIRST OUTPUT








# CREATION OF DF2 (VERTICES)

RES_M


dflhs <- RES_M[,c(1,3:ncol(RES_M))]
dfrhs <- RES_M[,c(2:ncol(RES_M))]

names(dflhs)[1] <- "UNIQUE_PROBES"
names(dfrhs)[1] <- "UNIQUE_PROBES"

dfmerge <- rbind(dflhs,dfrhs)

dfmerge[,1] <- as.factor(dfmerge[,1])



dfmerge_means <- dfmerge 
dfmerge_means <- group_by(dfmerge_means,UNIQUE_PROBES)
dfmerge_means <- summarise_all(dfmerge_means,"mean")
names(dfmerge_means)[2:ncol(dfmerge_means)] <- paste(names(dfmerge_means)[2:ncol(dfmerge_means)],"_mean",sep="")
print(dfmerge_means, n=nrow(dfmerge_means))
dfmerge_means <- as.data.frame(dfmerge_means)





PROBES <- gsub("\\/.*","",dfmerge_means$UNIQUE_PROBES)
PROBES
UNIQUE_GENES <- str_match(dfmerge_means$UNIQUE_PROBES, "/(.*?)=")[,2]
UNIQUE_GENES
REGULATION_TYPE <- sub(".*=", "", dfmerge_means$UNIQUE_PROBES)
UNIQUE_GENES <- gsub("\\|.*","",UNIQUE_GENES)
UNIQUE_GENES


DF2 <- data.frame(dfmerge_means[,1],PROBES,UNIQUE_GENES,REGULATION_TYPE,dfmerge_means[,c(2:ncol(dfmerge_means))])
names(DF2)[1] <- "UNIQUE_PROBES"

DF2 <- DF2[order(DF2$REGULATION_TYPE,DF2$UNIQUE_PROBES),] # Customisable: In the case of changing the grouping variable, the colname "REGULATION_TYPE" will change according to the new selected grouping variable.



ngrp <- 2 # Customisable: This argument will change according to the selected grouping variable. Here it's has been set as 2, since we have chosen the type of gene expression change (upregulation/downregulation) as grouping variable.
matrix_of_groups <- as.data.frame(matrix(NA,ngrp+1,dim(DF2)[2]))
names(matrix_of_groups) <- names(DF2)

DF2_ <- rbind(matrix_of_groups, DF2)
DF2_[1:(ngrp+1),1] <- c(0,1:ngrp)
DF2_[2:(ngrp+1),4] <- c(0,1) # Customisable: The number four refers to the selected column. 
DF2_$UNIQUE_PROBES <- as.factor(DF2_$UNIQUE_PROBES)
DF2_$REGULATION_TYPE <- as.factor(DF2_$REGULATION_TYPE)

str(DF2_)

GROUPs <- seq(0, 5, by = 0.10)  #INTERVALS FOR GROUPS BASED ON BIOLOGICAL MEASURES.

MEASURES_CAT <- as.data.frame(lapply(DF2_[,c("BP_mean","MF_mean","CC_mean","SP_mean","TF_mean")],function(x) cut2(x,GROUPs)))
MEASURES_CAT <- lapply(MEASURES_CAT,droplevels)
 
names(MEASURES_CAT) <- paste(names(MEASURES_CAT),"_cat",sep="")

DF2_def <- data.frame(DF2_,MEASURES_CAT)

DF2_def #DF2 (VERTICES) ***************** SECOND OUTPUT

str(DF2_def)

vertices <- DF2_def


#CREATION OF DF3 (CONNECTIONS)

RES_M

GROUPs <- seq(0, 5, by = 0.25)  

MEASURES_CAT <- as.data.frame(lapply(RES_M[,c("BP","MF","CC","SP","TF")],function(x) cut2(x,GROUPs)))
MEASURES_CAT <- lapply(MEASURES_CAT,droplevels)
 
names(MEASURES_CAT) <- paste(names(MEASURES_CAT),"_cat",sep="")

DF1_def <- data.frame(RES_M,MEASURES_CAT)

DF1_def #DF2 (VERTICES) ***************** THIRD OUTPUT

str(DF1_def)

connect <- DF1_def










#LAST CHANGES

names(vertices)[which(names(vertices) == "REGULATION_TYPE")] <- "group"

edges[,2] <- as.factor(edges[,2])
names(edges) <- c("from","to")


vertices[,1] <- as.factor(vertices[,1])
names(vertices)[1] <- c("name")

vertices$id <- c(rep(NA,ngrp),1:(nrow(vertices)-ngrp))

connect[,1] <- as.factor(connect[,1])
connect[,2] <- as.factor(connect[,2])
names(connect)[1:2] <- c("from","to")



str(edges)
str(vertices)
str(connect)








connect <- connect[order(connect$BP, decreasing=T),]  # Order the dataframe according to a variable of interest in such way that rules with best metric values are plotted in the external layers of the graph. 



#VARIABLE for edge WIDTH

width_var <- factor(connect$MF_cat,levels(connect$MF_cat)[c(5,4,3,2,1)]) # Customisable: MF_cat may be changed by the variable selected for controlling edge width.
width_var <- as.numeric(width_var)

#SCALE TO [0.3,1]
width_var <-  (1 - 0.3) * ((width_var - min(width_var))/(max(width_var) - min(width_var))) + 0.3





##############################

#PLOT Resulting dataframes- V1 Orange palette

##############################

# The colour of the edges here is controlled by the variable BP_cat
# Vertex size is controlled by the variable LIFT.
# Edge width is controlled by the variable MF_cat

library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)



nleaves <- length(1:(nrow(vertices)-ngrp)) 
vertices$angle= round((90 - 360 * vertices$id / nleaves),1)
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
vertices$hjust<-ifelse( vertices$angle < -90, 0, 1) ####NOTA (A DIFERENCIA DEL LCD), AQUI HE CAMBIADO ESTO A LA FUERZA
# flip angle BY to make them readable
vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)


l <- length(levels(connect$BP_cat))

coul = brewer.pal(5, "Oranges") #You can chose a different colour palette (for example: "Spectral")  or "BuPu")
coul = colorRampPalette(coul)(l)
connect$colors = rep(NA,nrow(connect))

#PARA VER LOS COLORES DE UNA PALETA EN R.
pie(rep(1, length(coul)), col = coul , main="") 

for (i in 1:l) {

level <-  rev(levels(connect$BP_cat))[i] # Customisable: change BP if you want another variable control edge colours
colour <- coul[i]

connect$colors[ which(connect$BP_cat == level) ] <- colour # Customisable: change BP if you want another variable control edge colours

}






# Create a graph object
mygraph <- graph_from_data_frame( edges, vertices=vertices )

# The connection object must refer to the ids of the leaves:
from = match( connect$from, vertices$name)
to = match( connect$to, vertices$name)




options(expression = 100000000000000000000000)

edge_bundle_plot <- ggraph(mygraph, layout="dendrogram", circular=TRUE) 

for (i in 1:nrow(connect)) { # Note: if there are too many rules, there could be memory overload problems. For solving the issue, please, plot only a subset of the best rules (e.g (i in 150:nrow(connect)) )
  fromm <- from[i]
  too <- to[i]
  colorin <- connect$colors[i]
  alpha <- width_var[i]
  edge_type <- ifelse( connect$TF > 0, 4, 1 )
  edge_type <- edge_type[i]
  edge_bundle_plot <- edge_bundle_plot +        
                      geom_conn_bundle(data=get_con(from = fromm, to = too), 
                                       colour=colorin,arrow=arrow(angle = 15, length = unit(0.07, "inches"), ends = "last", type = "closed"),
                                       edge_width=alpha,edge_linetype=edge_type,
                                       tension = 1) 
} 



# If you wan to change the selected colours, choose new identifiers in "https://htmlcolorcodes.com/" and change the code scale_colour_manual



edge_bundle_plot +   geom_node_text(aes(x = x*1.09, y=y*1.09, filter = leaf, label=name, angle = angle, hjust=hjust, colour=group,  fontface = "bold"), size=2.6, alpha=1) +
  
  geom_node_point(aes(filter = leaf, x = x*1.04, y=y*1.04, colour=group, size=CF_mean, alpha=0.2)) +
  scale_colour_manual(values= c("#3498db","#154360")) +
  scale_size_continuous( range = c(2,11) ) +
  
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))



ggsave("./Output_Plots/NetworkGraph_gene_gene_relations/colconexion_BPcat_widthconex_MF_sizePoint_CFmean_colpoint_upregulation_GROUP_X_oranges.pdf")

















