################################################################################################
## <---------------------------- Script Width for Screen -----------------------------------> ##

########################################## Contents by line ####################################
05  Contents
23  Citation and Authorship
29  Package Intall 
35  Remove previous Variables
41  Set Directory
35  Clear Previous Variables
49  Read Files
72  Data Normalization/Normalisation
76  Data Testing
88  Box Plotting of inputs
127 Dendrograms
167 Heatmap Correlations
189 PCA
207 Pairwise Variance

################################################################################################
################################# Citation and Authorship ######################################
################################################################################################
# Author: Jared Streich 12th of Feb this foul year of our lord, 2015
# Citation: The interwebz, Short Bursts of Continuity in Neural Conductivity

################################################################################################
################################# Call Packages as Needed ######################################
################################################################################################
library(cluster) 
library(dynamicTreeCut)
library(ape)

########### Clear all remaining variables in R, 
########### CAREFULL... you don't ruin another script running!!! 
rm(list = ls())


################################################################################################
############################## Change Working Directory ########################################
################################################################################################
setwd("/Users/u5212257/Desktop/")



################################################################################################
#################################### Read Data Files ###########################################
################################################################################################

startF <- 0
startF
startF <- read.delim("365_BrachyClimates.txt")
summary(startF)
dim(startF)


################################################################################################
############################# Z-Scores of Column Variables #####################################
################################################################################################

######################## Define Data set by Columns ###########################
row.names(startF)=startF[,1]

###### First and last Data Columns are.....
f=4
l=22

###### Define Column Data 
startF1 = startF[,f:l]


############################# Normalize/Normalise Data ########################
startF1 <- scale(startF1, center = TRUE, scale = TRUE)


################################ Test Data File ###############################
image(startF1)
summary(startF1)
boxplot(startF1)


############################## Write Scaled Data to file ######################
DataSheet <- file("Scaled_Summary_Data.txt")
writeLines(summary(startF1), DataSheet)
close(DataSheet)


########################## Test boxplot of variables #######################
GhE <- scale(startF1)
summary(GhE)
boxplot(GhE)
dim(GhE)

################# PDF of boxplot of Climate/Environmental Variables #############
pdf(file = "Boxplots.pdf", width = 20, height = 12)
par(mar=c(20,10,10,5) + 0.1)
boxplot(GhE, cex.axis = 0.8,las = 2)
title(main = "BoxPlot of Environmental Variables", ylab = "Z-Scores", cex.axis = "0.3")
dev.off()


####################### Test All Pairs of Variables #########################
pairs(startF1, cex = 0.1, fig=TRUE)


############################## Make Pairs Fig PDF ###########################
pdf(file = "Variable_Pairs.pdf", width = 50, height = 50)
pairs(startF1, cex = 1.7, fig=TRUE)
dev.off()

############################ Test Specific Pairs ############################
pairs(startF1[,4:8], cex = 0.1, fig=TRUE)

########################## Make Specific Pairs PDF ##########################
pdf(file = "SPecific_Variable_Pairs.pdf", width = 50, height = 50)
pairs(startF1[,4:8], cex = 1.7, fig=TRUE)
dev.off()


############################# Check Dimensions of GhE #######################
dim(GhE)

################################################################################################
############################# Various Dendrogram Outputs #######################################
################################################################################################

###################### Convert matrix to clusterable data #####################
Ghc <- hclust(dist(GhE))


####################### Test Groups in Unrooted tree ##########################
plot(as.phylo(Ghc), type = "unrooted")


################################ Define Groups ################################
groups <- 10


############################## PdF Unrooted tree ##############################
pdf(file="Unrooted_Tree.pdf")
plot(as.phylo(Ghc), type = "unrooted")
dev.off()


############### Check for groups in Typical plotting Dendrogram ###############
plot((Ghc), cex=0.5)
rect.hclust(Ghc,groups)


##################### Create PDF of Typical plotting Dendrogram ###############
pdf(file = "Resolution_Distance_Dendro.pdf", width = 20, height = 10)
plot((Ghc), cex=0.5)
rect.hclust(Ghc,groups)
dev.off()


################################ PDF of Grouped ###############################
pdf(file = "NoHang_Dendrogram.pdf", width = 20, height = 10)
plot((Ghc), cex=0.8, hang = -1)
rect.hclust(Ghc, groups)
dev.off()

################################################################################################
############################# Heatmaps with Dendrograms ########################################
################################################################################################

################################# Test heatmap #################################
x  <- as.matrix(GhE)
hv <- heatmap(x, col = topo.colors(512), scale = "column", cexCol=0.5, cexRow = 0.1,
              margins = c(10,20),
              xlab = "Column_Input", ylab =  "Row_Input", cex = 1.5,
              main = "Correlation of Variables")
utils::str(hv) # the two re-ordering index vectors

################################ Create Heatmap PDF ############################
pdf(file="Cluster_HeatMap.pdf", width = 15, height = 25)
par(mar=c(30,10,20,5) + 0.1)
x  <- as.matrix(GhE)
hv <- heatmap(x, col = topo.colors(512), scale = "column", cexCol=1.3, cexRow = 1.8,
              margins = c(10,20),
              xlab = ".", ylab =  "Populations", cex = 1.5,
              main = "Correlation of Variables")
utils::str(hv) # the two re-ordering index vectors
dev.off()


################################################################################################
########################### Principle Component Analysis #######################################
################################################################################################
PCA.res<-prcomp(GhE, scale.=F, axes = F)
biplot(PCA.res)
summary(PCA.res)


##################################### Make PCA PDF ############################
pdf(file="PCA.pdf", width = 20, height = 8)
PCA.res<-prcomp(GhE, scale.=F, axes = F)
biplot(PCA.res)
summary(PCA.res)
dev.off()


################################################################################################
############################ Pairwises Variance Analysis #######################################
################################################################################################

#  Transpose, if so remove sharp sign on next line. Probably don't need to transpose the matrix
ET <- t(GhE)

ET <- GhE

###################################### Test of Normal Variance ##############################
V <- var(ET, y = ET, na.rm = FALSE, use = "everything")
x  <- as.matrix(V)
hv <- heatmap(x, col = heat.colors(512), scale = "column",
              margins = c(10,20),
              xlab = ".", ylab =  "Pairwise Comparison",
              main = "Pairwise Comparison")
utils::str(hv) # the two re-ordering index vectors


####################################### PDF of Normal Variance ##############################
pdf(file="Variance_Kinship_Variance.pdf", width = 20, height = 20)
x  <- as.matrix(V)
hv <- heatmap(x, col = heat.colors(512), scale = "column",
              margins = c(10,20),
              xlab = ".", ylab =  "Pairwise Comparison",
              main = "Pairwise Comparison")
utils::str(hv) # the two re-ordering index vectors
dev.off()


################################################################################################
###################################### Alternate Variance Equations ############################
################################################################################################


################################ Test Pearson Variance ################################
# Pearson Variance: Least Strict, matrix signal is more weak, low matrix image contrast
P <- cor(ET, y = ET, use = "everything", method = c("pearson"))
x  <- as.matrix(P)
hv <- heatmap(x, col = heat.colors(512), scale = "column",
              margins = c(10,20),
              xlab = ".", ylab =  "Pairwise Comparison",
              main = "Pairwise Comparison")
utils::str(hv) # the two re-ordering index vectors


##################################### Pearson Variance PDF ############################
pdf(file="Pearson_Kinship_Variance.pdf", width = 20, height = 20)
x  <- as.matrix(P)
hv <- heatmap(x, col = heat.colors(512), scale = "column",
              margins = c(10,20),
              xlab = ".", ylab =  "Pairwise Comparison",
              main = "Pairwise Comparison")
utils::str(hv) # the two re-ordering index vectors
dev.off()


################################ Test Kendall Variance ################################
# Kendall Variance: Medium Strictness, matrix signal is more neutral, medium matrix image contrast
K <- cor(ET, y = ET, use = "everything", method = c("kendall"))
x  <- as.matrix(K)
hv <- heatmap(x, col = heat.colors(512), scale = "column",
              margins = c(10,20),
              xlab = ".", ylab =  "Pairwise Comparison",
              main = "Pairwise Comparison")
utils::str(hv) # the two re-ordering index vectors


##################################### Kendall Variance PDF ############################
pdf(file="Kendall_Kinship_Variance.pdf", width = 20, height = 20)
x  <- as.matrix(K)
hv <- heatmap(x, col = heat.colors(512), scale = "column",
              margins = c(10,20),
              xlab = ".", ylab =  "Pairwise Comparison",
              main = "Pairwise Comparison")
utils::str(hv) # the two re-ordering index vectors
dev.off()


################################ Test Spearman Variance ################################
# Spearman Variance: Most Strict, matrix signal is more dynamic, high matrix image contrast
S <- cor(ET, y = ET, use = "everything", method = c("spearman")) 
x  <- as.matrix(S)
hv <- heatmap(x, col = heat.colors(512), scale = "column",
              margins = c(10,20),
              xlab = ".", ylab =  "Pairwise Comparison",
              main = "Pairwise Comparison")
utils::str(hv) # the two re-ordering index vectors


##################################### Spearman Variance PDF ############################
pdf(file="Spearman_Kinship_Variance.pdf", width = 20, height = 20)
x  <- as.matrix(S)
hv <- heatmap(x, col = heat.colors(512), scale = "column",
              margins = c(10,20),
              xlab = ".", ylab =  "Pairwise Comparison",
              main = "Pairwise Comparison")
utils::str(hv) # the two re-ordering index vectors
dev.off()

################################################################################################
##################################### End of Script Joke #######################################
################################################################################################
joke <- smart_people_in_room/penguins_who_can_not_swim_so_well
