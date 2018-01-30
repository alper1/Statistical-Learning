#Statistical Learning Project Part 3#

setwd("C:/Users/alper/Desktop")
data_credit <- read.csv("data_credit_without_outliers.csv")
head(data_credit)

##################################################################################################################
#Prepare data without qualitative variables and response variable
data_credit_for_pca <- data_credit
names(data_credit_for_pca)
data_credit_for_pca<-data_credit_for_pca[,-c(2,3,4,5,24)]
data_credit_response <- data_credit
data_credit_response <- data_credit_response[,24]

data_credit_response <- as.factor(data_credit_response)
colors.credit <- c("deepskyblue4","firebrick4")[data_credit_response]
table(colors.credit)

########################################################################################################
#Correlation Relations

#install.packages("corrplot")
library(corrplot)
m_cor<-cor(data_credit_for_pca[,-ncol(data_credit_for_pca)])
corrplot(m_cor,order="hclust",method="pie",type="full",tl.cex=0.6)

#dimension before reduction
dim(data_credit_for_pca)

##################################################################################################################
# Principal Component Analysis with the College data set
##################################################################################################################
##################################################################################################################

PCS.data_credit <- prcomp(data_credit_for_pca)

# Have a look at the outputs
names(PCS.data_credit)

# The PC scores is are in x.
head(PCS.data_credit$x)

# Make a plot of the first two PCs
plot(PCS.data_credit$x[,1:2],pch=20,col=colors.credit)

# We can add the names of the cancer
#***text(PCS.data_credit$x[,1:2],labels=data_credit_response,pos = 1,col="firebrick4",cex=0.5)

# Interpretation of the first PC: Weights for the first PC

plot(1:ncol(data_credit_for_pca),PCS.data_credit$rotation[,1],pch=20,col="deepskyblue4",main="Weights for the first PC",xlab = "Variables",ylab="Coefficients")
abline(h=0)
text(1:ncol(data_credit_for_pca),PCS.data_credit$rotation[,1],labels=colnames(data_credit_for_pca),pos=1,col="firebrick4",cex=0.5)

##################################################################################################################
# Interpretation of the second PC: Weights for the second PC

plot(1:ncol(data_credit_for_pca),PCS.data_credit$rotation[,2],pch=20,col="deepskyblue4",main="Weights for the second PC",xlab = "Variables",ylab="Coefficients")
abline(h=0)
text(1:ncol(data_credit_for_pca),PCS.data_credit$rotation[,2],labels=colnames(data_credit_for_pca),pos=1,col="firebrick4",cex=0.5)

# Interpretation of the second PC: Weights for the third PC
plot(1:ncol(data_credit_for_pca),PCS.data_credit$rotation[,3],pch=20,col="deepskyblue4",main="Weights for the third PC",xlab = "Variables",ylab="Coefficients")
abline(h=0)
text(1:ncol(data_credit_for_pca),PCS.data_credit$rotation[,3],labels=colnames(data_credit_for_pca),pos=1,col="firebrick4",cex=0.5)


# How many PCs?
# Screeplot with the 19 eigenvalues
screeplot(PCS.data_credit,npcs=19,main="Screeplot",col="deepskyblue4",type="lines",pch=19)

# Have a look at the proportion of explained variance and the cumulative proportion of explained variance
summary(PCS.data_credit)

# We need 2 PCs to have the 91% of the total variability of the data set.
# The mean of the eigenvalues can be obtained as follows
EVAL.data_credit <- PCS.data_credit$sdev^2
mean(EVAL.data_credit)

# The number of eigenvalues larger than the mean of them is 
sum(EVAL.data_credit>mean(EVAL.data_credit))

# In any case, we reduce the dimension of the data set from 19 to 2. That is we consider
# either around 10% of the number of variables in the data set keeping around the 91%
# of the information inside

# Make a plot of the first two PCs
plot(PCS.data_credit$x[,1:2],pch=20,col=colors.credit,main="First Two PCs")

# Plot the scores (the four PCs)
pairs(PCS.data_credit$x[,1:4],col=colors.credit,pch=20,main="The first four PCs")

#nrow(PCS.data_credit$x)

##################################################################################################################
# Sparse Principal Component Analysis with the College data set
##################################################################################################################

install.packages("PMA")
source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library("PMA")

# Define X as a matrix object
is(data_credit_for_pca)
data_credit_for_pca.mat <- as.matrix(data_credit_for_pca)
head(data_credit_for_pca.mat)

# Then, obtain the scaled variables
scaled <- scale(data_credit_for_pca.mat)
head(scaled)

SPC.s <- SPC(scaled,sumabsv=2.5,K=ncol(data_credit_for_pca),orth=TRUE)

# First, we see the proportion of variance explained by the SPCs.
# The formulas are different than those from the standard PCs because the SPC
# are not obtained from eigenvectors and eigenvalues (not entering into details)

prop.var.expl <- SPC.s$prop.var.explained
prop.var.expl

# As you can see, the variance explained by the first 5 SPCs is smaller than the corresponding 
# to the first 5 PCs obtained previously. Anyway, for comparison purposes, we take 5 SPCs
# The weights are given by
V <- SPC.s$v[,1:2]
V

plot(1:ncol(data_credit_for_pca),SPC.s$v[,1],pch=20,col="deepskyblue4",main="Weights for the first SPC",xlab = "Variables",ylab="Coefficients")
abline(h=0)
text(1:ncol(data_credit_for_pca),SPC.s$v[,1],labels=colnames(data_credit_for_pca),pos=1,col="firebrick4",cex=0.5)

plot(1:ncol(data_credit_for_pca),SPC.s$v[,2],pch=20,col="deepskyblue4",main="Weights for the second SPC",xlab = "Variables",ylab="Coefficients")
abline(h=0)
text(1:ncol(data_credit_for_pca),SPC.s$v[,2],labels=colnames(data_credit_for_pca),pos=1,col="firebrick4",cex=0.5)


# As it can be seen, most of the weights are equal to 0.
# Therefore, PCs are associated only to some of the variables favouring the interpretation.
# Check the important variables
# Obtain the scores for these 5 SPCs
Z.S <- scaled %*% V
head(Z.S)
colnames(Z.S) <- c("SPC1","SPC2","SPC3","SPC4","SPC5")

# Note that, altough we do impose orthogonality in the optimization problem, we do not obtain 
# orthogonal sparse PCs as in the case of the standard PCs
cor(Z.S)

# Plot the scores
pairs(Z.S,col=colors.credit,main="The first 5 Sparse PCs")

# Compare the SPC scores with the sparse PC scores
par(mfrow=c(2,2))
plot(PCS.data_credit$x[,1],Z.S[,1],main="Comparison of first scores",xlab="Fist PC",ylab="First SPC")
plot(PCS.data_credit$x[,2],Z.S[,2],main="Comparison of second scores",xlab="Second PC",ylab="Second SPC")
plot(PCS.data_credit$x[,3],Z.S[,3],main="Comparison of third scores",xlab="Third PC",ylab="Third SPC")
plot(PCS.data_credit$x[,4],Z.S[,4],main="Comparison of fourth scores",xlab="Fourth PC",ylab="Fourth SPC")


##################################################################################################################
# Independent Components Analysis with the College data set
##################################################################################################################

install.packages("fastICA")
require("fastICA")

# Obtain the ICs with the College data set. We use 5 ICs as in PCA and standarize the data
FICA.data_credit <- fastICA(data_credit_for_pca,n.comp=2,row.norm=TRUE)

# The ICA weights can be found here
FICA.data_credit$K


plot(1:ncol(data_credit_for_pca),FICA.data_credit$K[,1],pch=20,col="deepskyblue4",main="Weights for the first IC",xlab = "Variables",ylab="Coefficients")
abline(h=0)
text(1:ncol(data_credit_for_pca),FICA.data_credit$K[,1],labels=colnames(data_credit_for_pca),pos=1,col="firebrick4",cex=0.5)

plot(1:ncol(data_credit_for_pca),FICA.data_credit$K[,2],pch=20,col="deepskyblue4",main="Weights for the second IC",xlab = "Variables",ylab="Coefficients")
abline(h=0)
text(1:ncol(data_credit_for_pca),FICA.data_credit$K[,2],labels=colnames(data_credit_for_pca),pos=1,col="firebrick4",cex=0.5)



# Compare it with the PCA weights
PCS.data_credit$rotation[,1:2]

# The ICA scores can be found here
head(FICA.data_credit$S)

# Have a look at the ICA scores obtained
pairs(FICA.data_credit$S,col=colors.credit,main="Four ICA components")

# Note that the ICs scores have zero mean and identity covariance matrix
colMeans(FICA.data_credit$S)
cov(FICA.data_credit$S)

# Compare PCA scores with ICA scores
par(mfrow=c(2,2))
plot(PCS.data_credit$x[,1],FICA.data_credit$S[,1],pch=19,xlab="PC",ylab="IC",main="First")
plot(PCS.data_credit$x[,2],FICA.data_credit$S[,2],pch=19,xlab="PC",ylab="IC",main="Second")
plot(PCS.data_credit$x[,3],FICA.data_credit$S[,3],pch=19,xlab="PC",ylab="IC",main="Third")
plot(PCS.data_credit$x[,4],FICA.data_credit$S[,4],pch=19,xlab="PC",ylab="IC",main="Fourth")

# The PCs and the ICs are quite different



