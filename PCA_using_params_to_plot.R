# pacman::p_load(pacman, tidyverse, GGally, taxize) 
# install.packages("tidyverse")
# install.packages("devtools")
# library(devtools)
# install_github("vqv/ggbiplot")
library(ggbiplot)
library(GGally)
rm(list = ls())

load(file = "params_to_plot.Rdata")
head(params_to_plot)
typeof(params_to_plot)
dim(params_to_plot)
unique(params_to_plot$countgroup_2)

# Generate numerical matrix (pca_mx) ####
pca_mx <- matrix(unlist(params_to_plot), ncol=12)
head(pca_mx)
typeof(pca_mx)
rownames(pca_mx) <- params_to_plot$countgroup_2
colnames(pca_mx) <- colnames(params_to_plot)
head(pca_mx)
pca_mx <- pca_mx[,-c(1:6)]
head(pca_mx)
class(pca_mx)<-"numeric"

# Perform the PCA analyses ####
## Deal with NAs ####
sum(is.na(pca_mx))

# Columns 5 and 6 have NA values and there are several ways to deal with this:
### 1 run the pca analysis on only the first 4 columns (with no missing data) ####
pc_first4 <- prcomp(pca_mx[,c(1,2,3,4)], center = TRUE, scale = TRUE)

#can graph this using the countgroup_2 categories (from params_to_plot) because the number of rows is the same 
ggbiplot(pc_first4)
ggbiplot(pc_first4, obs.scale=1, var.scale = 1, groups = params_to_plot$countgroup_2, ellipse = TRUE, circle = TRUE)+
  ggtitle("all entries from params_to_plot, first 4 variables")

### 2 replace NAs with 0 ####
# pca_mx0 <- pca_mx
# pca_mx0[is.na(pca_mx0)] <- 0
# pc_0 <- prcomp(pca_mx0, center = TRUE, scale = TRUE)

### 3 run an imputation to replace NAs with mean values from complete data (but this may be no good - differences between taxnomic groups!) ####
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("pcaMethods")
#browseVignettes("pcaMethods")
# library(pcaMethods)
# 
# pc <- pca(pca_mx, method="ppca",center = TRUE)
# imputed <- completeObs(pc)
# head(imputed)

# #write table to file so we can see what the NAs were replaced with 
# library(MASS)
# write.table(imputed, file = "imputed-matrix.csv", sep = ",", row.names = TRUE)

##run a new pca analysis on the imputed data matrix
#pc_imputed <- prcomp(imputed, center = TRUE, scale = TRUE)

##plot this with groups coloured
#ggbiplot(pc_imputed)
#ggbiplot(pc_imputed, obs.scale=1, var.scale = 1, groups = params_to_plot$countgroup_2, ellipse = TRUE, circle = TRUE)+
#   ggtitle("all")

### 4 omit cases with missing data (NAs) ####
pc_omit <- prcomp(na.omit(pca_mx), center = TRUE, scale = TRUE)


# Results ####
#So, we have 4 possible pca analyses that we can plot/find out more. NAs omitted is best for now: 
ggbiplot(pc_omit, obs.scale=1, var.scale = 1)+
  ggtitle("6 variables, NAs omitted")

summary(pc_omit) #Get summary stats
str(pc_omit) 
plot(pc_omit) #Screeplot for number of components
pc_omit #Get standard deviations and rotation

predict(pc_omit) %>% round(2) # See how cases load on PCs
biplot(pc_omit) # Biplot of first two components

## plot that shows the PCs and the variation:
plot(pc_omit$x[,1], pc_omit$x[,2])

#calculate how much variation in the original data the pcs count for
pca.var <- pc_omit$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

library(ggplot2)
pca.data <- data.frame(Sample=rownames(pc_omit$x),
                       X=pc_omit$x[,1],
                       Y=pc_omit$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))
#+ggtitle("")

#to determine which variables have the largest effect on where samples are plotted
loading_scores <- pc_omit$rotation[,1]
var_scores <- abs(loading_scores) ## get the magnitudes
var_scores_ranked <- sort(var_scores, decreasing=TRUE)
top_3_variables <- names(var_scores_ranked[1:3])
top_3_variables ## show the names of the top 3 variables

####Stop here for now ####
##generate data frame (so taxa name can be grouping variable)
pca_df <- data.frame(matrix(unlist(params_to_plot), ncol=12))
head(pca_df)
colnames(pca_df) <- colnames(params_to_plot)
head(pca_df)

pca_df <- pca_df[,-c(1:3, 5,6)]

class(pca_df) <- "numeric"
pc <- prcomp(pca_df[,c(1:4)], center = TRUE, scale = TRUE)


##group the rows by the taxonomic group (x9)
by_taxgroup <- params_to_plot %>% group_by(countgroup_2)
by_taxgroup
ungroup(params_to_plot)

##now, try to make the matrix again but only with one group at a time

# install.packages("dplyr")
# library(dplyr)
df %>% group_by()

install.packages("plyr")
library(plyr)

split(params_to_plot, params_to_plot$countgroup_2)
params_to_plot$Tardigrada
