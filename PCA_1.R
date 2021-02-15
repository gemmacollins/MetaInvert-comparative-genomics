library(tidyverse)
library(ggbiplot)

rm(list = ls())

# Read data ####
seq_success <- read_tsv(file = "Report_Full.tsv")

names(seq_success) <- c("libid", "batch", "phylum", "countgroup", "taxon", "taxid", 
                        "read_length", "reads_max", "Reads median seq:", 
                        "Reads min seq:", "Reads N 50:", "total_reads", "total_bases",
                        "n_contigs", "largest_contig", "assembly_length", "GC", 
                        "N50", "N75", "L50", "L75", "bag", "human_reads",
                        "busco_lineage", "Busco_C", "Busco_S", "Busco_D", "Busco_F",
                        "Busco_M", "Busco_set_size")

# Select only bag2, no controls, otherwise unfiltered
seq_success %>%
  filter(bag == "bag2",
         !(phylum %in% c("Chordata", "CONTROL"))) ->
  bag_2

#select columns of interest
bag_2 %>%
  dplyr::select(libid, batch, phylum, countgroup, taxon, taxid, # read_length, total_reads, total_bases, 
                total_reads, n_contigs, assembly_length,
                N50, human_reads, Busco_M) ->
  params_to_plot_GC

# head(params_to_plot_GC)
# typeof(params_to_plot_GC)
#dim(params_to_plot_GC)
#unique(params_to_plot_GC$countgroup)

# Generate numerical matrix (pca_mx) ####
pca_mx <- matrix(unlist(params_to_plot_GC), ncol=12)
# head(pca_mx)
# typeof(pca_mx)
rownames(pca_mx) <- params_to_plot_GC$countgroup
colnames(pca_mx) <- colnames(params_to_plot_GC)
# head(pca_mx)
pca_mx <- pca_mx[,-c(1:6)] 
# head(pca_mx)
class(pca_mx)<-"numeric"

# Perform the PCA analyses ####
## Deal with NAs ####
sum(is.na(pca_mx))
pc_omit <- prcomp(na.omit(pca_mx), center = TRUE, scale = TRUE)
sum(is.na(pc_omit))

## Results ####
ggbiplot(pc_omit, obs.scale=1, var.scale = 1)+
  ggtitle("PCA_1: 6 variables, NAs omitted")

summary(pc_omit) #Get summary stats
#str(pc_omit) 
#plot(pc_omit) #Screeplot for number of components
#pc_omit #Get standard deviations and rotation

