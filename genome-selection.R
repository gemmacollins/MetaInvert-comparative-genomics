library(tidyverse)
library(ggbiplot)
library(GGally)
library(taxize)

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
  bag2
  
# Fix Acari ####
# bag2 %>%
#   filter(countgroup == "Acari") %>%
#   group_by(taxon) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n)) %>%
#   mutate(genus = str_split_fixed(taxon, " ", 2)[,1]) ->
#   Acari
# 
# Acari_orders <- tax_name(query = Acari$genus, get = c("suborder","order"), db = "itis")
# 
# Acari %>%
#   left_join(Acari_orders, by = c("genus" = "query")) ->
#   Acari
# 
# Acari[8,5] <- c("Sarcoptiformes")
# Acari[8,6] <- c("Oribatida")
# Acari[13,5] <- c("Sarcoptiformes")
# Acari[13,6] <- c("Oribatida")
# 
# save(file = "Acari_orders_all.Rdata", Acari)
load(file = "Acari_orders_all.Rdata")

Acari %>%
  filter(!is.na(suborder)) ->
  oribatids

bag2 %>%
  mutate(countgroup_2 = case_when(
    taxon %in% oribatids$taxon ~ "Oribatida",
    TRUE ~ countgroup
  )) ->
  bag2

rm(Acari, oribatids) # Clean up
write.table(bag2, "bag2.csv", sep = ",", col.names=NA)
#save(file = "bag2.Rdata", bag2)
load(file = "bag2.Rdata")

# Explore the data further with some basic visual plots ####
# install.packages("Rtools","ggthemes")
# library(readxl)
library(ggplot2)
library(ggthemes)

datafile <- bag2

ggplot(datafile) + geom_point(aes(x=total_reads, y=n_contigs, col=Busco_M))+
  ggtitle("bag2, n = 412-104NAs = 308")

ggplot(bag2, aes(x = total_reads, y = n_contigs, col=Busco_M))+
  geom_point()+
  facet_wrap(~countgroup_2)+
  theme_bw()+
  ggtitle("bag2, n = 412-104NAs = 308")



#to do a box plot:
# ggplot(datafile, aes(x=countgroup_2, y=total_reads, fill=countgroup_2,))+
#   geom_boxplot()+
#   stat_boxplot(geom="errorbar")

# Select columns of interest
bag2 %>%
  dplyr::select(libid, batch, phylum, countgroup_2, taxon, taxid, # read_length, total_bases, 
                total_reads, n_contigs, assembly_length,
                N50, human_reads, Busco_M, Busco_C, Busco_F) ->
  params_to_plot_GC

params_to_plot_GC %>%
  pivot_longer(total_reads:Busco_M, names_to = "parameter") %>%
  ggplot(mapping = aes(x = value + 1)) +
  geom_histogram(bins=50) + 
  scale_x_log10() + 
  scale_y_log10() +
  facet_grid(countgroup_2 ~ parameter, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 0)) ->
  parameter_distribution
parameter_distribution

# pdf('plots.pdf')
# par(mfrow=c(2,2))
params_to_plot_GC %>% 
  ggplot(mapping = aes(x = total_reads)) + 
  ggtitle("total_reads")+
  geom_histogram(bins=100) ->
  total_reads
total_reads

params_to_plot_GC %>% 
  ggplot(mapping = aes(x = n_contigs)) + 
  ggtitle("n_contigs")+
  geom_histogram(bins = 100) ->
  n_contigs
n_contigs

params_to_plot_GC %>% 
  ggplot(mapping = aes(x = assembly_length)) + 
  ggtitle("assembly_length")+
  geom_histogram(bins=100) ->
  assembly_length
assembly_length

params_to_plot_GC %>% 
  ggplot(mapping = aes(x = N50)) + 
  ggtitle("N50")+
  geom_histogram(bins=100) ->
  N50
N50

params_to_plot_GC %>% 
  ggplot(mapping = aes(x = human_reads)) + 
  ggtitle("human_reads")+
  geom_histogram(bins=100) ->
  human_reads
human_reads

params_to_plot_GC %>% 
  ggplot(mapping = aes(x = Busco_M)) + 
  ggtitle("Busco_M")+
  geom_histogram(bins=100) ->
  Busco_M
Busco_M

pdf(file = "assembly_statistics_plot_variables.pdf")
total_reads
n_contigs
assembly_length
N50
human_reads
Busco_M
dev.off()



# Generate numerical matrix (pca_mx) ####
dim(params_to_plot_GC)
pca_mx <- matrix(unlist(params_to_plot_GC), ncol=14)
#rownames(pca_mx) <- params_to_plot_GC$countgroup_2
colnames(pca_mx) <- colnames(params_to_plot_GC)

pca_mx <- pca_mx[,-c(1:6)] # exclude columns with character entries
class(pca_mx)<-"numeric" # specify data as numeric (required for pca)

Busco_CF <- rowSums(pca_mx[, c(7,8)]) # sum Busco_C + Busco_F
pca_mx <- cbind(pca_mx, Busco_CF) # add the Busco_CF column 
head(pca_mx) # visible check

#write.table(pca_mx, "pca_mx.csv", sep = ",", col.names=NA)

# Perform the PCA analyses ####
## Deal with NAs ####
sum(is.na(pca_mx))
pc_omit <- prcomp(na.omit(pca_mx), center = TRUE, scale = TRUE)
sum(is.na(pc_omit))

## Results ####
ggbiplot(pc_omit, obs.scale=1, var.scale = 1)+
  ggtitle("PCA2: bag2_improved, 9 variables, NAs omitted")

summary(pc_omit) #Get summary stats
str(pc_omit) 
plot(pc_omit) #Screeplot for number of components
pc_omit #Get standard deviations and rotation


install.packages("factoextra")
library(factoextra)
fviz_pca_ind(pc_omit,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(pc_omit,
             col.var = "black", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_biplot(pc_omit,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)






# Manually removed some rows in excel then imported back here ####
install.packages("readxl")
library(readxl)
bag2_improved <- read_excel("bag2_improved.xlsx")

datafile <- bag2_improved

bag2_improved %>%
  dplyr::select(rownumber, libid, batch, phylum, countgroup_2, taxon, taxid, # read_length, total_bases, 
                total_reads, n_contigs, assembly_length,
                N50, human_reads, Busco_M, Busco_C, Busco_F) ->
  params_to_plot_GC

dim(params_to_plot_GC)
pca_mx <- matrix(unlist(params_to_plot_GC), ncol=15)
rownames(pca_mx) <- params_to_plot_GC$rownumber
colnames(pca_mx) <- colnames(params_to_plot_GC)
head(pca_mx)
pca_mx <- pca_mx[,-c(1:7)] # exclude columns with character entries
class(pca_mx)<-"numeric" # specify data as numeric (required for pca)

Busco_CF <- rowSums(pca_mx[, c(7,8)]) # sum Busco_C + Busco_F
pca_mx <- cbind(pca_mx, Busco_CF) # add the Busco_CF column 
head(pca_mx) # visible check

sum(is.na(pca_mx))
pc_omit <- prcomp(na.omit(pca_mx), center = TRUE, scale = TRUE)
sum(is.na(pc_omit))

## Results ####
ggbiplot(pc_omit, obs.scale=1, var.scale = 1)+
  ggtitle("PCA2: bag2_improved, 9 variables, NAs omitted")


