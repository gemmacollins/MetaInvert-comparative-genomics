library(tidyverse)
library(ggbiplot)
library(GGally)
library(taxize)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(Manu)

##install Manu (NZ bird) colour palette
#devtools::install_github("G-Thomson/Manu")
#library(devtools)
#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")

rm(list = ls())
#
# Read data and create bag2 from the tsv file ####
seq_success <- read_tsv(file = "Report_Full.tsv")

names(seq_success) <- c("libid", "batch", "phylum", "countgroup", "taxon", "taxid","sample_number", 
                        "read_length", "reads_max", "Reads median seq:", 
                        "Reads min seq:", "Reads N 50:", "total_reads", "total_bases",
                        "n_contigs", "largest_contig", "assembly_length", "GC", 
                        "N50", "N75", "L50", "L75", "bag", "human_reads",
                        "busco_lineage", "Busco_C", "Busco_S", "Busco_D", "Busco_F",
                        "Busco_M", "Busco_set_size", "mode_cov_dist", "mapped_bases", 
                        "genome_size", "COI_Order", "COI_Genus", "COI_Species")

# Select only bag2, no controls, otherwise unfiltered
seq_success %>%
  filter(bag == "bag2",
         !(phylum %in% c("Chordata", "CONTROL"))) ->
  bag2
  
# Fix Acari ####
bag2 %>%
  filter(countgroup == "Acari") %>%
  group_by(taxon) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(genus = str_split_fixed(taxon, " ", 2)[,1]) ->
  Acari

Acari_orders <- tax_name(query = Acari$genus, get = c("suborder","order"), db = "itis")

Acari %>%
  left_join(Acari_orders, by = c("genus" = "query")) ->
  Acari

#Eupelops 
Acari[12,5] <- c("Sarcoptiformes")
Acari[12,6] <- c("Oribatida")

#Hermannia
Acari[33,5] <- c("Sarcoptiformes")
Acari[33,6] <- c("Oribatida")

#Nothrus
Acari[41:44,5] <- c("Sarcoptiformes")
Acari[41:44,6] <- c("Oribatida")

# Not oribatid: Tyrophagus, Platinothrus, Hypoaspis 
save(file = "Acari_orders_all.Rdata", Acari)

# Add fixed Acari back to bag2 and save as bag2.Rdata ####
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
save(file = "bag2.Rdata", bag2)
# Tidy up data #####
load(file = "bag2.Rdata")

bag2 %>%
  dplyr::select(libid, batch, phylum, countgroup_2, taxon, taxid, Busco_C, Busco_F,#read_length, total_bases, 
                total_reads, n_contigs, assembly_length, GC,
                N50, human_reads, Busco_M, genome_size) %>%
  #filter(total_reads > 1e+06)%>%
  filter_all(all_vars(!is.infinite(.))) %>%
  na.omit() -> 
  bag2_selection

dim(bag2_selection)

#write.table(bag2_selection, "bag2_selection.csv", sep = ",", col.names=NA)

## Manu colour palette (NZ native birds) ####
names(manu_palettes)
colourmanu <- get_pal("Hoiho")
print_pal(colourmanu)
#"Hihi" "Hoiho" "Kaka" "Kakapo" "Kakariki" "Kea" "Kereru" "Kereru_orig" 
#"Korimako" Korora" "Kotare" "Putangitangi" "Takahe" "Takapu" "Titipounamu"
#"Tui" "Pepetuna" "Pohutukawa" "Gloomy_Nudi"

#scale_color_manual(values = colorRampPalette(Pohutukawa)(16)) +
#scale_color_gradient(low=get_pal("Takahe")[1], high = get_pal("Takahe")[2])+

# Explore the data further with some basic visual plots ####
## Plot relationships between two variables
ggplot(bag2_selection) + geom_point(aes(x=total_reads, y=n_contigs, col=countgroup_2), size = 2)+
  scale_color_manual(values = colorRampPalette(colourmanu)(9)) +
  ggtitle("bag2_25-02-2021")

## Subset bag2_selected_omit by taxonomic group
unique(bag2_selection$countgroup_2)

# "Collembola" "Oribatida" "Enchytraeidae" "Acari" "Myriapoda" 
# "Lumbricina" "Lumbricidae"  "Nematoda"  "Diplopoda" "Chilopoda"

## Taxa removed due to NAs "Gamasina"  "Tardigrada" "Symphyla" "Pauropoda" "Isopoda" "Diplura"

Taxon <- "Collembola"

ggplot(subset(bag2_selection, countgroup_2==Taxon)) + 
  geom_point(aes(x=total_reads, y=n_contigs, col=Busco_M), size = 2)+
  scale_color_gradient2(mid=colourmanu[1], high = colourmanu[4])+
  ggtitle("bag2_25-02-2021")

ggplot(subset(bag2_seleciton_omit, countgroup_2==Taxon)) + 
  geom_point(aes(x=log(N50), y=Busco_C, col=assembly_length), size = 2)+
  scale_color_gradient2(mid = colourmanu[4], high = colourmanu[1])+
  ggtitle("bag2_25-02-2021")


ggplot(bag2_selection, aes(x = total_reads, y = n_contigs, col=Busco_M))+
  scale_color_gradient2(mid = colourmanu[1], high = colourmanu[4])+
  geom_point(size=2)+
  facet_wrap(~countgroup_2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ggtitle("bag2_25-02-2021")

## to do a box plot
ggplot(bag2_selection, aes(x=countgroup_2, y=log(N50), fill=countgroup_2))+
  geom_boxplot()+
  scale_fill_manual(values = colorRampPalette(colourmanu)(9)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


# Simple plots for each parameter ####
bag2_selection %>%
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

# par(mfrow=c(2,2))
bag2_selection %>% 
  ggplot(mapping = aes(x = total_reads)) + 
  ggtitle("total_reads")+
  geom_histogram(bins=100) ->
  total_reads
total_reads

pdf(file = "assembly_statistics_plot_variables.pdf")
total_reads
dev.off()

# Generate numerical matrix pca_mx) ####
dim(bag2_selection)
pca_mx <- matrix(unlist(bag2_selection), ncol=16)
rownames(pca_mx) <- bag2_selection$taxon
colnames(pca_mx) <- colnames(bag2_selection)

pca_mx <- pca_mx[,-c(1:6)] # exclude columns with character entries
class(pca_mx)<-"numeric" # specify data as numeric (required for pca)

Busco_CF <- rowSums(pca_mx[, c("Busco_C","Busco_F")]) # sum Busco_C + Busco_F
pca_mx <- cbind(pca_mx, Busco_CF) # add the Busco_CF column 
head(pca_mx) # visible check

#write.table(pca_mx, "pca_mx.csv", sep = ",", col.names=NA)

# PCA analyses ####
pca <- prcomp(pca_mx, center = TRUE, scale = TRUE)

selected_colours <- get_pal("Hihi")
ggbiplot(pca, obs.scale=1, var.scale = 1, circle=TRUE, ellipse = F, groups = rownames(pca_mx))+
  scale_color_manual(values = colorRampPalette(selected_colours)(11)) +
  ggtitle("Updated table 24-02-2021, 10 variables, NAs omitted")

summary(pca) #Get summary stats
str(pca) 
plot(pca) #Screeplot for number of components
pca #Get standard deviations and rotation

#install.packages("factoextra")
library(factoextra)
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(pca,
             col.var = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_biplot(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)