library(tidyverse)
library(ggbiplot)
 # library(GGally)
 # library(taxize)
# library(ggplot2)
# library(ggthemes)
# library(dplyr)
library(Manu)

# Read data and create bag2 object from the tsv file ####
full_view <- read_tsv(file = "full_view.tsv")

full_view %>%
  filter(bag == "bag2",
         !(phylum %in% c("Chordata", "CONTROL")))%>%
  select(phylum, countgroup, taxon, co1_Species, co1_Order, co1_Class, co1_length, 
         reads_avg_seq, reads_total_n, perc_reads_human, 
         n_contigs, largest_contig, total_length, GC_perc, N50, N75, L50, L75, 
         busco_perc_C, busco_perc_S,busco_perc_D, busco_perc_F, busco_perc_M, 
         mapped_bases, genome_size_est) %>%
  # filter(total_length > 1000000) %>%
  # filter(reads_total_n > 1000000)%>%
  # filter(n_contigs > 1000)%>%
  #filter(N50 > 1000) %>%
  mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) ->
  bag2_selection

write.table(bag2_selection, "bag2_selection.csv", sep = ",", col.names=NA)
sum(is.na(bag2_selection))
names(which(colSums(is.na(bag2_selection)) > 0))

busco_perc_CF <- rowSums(bag2_selection[, c("busco_perc_C","busco_perc_F")]) # sum Busco_C + Busco_F
bag2_selection <- cbind(bag2_selection, busco_perc_CF) # add the Busco_CF column

## Manu colour palette (NZ native birds) ####
names(manu_palettes)
colourmanu <- get_pal("Takahe")
print_pal(colourmanu)
#"Hihi" "Hoiho" "Kaka" "Kakapo" "Kakariki" "Kea" "Kereru" "Kereru_orig" 
#"Korimako" Korora" "Kotare" "Putangitangi" "Takahe" "Takapu" "Titipounamu"
#"Tui" "Pepetuna" "Pohutukawa" "Gloomy_Nudi"
#scale_color_manual(values = colorRampPalette(Pohutukawa)(16)) +
#scale_color_gradient(low=get_pal("Takahe")[1], high = get_pal("Takahe")[2])+


# Explore the data further with some basic visual plots ####
## Two-variable comparisons

ggplot(bag2_selection)+
  geom_point(aes(x=log(N50), y=busco_perc_C, colour = countgroup), size = 2)+
  #geom_point(aes(x=log(N50), y=Busco_C, col=countgroup), size = 2)+
  #scale_color_gradient2(high=colourmanu[2], mid = colourmanu[1])+
   # scale_x_log10()+
   # scale_y_log10()+
  #scale_color_manual(values = colorRampPalette(get_pal("Takahe"))(3))+
  ggtitle("bag2_12-03-2021")
 #scale_x_continuous(breaks=,,1000000)+
 #scale_y_continuous(breaks=,,1000)
#  expand_limits(x=0,y=0)+


## Subset by taxonomic group
ggplot(bag2_selection) + 
  geom_point(aes(x=reads_total_n, y=N50, col=busco_perc_M), size = 2)+
  scale_x_log10()+
  scale_y_log10()+
  # scale_x_log10(breaks=,,1000000)+
  # scale_y_log10(breaks=,,1000)+
  facet_wrap(~countgroup)+
  scale_color_gradient(low=get_pal("Takahe")[2], high = get_pal("Takahe")[1])+
  ggtitle("bag2_11-03-2021")

Taxon <- "Astigmata"

ggplot(subset(bag2_selection, countgroup==Taxon)) + 
  geom_point(aes(x=reads_total_n, y=n_contigs, col=busco_perc_M<75), size = 2)+
  #scale_color_gradient2(mid=colourmanu[1], high = colourmanu[4])+
  #scale_x_continuous(breaks=,,16000000)+
  ggtitle("Astigmata, 11-03-2021")

ggplot(subset(bag2_selection, countgroup==Taxon)) + 
  geom_point(aes(x=reads_total_n, y=N50, col=busco_perc_M<75), size = 2)+
  #scale_color_gradient2(mid=colourmanu[1], high = colourmanu[4])+
  #scale_x_continuous(breaks=,,16000000)+
  ggtitle("Astigmata, 11-03-2021")


## Box plot
ggplot(bag2_selection, aes(x=countgroup, y=total_length, fill=countgroup))+
  geom_point()+
  geom_boxplot()+
  scale_fill_manual(values = colorRampPalette(colourmanu)(15)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_y_log10()+
  expand_limits(y = 0)
  ggtitle("bag2_12-03-2021")


# Simple plots for each parameter ####
bag2_selection %>%
  pivot_longer(reads_total_n:Busco_M, names_to = "parameter") %>%
  ggplot(mapping = aes(x = value + 1)) +
  geom_histogram(bins=50) + 
  scale_x_log10() + 
  scale_y_log10() +
  facet_grid(countgroup ~ parameter, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 0)) ->
  parameter_distribution
parameter_distribution


bag2_selection %>% 
  ggplot(mapping = aes(x = (n_contigs))) + 
  ggtitle("reads_total_n")+
  # scale_x_log10()+
  # scale_y_log10()+
  # annotation_logticks()+
  geom_histogram() ->
  n_contigs
n_contigs

# PCA analyses ####
## 1 Prepare a numerical matrix (pca_mx) with no NAs
bag2_selection %>%
  select(phylum, countgroup, reads_total_n, perc_reads_human, n_contigs, 
         total_length, GC_perc, N50, busco_perc_M) %>%
  mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))%>%
  na.omit()->
  bag2_selection_omit
# 
# sum(is.na(bag2_selection_omit))
# names(which(colSums(is.na(bag2_selection_omit)) > 0))
# bag2_selection_omit <- na.omit(bag2_selection_omit)
# pca_phylum <- bag2_selection_omit$phylum
# pca_countgroup <- bag2_selection_omit$countgroup

bag2_selection_omit2 <- subset(bag2_selection_omit, select=-c(phylum, countgroup))
dim(bag2_selection_omit2)
pca_mx <- matrix(unlist(bag2_selection_omit2), ncol=7)

rownames(pca_mx) <- bag2_selection_omit$phylum # or change this to countgroup/phylum
#rownames(pca_mx) <- 1:235 #need unique rownames for fviz_pca below
colnames(pca_mx) <- colnames(bag2_selection_omit2)
class(pca_mx) <-"numeric" 

## 2 Perform the pca analysis and plot
pca <- prcomp(pca_mx, center = TRUE, scale = TRUE)

n_distinct(rownames(pca_mx))

ggbiplot(pca, obs.scale = 1, var.scale = 1, circle = F, ellipse = F, groups = rownames(pca_mx))+
  scale_color_manual(values = colorRampPalette(get_pal("Takahe"))(3)) +
  #theme_classic()+
  ggtitle("12-03-2021_filtered")

## 3 Get PCA summary stats and contributions to loading vars
summary(pca) #Get summary stats
str(pca) 
plot(pca) #Screeplot for number of components
pca #Get standard deviations and rotation
var_explained <- pca$sdev^2/sum(pca$sdev^2)
pca$x %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(size=2) +
  theme_bw(base_size=32) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")


library(factoextra)
fviz_contrib(pca, choice = "var", axes = 1, top = 10) # Contributions of variables to PC1
fviz_contrib(pca, choice = "var", axes = 2, top = 10) # Contributions of variables to PC2
fviz_contrib(pca, choice = "var", axes = 1:2, top = 10) # Contributions of variables to both PCs

# loading_scores <- pca$rotation[,1] #PCA1
# loading_scores <- pca$rotation[,2] #PCA2
# var_scores <- abs(loading_scores) ## get the magnitudes
# var_scores_ranked <- sort(var_scores, decreasing=TRUE)
# top_3_variables <- names(var_scores_ranked[1:3])
# top_3_variables ## show the names of the top 3 variables (just contributing to one PC)

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



#rm(list = ls())
# par(mfrow=c(2,2))
#pdf(file = "assembly_statistics_plot_variables.pdf")
#reads_total_n
#dev.off()
#write.table(bag2_selection, "bag2_selection.csv", sep = ",", col.names=NA)
# save(file = "bag2.Rdata", bag2)
#load(file = "bag2.Rdata")
