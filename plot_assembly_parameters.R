library(tidyverse)
library(GGally)
library(taxize)

rm(list = ls())

# Read data ####
seq_success <- read_tsv(file = "Report_Full.tsv")
seq_success$batch <- factor(seq_success$batch, 
                            levels = c("b0","b1", "b2", "b3_b5", "b4", 
                                       "b7", "b8", "b9", "b10", "b11"))

names(seq_success) <- c("libid", "batch", "phylum", "countgroup", "taxon", "taxid", 
                        "read_length", "reads_max", "Reads median seq:", 
                        "Reads min seq:", "Reads N 50:", "total_reads", "total_bases",
                        "n_contigs", "largest_contig", "assembly_length", "GC", 
                        "N50", "N75", "L50", "L75", "bag", "human_reads",
                        "busco_lineage", "Busco_C", "Busco_S", "Busco_D", "Busco_F",
                        "Busco_M", "Busco_set_size")

# Select genomes ####
# bag2, no controls

seq_success %>% 
  ggplot(mapping = aes(x = total_reads)) + 
  geom_histogram()  
  # scale_x_log10() +
  # scale_y_log10()

seq_success %>%
  filter(bag == "bag2",
         !(phylum %in% c("Chordata", "CONTROL")),
         total_reads > 5e+06) ->
  bag2

# seq_success %>%
#   filter(bag == "bag2",
#          !(phylum %in% c("Chordata", "CONTROL")),
#          total_reads > 5e+06) ->
#   deeper_sequenced

# Fix Acari ####
# deeper_sequenced %>%
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
# save(file = "Acari_orders.Rdata", Acari)
load(file = "Acari_orders.Rdata")

Acari %>%
  filter(!is.na(suborder)) ->
  oribatids

deeper_sequenced %>%
  mutate(countgroup_2 = case_when(
    taxon %in% oribatids$taxon ~ "Oribatida",
    TRUE ~ countgroup
  )) ->
  deeper_sequenced

# Clean up
rm(Acari, oribatids)

# Plot ####
very_few <- c("Lumbricidae", "Lumbricina", "Symphyla", "Pauropoda", "Diplura", "Acari")
deeper_sequenced %>%
  dplyr::select(libid, batch, phylum, countgroup_2, taxon, taxid, # read_length, total_reads, total_bases, 
                total_reads, n_contigs, assembly_length,
                N50, human_reads, Busco_M) %>%
  filter(!(countgroup_2 %in% very_few)) ->
  params_to_plot

# Histograms
params_to_plot %>%
  pivot_longer(total_reads:Busco_M, names_to = "parameter") %>%
  ggplot(mapping = aes(x = value + 1)) +
  geom_histogram() + 
  scale_x_log10() + 
  scale_y_log10() +
  facet_grid(countgroup_2 ~ parameter, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 0)) ->
  parameter_distribution
parameter_distribution

# Busco vs assembly
params_to_plot %>%
  ggplot(mapping = aes(x = assembly_length, y = Busco_M, color = log(N50))) +
  geom_smooth() +
  geom_point(size = 1) +
  ylim(0,100) +
  facet_wrap(~countgroup_2) +
  scale_color_gradient(low="orange", high="darkblue", space ="Lab") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Assembly length VS missing BUSCOs") ->
  busco_assembly_length
busco_assembly_length

# N50 - BUSCO
params_to_plot %>%
  ggplot(mapping = aes(x = N50, y = Busco_M, color = assembly_length)) +
  geom_smooth() +
  geom_point() +
  ylim(0,100) +
  scale_x_log10() +
  facet_wrap(~countgroup_2) +
  scale_color_gradient(low="orange", high="darkblue", space ="Lab") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("N50 VS missing BUSCOs") ->
  busco_N50
busco_N50

# N contigs VS N50
params_to_plot %>%
  ggplot(mapping = aes(x = n_contigs, y = N50, color = Busco_M)) +
  geom_smooth() +
  geom_point() +
  # ylim(0,100) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~countgroup_2) +
  scale_color_gradient(low="orange", high="darkblue", space ="Lab") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Number of contigs VS N50") ->
  n_contigs_N50
n_contigs_N50

# Assembly length VS N50 VS BUSCO
params_to_plot %>%
  ggplot(mapping = aes(x = assembly_length, y = N50, color = Busco_M)) +
  geom_smooth() + # method = "gam", formula = y ~ s(x, bs = "cs")
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~countgroup_2) +
  scale_color_gradient(low="orange", high="darkblue", space ="Lab") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Assembly length VS N50") ->
  assembly_length_N50
assembly_length_N50


pdf(file = "assembly_statistics_plots.pdf")
parameter_distribution
busco_assembly_length
busco_N50
n_contigs_N50
assembly_length_N50
dev.off()




# Busco vs seq. depth

# Busco VS 

# To Do ####
# fix countgrous: Acari -> Oribatida
# filter

# deeper_sequenced$countgroup %>%
#   unique()
# 
# # Count genomes



# oribatid_names <- c("Acrotritia duplicata", )
# other_mites <- c("Tyrophagus putrescentiae", "Allosuctobelba grandis", )
# # 
# deeper_sequenced %>%
#   group_by(countgroup) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))
# 
# seq_success %>%
#   group_by(phylum) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))
# 
# seq_success %>%
#   filter(bag == "bag2",
#          !(phylum %in% c("Aves", "CONTROL")),
#          total_reads < 2.5e+07) %>%
#   group_by(countgroup) %>%
#   summarize(n = n()) %>%
#   arrange(desc(n))
# 
# seq_success %>%
#   filter(bag == "bag2",
#          countgroup == "Acari") %>%
#   ggplot(mapping = aes(x = total_bases)) +
#   geom_histogram() +
#   scale_x_log10()
# 
# seq_success %>%
#   group_by(batch, countgroup) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))
# 
# 
# 
# 
# # bag3
# seq_success %>%
#   filter(bag == "bag3",
#          !(phylum %in% c("Chordata", "CONTROL")),
#          `Reads total seq:` > 2000000000) %>% 
#   dplyr::select(libid, `Total length`) ->
#   bag3
# names(bag3) <- c("libid", "Total length bag3")
# 
# seq_success %>%
#   filter(bag == "bag2",
#          !(phylum %in% c("Chordata", "CONTROL")),
#          `Reads total seq:` > 2000000000) %>% 
#   # arrange(desc(N50)) %>%
#   # group_by(taxon) %>%
#   # slice_head() %>% # keep assembly with the best N50 for a taxon
#   # ungroup %>%
#   dplyr::select(libid, batch, phylum, countgroup, `Reads total seq:`, 
#                 `Largest contig`, `Total length`, `N50`, 
#                 perc_reads_human, `Busco % M`) %>% 
#   left_join(bag3, by = "libid") %>%
#   mutate(Total_NonTarget = `Total length bag3` - `Total length`) %>% 
#   dplyr::select(libid, batch, phylum, countgroup, `Reads total seq:`, 
#                 `Largest contig`, `Total length`, `N50`, 
#                 perc_reads_human, `Busco % M`, Total_NonTarget) ->
#   for_insigts
# 
# names(for_insigts) <- c("libid", "batch", "phylum", "countgroup",
#                         "Bp sequenced", "Largest contig", "Assembly length",
#                         "N50", "Human %", "Missing BUSCO %", "Bp NonTarget")
# 
# # Plots ####
# log10_vars <- c("Bp sequenced", "Largest contig", "Assembly length",
#                 "N50", "Bp NonTarget")
# 
# pdf(file = "pairwise_parameters.pdf", height = 15, width = 15)
# for_insigts %>%
#   dplyr::select(phylum, `Bp sequenced`:`Bp NonTarget`) %>%
#   mutate_at(log10_vars, log10) %>% # convert variables into log10
#   rename_at(log10_vars, sprintf, fmt = "log10 %s") %>%
#   ggpairs(mapping = ggplot2::aes(color = phylum, alpha = 0.8),
#           lower = list(continuous = wrap("smooth", alpha = 0.7, size=0.4))) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# dev.off()
# 
#              # define variables to be transformed
# 
# iris %>%                                             # use standard R example dataframe
#   mutate_at(log10_vars, log10) %>%                   # log10 transform selected columns
#   rename_at(log10_vars, sprintf, fmt="log10 %s") %>% # rename variables accordingly
# GGally::ggpairs(aes(color=Species))
# 
#   
# pdf(file = "soil_genome_summaries.pdf")
# # seq_success %>%
# #   filter(bag == "bag2",
# #          !(phylum %in% c("Chordata", "CONTROL"))) %>%  # more complete bag only, 
# #   ggplot(mapping = aes(y = N50)) +
# #   geom_histogram() +
# #   scale_y_log10() +
# #   facet_grid(~batch) +
# #   theme(axis.text.x = element_text(angle = 90))
# 
# # N50 ~ taxon ####
# seq_success %>%
#   filter(bag == "bag2",
#          !(phylum %in% c("Chordata", "CONTROL")),  
#          `Reads total seq:` > 2000000000) %>% # more complete bag only
#   # 
#   arrange(desc(N50)) %>%
#   group_by(taxon) %>%
#   slice_head() %>% # keep assembly with the best N50 for a taxon
#   ungroup %>%
#   ggplot(mapping = aes(y = N50)) +
#   geom_histogram() +
#   scale_y_log10() +
#   facet_grid(~countgroup) +
#   # theme(axis.text.x = element_text(angle = 90)) +
#   theme(strip.text.x = element_text(angle = 90))
# 
# seq_success %>%
#   filter(bag == "bag2",
#          !(phylum %in% c("Chordata", "CONTROL")),
#          `Reads total seq:` > 2000000000) %>% # more complete bag only
#   arrange(desc(N50)) %>%
#   group_by(taxon) %>%
#   slice_head() %>% # keep assembly with the best N50 for a taxon
#   ungroup %>%
#   ggplot(mapping = aes(y = `Busco % M`)) +
#   geom_histogram() +
#   scale_y_log10() +
#   facet_grid(~countgroup) +
#   theme(strip.text.x = element_text(angle = 90)) +
#   ylab("Missing BUSCO %")
# dev.off()
