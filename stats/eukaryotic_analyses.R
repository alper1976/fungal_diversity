# Statistical analysis of alpha and beta diversity in relation to meta data



################## Work in progress ################################

## saga
#module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1 # in addition all bioconductor packages will be available
#module load GDAL/3.5.0-foss-2022a # this is needed to make maps
#module load MariaDB-connector-c/3.1.7-GCCcore-9.3.0 # this is needed to make maps




.libPaths(c("/cluster/projects/nn9745k/Rpackages_4_2_1", .libPaths()))
norway <- file.path("/cluster/projects/nn9744k/02_results/29_euk_20220325/Pedro/final_tax_classification") # <---------- Change to your path

sweden <- file.path("/cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/final_tax_classification") # <---------- Change to your path

figsPath <- file.path("/cluster/projects/nn9744k/02_results/30_scandinavian_fungi/figs")

# load(file=file.path(figsPath, "stats_eukaryotic.RData"))

# image
# save.image(file=file.path(figsPath, "stats_eukaryotic.RData"))

##load packages
if (!require("xlsx")) {
   install.packages("xlsx", dependencies = TRUE)
   library(xlsx)
   }
if (!require("vegan")) {
   install.packages("vegan", dependencies = TRUE)
   library(vegan)
   }
if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2)
   }
if (!require("ggrepel")) {
   install.packages("ggrepel", dependencies = TRUE)
   library(ggrepel)
   }
if (!require("ape")) {
   install.packages("ape", dependencies = TRUE)
   library(ape)
   }
if (!require("reshape2")) {
   install.packages("reshape2", dependencies = TRUE)
   library(reshape2)
   }
library(phyloseq)

if (!require("Hmisc")) {
   install.packages("Hmisc", dependencies = TRUE)
   library(Hmisc)
   }
if (!require("devtools")) {
   install.packages("devtools", dependencies = TRUE)
   library(devtools)
   }
if (!require("corrplot")) {
  install.packages("corrplot", dependencies = TRUE)
  library(corrplot)
  }
if (!require("gridExtra")){
  install.packages("gridExtra")
  library(gridExtra)
  }
if (!require("dplyr")){
  install.packages("dplyr")
  library(dplyr)
  }
if (!require("mice")){
  install.packages("mice")
  library(mice)
  }
# if (!require("EcoSimR")) {
#    install.packages("EcoSimR", dependencies = TRUE)
#    library(EcoSimR)
#    }
if (!require("Biostrings")) {
   install.packages("Biostrings", dependencies = TRUE)
   library(Biostrings)
   }
   library(dada2)

if (!require("cowplot")) {
   install.packages("cowplot", dependencies = TRUE)
   library(cowplot)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }
if (!require("R.utils")) {
   install.packages("R.utils", dependencies = TRUE)
   library(R.utils)
   }
if (!require("ade4")) {
   install.packages("ade4", dependencies = TRUE)
   library(ade4)
   }
if (!require("stringr")) {
   install.packages("stringr", dependencies = TRUE)
   library(stringr)
   }
if (!require("data.table")) {
   install.packages("data.table", dependencies = TRUE)
   library(data.table)
   }

## Libraries to make map


################ Load local packages #####################################

sourceDirectory('/cluster/projects/nn9745k/scripts/amplicon_analysis_package')

# install_github("GotelliLab/EcoSimR")

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color_palette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

color_palette_class = c("#D1BBD7", "#AE76A3", "#882E72", "#1965B0",
                        "#5289C7", "#7BAFDE", "#4EB265", "#90C987",
                        "#CAE0AB", "#F7F056", "#F6C141", "#F1932D",
                        "#E8601C", "#DC050C", "#72190E")
color_palette_phylum = c("#114477", "#4477AA", "#77AADD", "#117755",
                          "#44AA88", "#99CCBB", "#777711", "#AAAA44",
                          "#DDDD77", "#771111", "#AA4444", "#DD7777",
                          "#771144", "#AA4477", "#DD77AA")


########################## Load metadata ############################

nor_meta = data.frame(read.xlsx(file.path(norway, "100lakes_norway_meta_sample.xlsx"), 1, header=T, row.names = "sample", check.names=FALSE))

summary(nor_meta)

swe_meta = data.frame(read.table(file.path("/cluster/projects/nn9745k/02_results/02_hundred_lakes/lake_metadata", "100lakes_metadata_final.csv"), sep = ";", row.names = "euk_seq_ids", header = TRUE))

swe_meta = swe_meta[,1:73]
summary(swe_meta)

######################### Eukaryotic data analysis #####################################

#NORWAY
euk_nor = read.table(file.path("/cluster/projects/nn9745k/02_results/03_hundred_lakes_norway/eukaryota", "ASV_table.tsv"))
euk_nor = euk_nor[c(6:79),]
euk_nor_names = str_remove(colnames(euk_nor), "([X])")
colnames(euk_nor) = euk_nor_names
sample_names = str_extract(rownames(euk_nor), "[^_]+")
rownames(euk_nor) = sample_names

tax_euk_nor = read.table(file.path("/cluster/projects/nn9745k/02_results/03_hundred_lakes_norway/eukaryota", "Taxonomy_table_pr2.tsv"))
nor_otu_euk <- otu_table(t(euk_nor), taxa_are_rows = TRUE)
dim(nor_otu_euk)
nor_tax_euk <- tax_table(as.matrix(tax_euk_nor))

nor_physeq_euk <- phyloseq(nor_otu_euk, nor_tax_euk, sample_data(data.frame(nor_meta)))

nor_physeq_euk_clean <- subset_taxa(nor_physeq_euk, Kingdom == "Eukaryota")

nor_physeq_euk_clean = prune_samples(sample_sums(nor_physeq_euk_clean)>=1000, nor_physeq_euk_clean)

saveRDS(nor_physeq_euk_clean, file = file.path(figsPath, "nor_physeq_euk_clean.rds"))

taxonomy = "Division"

classGlom = tax_glom(nor_physeq_euk_clean, taxrank = taxonomy)

taxon_table = otu_table(classGlom)
tax_matrix = as(tax_table(classGlom), 'matrix')
rownames(taxon_table) = tax_matrix[,taxonomy]
tax_table = prop.table(taxon_table, margin = 2)*100
tax_table <- tax_table[order(rowSums(-taxon_table)),]
tax_table <- tax_table[1:15,]
rownames(tax_table)
summary(t(tax_table))


melted_tax_table = melt(tax_table)
melted_tax_table <- arrange(melted_tax_table, Var1, desc(value))
melted_tax_table$Var2 <- factor(melted_tax_table$Var2 , levels = unique(melted_tax_table$Var2))

class_colors <- setNames(color_palette_class, levels(melted_tax_table$Var1))

cairo_ps(file.path(figsPath, "nor_eukaryotes.eps"), width=120/25.4, height=160/25.4, pointsize = 6, bg = FALSE, fallback_resolution = 300)
  ggplot(melted_tax_table, aes(x = Var2, y = value, fill = Var1))+
     geom_bar(stat = "identity", position = "stack", alpha = .5) +
     guides(fill = guide_legend(title = taxonomy)) +
     coord_flip() +
     theme(axis.text = element_text(size=6),
       axis.title = element_text(size=10, face="bold"),
       legend.text = element_text(size=8),
       plot.background = element_rect(fill = "white"),
       panel.background = element_rect(fill = "white"),
       axis.line.x = element_line(color = "grey")) +
     xlab("lake systems") +
     ylab("proportion of reads [%]") +
     scale_fill_manual(values = class_colors) +
     scale_x_discrete(limits = rev(levels(melted_tax_table$Var2))
     )
dev.off()

metadata4 = data.frame(sample_data(classGlom))
metadata4 = metadata4[,c("Cond", "TP", "TOC", "pH", "Latitude", "Forest")]

rda_div_meta <- rda(decostand(t(tax_table), method="hellinger"), as.matrix(metadata4), na.rm = na.rm)

rda_div_meta
summary_rda_div_meta <- summary(rda_div_meta)
RsquareAdj(rda_div_meta)

metadata_class = cbind(metadata4, t(taxon_table))

rda_scores_env = vegan::scores(rda_div_meta, display = "bp")
rda_scores_species = vegan::scores(rda_div_meta, display = "sp")

names = rownames(rda_scores_species)
names
melted_tax_table$Var1


rda_plot_species <- ggplot(data.frame(rda_scores_species), aes(x = RDA1, y = RDA2, color = names)) +
  geom_point(size = 7, alpha = .5) +
  scale_color_manual(values = class_colors)

mult = 0.2

rda_biplot_class <- rda_plot_species +
  geom_segment(data = data.frame(rda_scores_env), aes(x = 0, xend = mult * RDA1,
                                                      y = 0, yend = mult * RDA2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey", alpha = .8) +
  geom_text_repel(data = data.frame(rda_scores_env),
            aes(x = mult * RDA1, y = mult * RDA2, label = rownames(rda_scores_env),
                hjust = 0.5 * (1-sign(RDA1)), vjust = 0.5 * (1-sign(RDA2))),
            color = "grey", size = 6, alpha = .8) +
  # coord_cartesian(xlim = c(-0.52, 0.25), ylim = c(-0.3, 0.2)) +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=18, face="bold"),
        legend.text = element_text(size=10),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey"),
        panel.grid = element_line(color = "grey"),
        legend.position = "None")

cairo_ps(file.path(figsPath, "nor_eukaryotic_rda.eps"), width=120/25.4, height=160/25.4, pointsize = 6, bg = FALSE, fallback_resolution = 300)
  rda_biplot_class
dev.off()


#SWEDEN
euk_swe = read.table(file.path("/cluster/projects/nn9745k/02_results/02_hundred_lakes/eukaryota", "ASV_table.tsv"))
euk_swe = euk_swe[c(1,11:100),]
euk_swe_names = str_remove(colnames(euk_swe), "([X])")
colnames(euk_swe) = euk_swe_names

tax_euk_swe = read.table(file.path("/cluster/projects/nn9745k/02_results/02_hundred_lakes/eukaryota", "Taxonomy_table_pr2.tsv"))

swe_otu_euk <- otu_table(t(euk_swe), taxa_are_rows = TRUE)
dim(swe_otu_euk)
swe_tax_euk <- tax_table(as.matrix(tax_euk_swe))

swe_physeq_euk <- phyloseq(swe_otu_euk, swe_tax_euk, sample_data(data.frame(swe_meta)))

swe_physeq_euk_clean <- subset_taxa(swe_physeq_euk, Kingdom == "Eukaryota")

swe_physeq_euk_clean = prune_samples(sample_sums(swe_physeq_euk_clean)>=1000, swe_physeq_euk_clean)
saveRDS(swe_physeq_euk_clean, file = file.path(figsPath, "swe_physeq_euk_clean.rds"))


taxonomy = "Division"

classGlom = tax_glom(swe_physeq_euk_clean, taxrank = taxonomy)

taxon_table = otu_table(classGlom)
tax_matrix = as(tax_table(classGlom), 'matrix')
rownames(taxon_table) = tax_matrix[,taxonomy]
tax_table = prop.table(taxon_table, margin = 2)*100
tax_table <- tax_table[order(rowSums(-taxon_table)),]
rowSums(taxon_table)
dim(taxon_table)
tax_table <- tax_table[1:15,]
rownames(tax_table)
summary(t(tax_table))

melted_tax_table = melt(tax_table)
sample_names = str_extract(melted_tax_table$Var2, "[^_]+")
sample_names = str_remove(sample_names, "(l100e)")
melted_tax_table$Var2 = sample_names
melted_tax_table <- arrange(melted_tax_table, Var1, desc(value))
melted_tax_table$Var2 <- factor(melted_tax_table$Var2 , levels = unique(melted_tax_table$Var2))


class_colors <- setNames(color_palette_class, levels(melted_tax_table$Var1))

cairo_ps(file.path(figsPath, "swe_eukaryotes.eps"), width=120/25.4, height=160/25.4, pointsize = 6, bg = FALSE, fallback_resolution = 300)
  ggplot(melted_tax_table, aes(x = Var2, y = value, fill = Var1))+
     geom_bar(stat = "identity", position = "stack", alpha = .5) +
     guides(fill = guide_legend(title = taxonomy)) +
     coord_flip() +
     theme(axis.text = element_text(size=6),
       axis.title = element_text(size=10, face="bold"),
       legend.text = element_text(size=8),
       plot.background = element_rect(fill = "white"),
       panel.background = element_rect(fill = "white"),
       axis.line.x = element_line(color = "grey")) +
     xlab("lake systems") +
     ylab("proportion of reads [%]") +
     scale_fill_manual(values = class_colors) +
     scale_x_discrete(limits = rev(levels(melted_tax_table$Var2))
     )
dev.off()

metadata4 = data.frame(sample_data(classGlom))
metadata4 = metadata4[,c("Kond_25", "TotP", "TOC", "pH", "Latitude", "Forest")]

rda_div_meta <- rda(decostand(t(tax_table), method="hellinger"), as.matrix(metadata4), na.rm = TRUE)

rda_div_meta
summary_rda_div_meta <- summary(rda_div_meta)
RsquareAdj(rda_div_meta)

metadata_class = cbind(metadata4, t(taxon_table))

rda_scores_env = vegan::scores(rda_div_meta, display = "bp")
rda_scores_species = vegan::scores(rda_div_meta, display = "sp")

names = rownames(rda_scores_species)
names
melted_tax_table$Var1


rda_plot_species <- ggplot(data.frame(rda_scores_species), aes(x = RDA1, y = RDA2, color = names)) +
  geom_point(size = 7, alpha = .5) +
  scale_color_manual(values = class_colors)

mult = 0.2

rda_biplot_class <- rda_plot_species +
  geom_segment(data = data.frame(rda_scores_env), aes(x = 0, xend = mult * RDA1,
                                                      y = 0, yend = mult * RDA2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey", alpha = .8) +
  geom_text_repel(data = data.frame(rda_scores_env),
            aes(x = mult * RDA1, y = mult * RDA2, label = rownames(rda_scores_env),
                hjust = 0.5 * (1-sign(RDA1)), vjust = 0.5 * (1-sign(RDA2))),
            color = "grey", size = 6, alpha = .8) +
  # coord_cartesian(xlim = c(-0.52, 0.25), ylim = c(-0.3, 0.2)) +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=18, face="bold"),
        legend.text = element_text(size=10),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey"),
        panel.grid = element_line(color = "grey"),
        legend.position = "None")

cairo_ps(file.path(figsPath, "swe_eukaryotic_rda.eps"), width=120/25.4, height=160/25.4, pointsize = 6, bg = FALSE, fallback_resolution = 300)
  rda_biplot_class
dev.off()

save.image(file=file.path(figsPath, "stats_eukaryotic.RData"))
