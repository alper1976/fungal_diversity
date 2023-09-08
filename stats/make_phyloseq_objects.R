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

# load(file=file.path(figsPath, "stats.RData"))

# image
# save.image(file=file.path(figsPath, "stats.RData"))

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

############################# 100 Lakes Sweden + Norway ##########################################

##Reading sequence table and taxa table
sweseqtab = as.matrix(read.xlsx(file.path(sweden, "OTU_table_lulu_curated.xlsx"), 1, header=T, row.names = "OTUid",check.names=FALSE))
sweseqtab = sweseqtab [,1:103]
seqtab_names = str_remove(colnames(sweseqtab), "([X])")
colnames(sweseqtab) = seqtab_names

swetaxa = as.matrix(read.xlsx(file.path(sweden, "lca_merged_taxonomy.xlsx"), 1, header=T, row.names = TRUE, check.names=FALSE))
swemeta = data.frame(read.csv(file.path(sweden, "swemeta_trans_interpolated.csv")))
rownames(swemeta) = swemeta$X

norseqtab = as.matrix(read.xlsx(file.path(norway, "OTU_table_lulu_curated_clean.xlsx"), 1, header=T, row.names = "OTUid",check.names=FALSE))

# remove one sample sa31

nortaxa = as.matrix(read.xlsx(file.path(norway, "taxonomy_clean.xlsx"), 1, header=T, row.names = TRUE, check.names=FALSE))
normeta = data.frame(read.csv(file.path(norway, "normeta_trans_interpolated.csv")))
rownames(normeta) = normeta$X
################# Align metadata and sequence data ###################
#norway

norseqtab <- t(norseqtab)
# norseqtab <- norseqtab[!rownames(norseqtab) %in% lakes_to_be_removed, ]

merged_nor <- merge(norseqtab, normeta, by="row.names")
metadata_nor <- merged_nor[, (ncol(norseqtab)+2):ncol(merged_nor)]
rownames(metadata_nor) = merged_nor$X

norseq <- merged_nor[, 2:(ncol(norseqtab)+1)]
rownames(norseq) = merged_nor$X

nor_otu <- otu_table(t(norseq), taxa_are_rows = TRUE)
dim(nor_otu)
nor_tax <- tax_table(nortaxa)
nor_physeq <- phyloseq(nor_otu, nor_tax, sample_data(data.frame(metadata_nor)))
nor_physeq

nor_physeq_clean <- subset_taxa(nor_physeq, kingdom == "Fungi")

nor_physeq_clean = prune_samples(sample_sums(nor_physeq_clean)>=1000, nor_physeq_clean)
saveRDS(nor_physeq_clean, file = file.path(figsPath, "nor_physeq_clean.rds"))

#NORWAY 61 SAMPLES


#sweden
sweseqtab <- t(sweseqtab)
merged_swe <- merge(sweseqtab, swemeta, by="row.names")
metadata_swe <- merged_swe[, (ncol(sweseqtab)+2):ncol(merged_swe)]
rownames(metadata_swe) = merged_swe$X

sweseq <- merged_swe[, 2:(ncol(sweseqtab)+1)]
rownames(sweseq) = merged_swe$X

swe_otu <- otu_table(t(sweseq), taxa_are_rows = TRUE)
dim(swe_otu)
swe_tax <- tax_table(swetaxa)

swe_physeq_full <- phyloseq(otu_table(t(sweseqtab), taxa_are_rows = TRUE), swe_tax)
swe_physeq_full

swe_physeq <- phyloseq(swe_otu, swe_tax, sample_data(data.frame(metadata_swe)))

swe_physeq_clean <- subset_taxa(swe_physeq, kingdom == "Fungi")

swe_physeq_clean = prune_samples(sample_sums(swe_physeq_clean)>=1000, swe_physeq_clean)

saveRDS(swe_physeq_clean, file = file.path(figsPath, "swe_physeq_clean.rds"))

# SWEDEN 83 SAMPLES

