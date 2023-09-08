## saga
#module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1 # in addition all bioconductor packages will be available
#module load GDAL/3.5.0-foss-2022a # this is needed to make maps
#module load MariaDB-connector-c/3.1.7-GCCcore-9.3.0 # this is needed to make maps




.libPaths(c("/cluster/projects/nn9745k/Rpackages_4_2_1", .libPaths()))

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
#####################################

nor_physeq_clean <- readRDS(file = "file.path(figsPath, nor_physeq_clean.rds)")

# get overall reads per phylum
phylum_table_nor <- phyloseq_summarize_taxa(nor_physeq_clean, 'phylum')

sort(taxa_sums(phylum_table_nor)/sum(taxa_sums(phylum_table_nor)))

# get overall reads per genus
genus_table_nor <- phyloseq_summarize_taxa(nor_physeq_clean, 'genus')

sort(taxa_sums(genus_table_nor)/sum(taxa_sums(genus_table_nor)))

phylum_table_nor <- phyloseq_summarize_taxa(nor_physeq_clean, 'phylum')

sort(taxa_sums(phylum_table_nor)/sum(taxa_sums(phylum_table_nor)))

# get overall reads per genus
genus_table_nor <- phyloseq_summarize_taxa(nor_physeq_clean, 'genus')

sort(taxa_sums(genus_table_nor)/sum(taxa_sums(genus_table_nor)))

# get richness stats from each phylum

taxa = get_taxa_unique(nor_physeq_clean, "phylum")

ntaxa(nor_physeq_clean)

tdt = data.table(tax_table(nor_physeq_clean),
                 TotalCounts = taxa_sums(nor_physeq_clean),
                 OTU = taxa_names(nor_physeq_clean))

taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]
# Cumulation curve
plot_cumsum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) +
                  geom_point() +
                  theme_bw() +
                  xlab("Abundance per ASV") +
                  ylab("Cumulative number of ASVs") +
                  ggtitle("Cumulation curve")
# plot_cumsum

png(file.path(figsPath, "pCumSum_nor.png"), width=86, height=86, units="mm", pointsize = 4, bg = FALSE, res = 300)
  plot_cumsum + xlim(0, 100)
dev.off()

## widely vs narrowly distributed taxa (genera)

# Prevalence
mdt = fast_melt(nor_physeq_clean)
prevdt = mdt[, list(Prevalence = sum(count > 0),
                    TotalCounts = sum(count)),
                    by = TaxaID]

over10 = prevdt[prevdt$Prevalence >10,]
max_prev = prevdt[prevdt$Prevalence == max(prevdt$Prevalence),]
tax_table(prune_taxa(over10$TaxaID, nor_physeq_clean))

plot_prev = ggplot(prevdt, aes(Prevalence)) +
                geom_histogram(binwidth = 5) +
                theme_bw() +
                ylab("Cumulative number of ASVs") +
                ggtitle("Histogram of ASV Prevalence")

png(file.path(figsPath, "hist_prev_nor.png"), width=86, height=86, units="mm", pointsize = 4, bg = FALSE, res = 300)
  plot_prev
dev.off()

prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]
plot_prevcumsum = ggplot(prevcumsum, aes(Prevalence, CumSum)) +
                    geom_point() +
                    xlab("Prevalence") +
                    ylab("Cumulative number of ASVs") +
                    theme_bw() +
                    ggtitle("Prevalence cumulation curve")

png(file.path(figsPath, "plot_prev_nor.png"), width=86, height=86, units="mm", pointsize = 4, bg = FALSE, res = 300)
  plot_prevcumsum
dev.off()


plot_prev_totcounts = ggplot(prevdt, aes(Prevalence, TotalCounts)) +
                        geom_point(size = 2, alpha = 0.75) +
                        theme_bw() +
                        scale_y_log10()

png(file.path(figsPath, "plot_prev_totcounts_nor.png"), width=86, height=86, units="mm", pointsize = 2, bg = FALSE, res = 300)
  plot_prev_totcounts
dev.off()

dim(prevdt)

ordered_prevdt = prevdt[order(prevdt$Prevalence, decreasing = TRUE),]
head(ordered_prevdt)

# prevalence plots - classification into rare, common, transient, etc.
total_reads = sum(colSums(otu_table(nor_physeq_clean)))

mean_rel_counts = mdt[, list(RelPrevalence = sum(count > 0)/72),
                         by = TaxaID]
mean_rel_occ = mdt[, list(RelCounts = sum(count)/total_reads),
                         by = TaxaID]
mdt = cbind(mdt, mean_rel_counts$RelPrevalence, mean_rel_occ$RelCounts)



plot_relcounts_relprev = ggplot(mdt, aes(V3, V2)) +
                        geom_point(size = 2, alpha = 0.75) +
                        theme_bw() +
                        labs(x = "Relative Counts", y = "Relative Prevalence") +
                        scale_y_log10()

png(file.path(figsPath, "plot_relcounts_relprev_nor.png"), width=86, height=86, units="mm", pointsize = 2, bg = FALSE, res = 300)
  plot_relcounts_relprev
dev.off()

## Lakes
# 
asv_table <- as.data.frame(t(otu_table(nor_physeq_clean)))
asv_table <- asv_table[,colSums(asv_table) > 0]
# make binary table
asv_table_binary <- asv_table
asv_table_binary[asv_table_binary > 0] <- 1

# get those ASVs occurring in a single system
asv_table_binary1 <- asv_table_binary[,colSums(asv_table_binary) == 1]

max(rowSums(asv_table_binary1))
min(rowSums(asv_table_binary1))
mean(rowSums(asv_table_binary1))

asv_unique_lake_percentage <- rowSums(asv_table_binary1)/rowSums(asv_table_binary)

max(asv_unique_lake_percentage)
min(asv_unique_lake_percentage)
mean(asv_unique_lake_percentage)

# get reads from ASV occurring in single system

asv_table_unique <- asv_table[,colSums(asv_table_binary) == 1]
sums_asv_table_unique <- rowSums(asv_table_unique)
asv_unique_read_percentage <- rowSums(asv_table_unique)/rowSums(asv_table)
ncol(asv_table_unique)/ncol(asv_table_binary) # percent of lake unique ASVs
ncol(asv_table_unique)


max(asv_unique_read_percentage)
min(asv_unique_read_percentage)
mean(asv_unique_read_percentage)


# with ASVs < 3

asv_table_binary23 <- asv_table_binary[,colSums(asv_table_binary) > 1 & colSums(asv_table_binary) < 4]
asv_23_lake_percentage <- rowSums(asv_table_binary23)/rowSums(asv_table_binary)
mean(asv_23_lake_percentage)
ncol(asv_table_binary23)/ncol(asv_table_binary) # percent of ASVs occurring in 2 to 3 systems

asv_table_2_3 <- asv_table[,colSums(asv_table_binary) > 1 & colSums(asv_table_binary) < 4]
asv_table_2_3_percentage <- rowSums(asv_table_2_3)/rowSums(asv_table)

dim(asv_table_2_3)
# with ASVs > 3
asv_table_binary_more3 <- asv_table_binary[,colSums(asv_table_binary) > 3]
asv_more3_lake_percentage <- rowSums(asv_table_binary_more3)/rowSums(asv_table_binary)
mean(asv_more3_lake_percentage)
ncol(asv_table_binary_more3)

ncol(asv_table_binary_more3)/ncol(asv_table_binary) # percent of ASVs occurring in more than 3 systems


asv_table_more3 <- asv_table[,colSums(asv_table_binary) > 3]

asv_table_more3_percentage <- rowSums(asv_table_more3)/rowSums(asv_table)

max(asv_table_more3_percentage)
min(asv_table_more3_percentage)
mean(asv_table_more3_percentage)

lake_occupancy <- rbind(asv_unique_read_percentage, asv_table_2_3_percentage, asv_table_more3_percentage)

melted_lake_occupancy <- melt(lake_occupancy)
colnames(melted_lake_occupancy) <- c("occupancy_type", "Lake", "percentage")

lake_occupancy_plot <- ggplot() +
                        geom_bar(aes(y = percentage, x = Lake, fill = occupancy_type),
                          data = melted_lake_occupancy,
                          stat="identity") +
                        labs(x = "Lakes", y = "Relative abundance", fill = "Occupancy type") +
                        scale_fill_discrete(name = "Occupancy type", labels = c("unique", "2-3 lakes", "> 3 lakes")) +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
                          legend.text=element_text(size = 6),
                          legend.key.size = unit(0.05, "in"),
                          legend.position = "top",
                          legend.box = "vertical")




## plot Figure S3

cairo_ps(file.path(figsPath, "Figure_S3A.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  lake_occupancy_plot
dev.off()

cairo_ps(file.path(figsPath, "Figure_S3B.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  plot_prevcumsum
dev.off()

cairo_ps(file.path(figsPath, "Figure_S3C.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  plot_relcounts_relprev
dev.off()



## get number of sequences not assigned to genera and family
prevdt_taxa <- mdt[, list(Prevalence = sum(count > 0),
                    TotalCounts = sum(count)),
                    by = c("kingdom", "phylum", "class", "order", "family", "genus", "species")]

abundorder_prevdt_taxa <- prevdt_taxa[order(prevdt_taxa$TotalCounts, decreasing = TRUE),]
head(abundorder_prevdt_taxa)
prevorder_prevdt_taxa <- prevdt_taxa[order(prevdt_taxa$Prevalence, decreasing = TRUE),]
head(prevorder_prevdt_taxa)

prevdt_genus <- mdt[, list(
                    TotalCounts = sum(count)),
                    by = genus]
abundorder_prevdt_genus <- prevdt_genus[order(prevdt_genus$TotalCounts, decreasing = TRUE),]
head(abundorder_prevdt_genus)

prevdt_family <- mdt[, list(
                    TotalCounts = sum(count)),
                    by = family]
abundorder_prevdt_family <- prevdt_family[order(prevdt_family$TotalCounts, decreasing = TRUE),]
head(abundorder_prevdt_family)


###### ADDITIONAL ANALYSES ######


###### FINAL STEPS ######
save.image(file=file.path(figsPath, "stats.RData"))
