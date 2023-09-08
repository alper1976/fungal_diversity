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

# load(file=file.path(figsPath, "stats_b_diversity.RData"))

# image
# save.image(file=file.path(figsPath, "stats_b_diversity.RData"))

##load packages
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

if (!require("devtools")) {
   install.packages("devtools", dependencies = TRUE)
   library(devtools)
   }
if (!require("gridExtra")){
  install.packages("gridExtra")
  library(gridExtra)
  }
if (!require("dplyr")){
  install.packages("dplyr")
  library(dplyr)
  }
if (!require("EcoSimR")) {
   install.packages("EcoSimR", dependencies = TRUE)
   library(EcoSimR)
   }
if (!require("Biostrings")) {
   install.packages("Biostrings", dependencies = TRUE)
   library(Biostrings)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }
if (!require("R.utils")) {
   install.packages("R.utils", dependencies = TRUE)
   library(R.utils)
   }
if (!require("stringr")) {
   install.packages("stringr", dependencies = TRUE)
   library(stringr)
   }
if (!require("corrplot")) {
   install.packages("corrplot", dependencies = TRUE)
   library(corrplot)
   }
if (!require("tidyverse")) {
   install.packages("tidyverse", dependencies = TRUE)
   library(tidyverse)
   }
if (!require("GGally")) {
   install.packages("GGally", dependencies = TRUE)
   library(GGally)
   }
if (!require("cluster")) {
   install.packages("cluster", dependencies = TRUE)
   library(cluster)
   }
if (!require("factoextra")) {
   install.packages("factoextra", dependencies = TRUE)
   library(factoextra)
   }
if (!require("Hmisc")) {
   install.packages("Hmisc", dependencies = TRUE)
   library(Hmisc)
   }
if (!require("phylin")) {
   install.packages("phylin", dependencies = TRUE)
   library(phylin)
   }


################ Load local packages #####################################

sourceDirectory('/cluster/projects/nn9745k/scripts/amplicon_analysis_package')


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
swe_physeq_clean <- readRDS(file = "file.path(figsPath, swe_physeq_clean.rds)")
nor_physeq_clean <- readRDS(file = "file.path(figsPath, nor_physeq_clean.rds)")

################################ rarefied data ################################

#NORWAY
rarefied_nor = rarefy_even_depth(nor_physeq_clean, sample.size = 1000,
  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

nor_otu = otu_table(rarefied_nor)
nor_otu = t(nor_otu[rowSums(nor_otu) != 0,])

# SWEDEN
rarefied_swe = rarefy_even_depth(swe_physeq_clean, sample.size = 1000,
  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

swe_otu = otu_table(rarefied_swe)
swe_otu = t(swe_otu[rowSums(swe_otu) != 0,])


#### Ecosim

nor_cscore <- cooc_null_model(nor_otu, algo = "sim9", metric = "c_score", nReps = 1000, burn_in = 1000)

summary(nor_cscore)

cairo_ps(file.path(figsPath, "nor_cscore_burn_in.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
plot(nor_cscore,type="burn_in")
dev.off()
cairo_ps(file.path(figsPath, "nor_cscore_hist.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
plot(nor_cscore,type="hist")
dev.off()
cairo_ps(file.path(figsPath, "nor_cscore_cooc.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
plot(nor_cscore,type="cooc")
dev.off()

nor_checker <- cooc_null_model(nor_otu, algo = "sim3", metric = "checker", nReps = 1000, burn_in = 1000)

summary(nor_checker)

cairo_ps(file.path(figsPath, "nor_checker_burn_in.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(nor_checker,type="burn_in")
dev.off()
cairo_ps(file.path(figsPath, "nor_checker_hist.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(nor_checker,type="hist")
dev.off()
cairo_ps(file.path(figsPath, "nor_checker_cooc.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(nor_checker,type="cooc")
dev.off()

nor_combo <- cooc_null_model(nor_otu, algo = "sim9", metric = "species_combo", nReps = 1000, burn_in = 10000)
summary(nor_combo)

nor_combo_2 <- cooc_null_model(nor_otu, algo = "sim2", metric = "species_combo", nReps = 1000, burn_in = 10000)
summary(nor_combo_2)


cairo_ps(file.path(figsPath, "nor_combo_burn_in.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(nor_combo,type="burn_in")
dev.off()
cairo_ps(file.path(figsPath, "nor_combo_hist.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(nor_combo,type="hist")
dev.off()
cairo_ps(file.path(figsPath, "nor_combo_cooc.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(nor_combo,type="cooc")
dev.off()

swe_cscore <- cooc_null_model(swe_otu, algo = "sim9", metric = "c_score", nReps = 1000, burn_in = 1000)
summary(swe_cscore)

swe_cscore_2 <- cooc_null_model(swe_otu, algo = "sim2", metric = "c_score", nReps = 1000, burn_in = 1000)
summary(swe_cscore_2)


cairo_ps(file.path(figsPath, "swe_cscore_burn_in.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
plot(swe_cscore,type="burn_in")
dev.off()
cairo_ps(file.path(figsPath, "swe_cscore_hist.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
plot(swe_cscore,type="hist")
dev.off()
cairo_ps(file.path(figsPath, "swe_cscore_cooc.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
plot(swe_cscore,type="cooc")
dev.off()

swe_checker <- cooc_null_model(swe_otu, algo = "sim3", metric = "checker", nReps = 1000, burn_in = 1000)
summary(swe_checker)

swe_checker_2 <- cooc_null_model(swe_otu, algo = "sim2", metric = "checker", nReps = 1000, burn_in = 1000)
summary(swe_checker_2)
swe_checker_9 <- cooc_null_model(swe_otu, algo = "sim9", metric = "checker", nReps = 1000, burn_in = 1000)
summary(swe_checker_9)

cairo_ps(file.path(figsPath, "swe_checker_burn_in.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(swe_checker,type="burn_in")
dev.off()
cairo_ps(file.path(figsPath, "swe_checker_hist.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(swe_checker,type="hist")
dev.off()
cairo_ps(file.path(figsPath, "swe_checker_cooc.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(swe_checker,type="cooc")
dev.off()

swe_combo <- cooc_null_model(swe_otu, algo = "sim9", metric = "species_combo", nReps = 1000, burn_in = 10000)

summary(swe_combo)

cairo_ps(file.path(figsPath, "swe_combo_burn_in.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(swe_combo,type="burn_in")
dev.off()
cairo_ps(file.path(figsPath, "swe_combo_hist.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(swe_combo,type="hist")
dev.off()
cairo_ps(file.path(figsPath, "swe_combo_cooc.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
   plot(swe_combo,type="cooc")
dev.off()

swe_combo_2 <- cooc_null_model(swe_otu, algo = "sim2", metric = "species_combo", nReps = 1000, burn_in = 10000)

summary(swe_combo_2)


####### Raup Crick ######


nor_raup_crick = raup_crick(nor_otu, plot_names_in_col1 = TRUE, classic_metric = FALSE, split_ties = TRUE, reps = 9999, set_all_species_equal = FALSE, as.distance.matrix = TRUE, report_similarity = FALSE)

swe_raup_crick = raup_crick(swe_otu, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE)

# # envfit
#
# nor_mds = metaMDS(nor_raup_crick, distance = FALSE, autotransform = FALSE, try = 200)
# nor_mdata <- data.frame(sample_data(rarefied_nor)[,-c(1,4:6,14)])
#
# envfit_results <- envfit(nor_mds, nor_mdata, na.rm = TRUE, permu = 999)
# envfit_table_nor <- data.frame(round((envfit_results$vectors)$arrows, 3), round((envfit_results$vectors)$r, 3), round((envfit_results$vectors)$pvals# , 3))
# colnames(envfit_table_nor) <- c("DIM 1", "DIM 2", "R", "p")
# write.csv(envfit_table_nor, file.path(figsPath, "raupcrickdimivReg_nor.csv"))
#
#
# swe_mds = metaMDS(swe_raup_crick, distance = FALSE, autotransform = FALSE, try = 200)
# swe_mdata <- data.frame(sample_data(rarefied_swe)[,-c(1)])
#
# envfit_results <- envfit(swe_mds, swe_mdata, na.rm = TRUE, permu = 999)
# envfit_table_swe <- data.frame(round((envfit_results$vectors)$arrows, 3), round((envfit_results$vectors)$r, 3), round((envfit_results$vectors)$pvals# , 3))
# colnames(envfit_table_swe) <- c("DIM 1", "DIM 2", "R", "p")
# write.csv(envfit_table_swe, file.path(figsPath, "raupcrickdimivReg_swe.csv"))

## cluster based on environmental properties

# visualize data
metadata_nor = data.frame(sample_data(rarefied_nor))[,-c(1,4:6,15,16,22:26)]

cairo_ps(file.path(figsPath, "metadata_nor.eps"), width=160/25.4, height=160/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
metadata_nor %>%
  gather(attributes, value, 1:13) %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill = 'lightblue2', color = 'black') +
  facet_wrap(~attributes, scales = 'free_x') +
  labs(x="Values", y="Frequency") +
  theme_bw()
dev.off()

metadata_swe = data.frame(sample_data(rarefied_swe))[,-c(1,38,39)]

cairo_ps(file.path(figsPath, "metadata_swe.eps"), width=160/25.4, height=160/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
metadata_swe %>%
  gather(attributes, value, 1:13) %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill = 'lightblue2', color = 'black') +
  facet_wrap(~attributes, scales = 'free_x') +
  labs(x="Values", y="Frequency") +
  theme_bw()
dev.off()

# build correlation plot
metadata_nor_spearman <- rcorr(as.matrix(metadata_nor), type = "spearman")
metadata_nor_spearman$P[is.na(metadata_nor_spearman$P)]<-1

cairo_ps(file.path(figsPath, "metadata_nor_spearman.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  corrplot::corrplot(metadata_nor_spearman$r, type="upper", p.mat = metadata_nor_spearman$P, sig.level = 0.05, insig="blank", order="hclust", addrect=2)
dev.off()

metadata_swe_spearman <- rcorr(as.matrix(metadata_swe), type = "spearman")
metadata_swe_spearman$P[is.na(metadata_swe_spearman$P)]<-1

cairo_ps(file.path(figsPath, "metadata_swe_spearman.eps"), width=80/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  corrplot::corrplot(metadata_swe_spearman$r, type="upper", p.mat = metadata_swe_spearman$P, sig.level = 0.05, insig="blank", order="hclust", addrect=2)
dev.off()
##### classify lakes into categories #####
set.seed(123)

#Norway
metadata_nor = data.frame(sample_data(rarefied_nor))[,c("Marsh", "TP", "pH", "Annual.temp", "TOC", "Forest")]

metadata_nor_K2 <- kmeans(metadata_nor, centers = 2, nstart = 25)
print(metadata_nor_K2)

cairo_ps(file.path(figsPath, "metadata_nor_clusters2.eps"), width=120/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
   fviz_cluster(metadata_nor_K2, data = metadata_nor)
dev.off()

metadata_nor_K2$centers
metadata_nor_K2$pointsize
metadata_nor_K2$betweenss
metadata_nor_K2$withinss
metadata_nor_K2$tot.withinss
metadata_nor_K2$totss

metadata_nor_K3 <- kmeans(metadata_nor, centers = 3, nstart = 25)
metadata_nor_K4 <- kmeans(metadata_nor, centers = 4, nstart = 25)
metadata_nor_K5 <- kmeans(metadata_nor, centers = 5, nstart = 25)

p1 <- fviz_cluster(metadata_nor_K2, geom = "point", data = metadata_nor) + ggtitle(" K = 2") + theme_bw()
p2 <- fviz_cluster(metadata_nor_K3, geom = "point", data = metadata_nor) + ggtitle(" K = 3") + theme_bw()
p3 <- fviz_cluster(metadata_nor_K4, geom = "point", data = metadata_nor) + ggtitle(" K = 4") + theme_bw()
p4 <- fviz_cluster(metadata_nor_K5, geom = "point", data = metadata_nor) + ggtitle(" K = 5") + theme_bw()


cairo_ps(file.path(figsPath, "metadata_nor_clusters.eps"), width=160/25.4, height=160/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
   grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

# optimal number of clusters
p5 <- fviz_nbclust(x = metadata_nor,FUNcluster = kmeans, method = 'wss' ) + theme_bw()
p6 <- fviz_nbclust(x = metadata_nor,FUNcluster = kmeans, method = 'silhouette' ) + theme_bw()

# compute gap stats
gap_stat_nor <- clusGap(x = metadata_nor, FUN = kmeans, K.max = 15, nstart = 25, B = 50 )

# Print the result
print(gap_stat_nor, method = "firstmax")

p7 <- fviz_gap_stat(gap_stat_nor) + theme_bw()

cairo_ps(file.path(figsPath, "metadata_nor_optclusters.eps"), width=160/25.4, height=160/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  grid.arrange(p2, p5, p6, p7, nrow = 2)
dev.off()

# summarize clusters
metadata_nor %>%
  mutate(Cluster = metadata_nor_K3$cluster) %>%
  group_by(Cluster) %>%
  summarize_all('median')

#Sweden
metadata_swe = data.frame(sample_data(rarefied_swe))[,c("Wetland", "TotP", "pH", "Vattentemperatur", "TOC", "Forest")]

metadata_swe_K2 <- kmeans(metadata_swe, centers = 2, nstart = 25)
print(metadata_swe_K2)
cairo_ps(file.path(figsPath, "metadata_swe_clusters2.eps"), width=120/25.4, height=80/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
   fviz_cluster(metadata_swe_K2, data = metadata_swe)
dev.off()

metadata_swe_K2$centers
metadata_swe_K2$pointsize
metadata_swe_K2$betweenss
metadata_swe_K2$withinss
metadata_swe_K2$tot.withinss
metadata_swe_K2$totss

metadata_swe_K3 <- kmeans(metadata_swe, centers = 3, nstart = 25)
metadata_swe_K4 <- kmeans(metadata_swe, centers = 4, nstart = 25)
metadata_swe_K5 <- kmeans(metadata_swe, centers = 5, nstart = 25)

p1 <- fviz_cluster(metadata_swe_K2, geom = "point", data = metadata_swe) + ggtitle(" K = 2") + theme_bw()
p2 <- fviz_cluster(metadata_swe_K3, geom = "point", data = metadata_swe) + ggtitle(" K = 3") + theme_bw()
p3 <- fviz_cluster(metadata_swe_K4, geom = "point", data = metadata_swe) + ggtitle(" K = 4") + theme_bw()
p4 <- fviz_cluster(metadata_swe_K5, geom = "point", data = metadata_swe) + ggtitle(" K = 5") + theme_bw()


cairo_ps(file.path(figsPath, "metadata_swe_clusters.eps"), width=160/25.4, height=160/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
   grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

# optimal number of clusters
p5 <- fviz_nbclust(x = metadata_swe,FUNcluster = kmeans, method = 'wss' ) + theme_bw()
p6 <- fviz_nbclust(x = metadata_swe,FUNcluster = kmeans, method = 'silhouette' ) + theme_bw()

# compute gap stats
gap_stat_swe <- clusGap(x = metadata_swe, FUN = kmeans, K.max = 15, nstart = 25, B = 50 )

# Print the result
print(gap_stat_swe, method = "firstmax")

p7 <- fviz_gap_stat(gap_stat_swe) + theme_bw()

cairo_ps(file.path(figsPath, "metadata_swe_optclusters.eps"), width=160/25.4, height=160/25.4,
  pointsize = 4, bg = FALSE, fallback_resolution = 300)
  grid.arrange(p1, p5, p6, p7, nrow = 2)
dev.off()

# summarize clusters
metadata_swe %>%
  mutate(Cluster = metadata_swe_K2$cluster) %>%
  group_by(Cluster) %>%
  summarize_all('median')

## Analyse separate clusters based on Raup Crick t-test statistics

# first check if groups separate in beta diversity
# run permanova on clusters using bray curtis distance (from beta_diversity.R)

load(file=file.path(figsPath, "stats_b_diversity.RData"))

adonis2(swe_dist ~ metadata_swe_K2$cluster)
adonis2(nor_dist ~ metadata_nor_K3$cluster)

# adonis2(formula = swe_dist ~ metadata_swe_K2$cluster)
#                         Df SumOfSqs      R2      F Pr(>F)
# metadata_swe_K2$cluster  1    0.451 0.01369 1.1242  0.263
# Residual                81   32.486 0.98631
# Total                   82   32.937 1.00000

# adonis2(formula = nor_dist ~ metadata_nor_K3$cluster)
#                         Df SumOfSqs   R2      F Pr(>F)
# metadata_nor_K3$cluster  1   0.5271 0.02 1.2042  0.139
# Residual                59  25.8266 0.98
# Total                   60  26.3537 1.00

summary(nor_raup_crick[metadata_nor_K3$cluster==3])
summary(nor_raup_crick[metadata_nor_K3$cluster==2])
summary(nor_raup_crick[metadata_nor_K3$cluster==1])

raup_crick_clusters = rbind(as.matrix(nor_raup_crick[metadata_nor_K3$cluster==1]),
                            as.matrix(nor_raup_crick[metadata_nor_K3$cluster==2]),
                            as.matrix(nor_raup_crick[metadata_nor_K3$cluster==3]))

cluster_nor = rbind(as.matrix(rep("E", length(nor_raup_crick[metadata_nor_K3$cluster==1]))),
                    as.matrix(rep("O", length(nor_raup_crick[metadata_nor_K3$cluster==2]))),
                    as.matrix(rep("M", length(nor_raup_crick[metadata_nor_K3$cluster==3]))))
raup_crick_clusters_nor = data.frame(cbind(raup_crick_clusters, cluster_nor))
colnames(raup_crick_clusters_nor) = c("raup_crick", "cluster")

kruskal.test(raup_crick ~ cluster, data = raup_crick_clusters_nor)


summary(swe_raup_crick[metadata_swe_K2$cluster==2])
summary(swe_raup_crick[metadata_swe_K2$cluster==1])

raup_crick_clusters = rbind(as.matrix(swe_raup_crick[metadata_swe_K2$cluster==1]),as.matrix(swe_raup_crick[metadata_swe_K2$cluster==2]))
cluster_swe = rbind(as.matrix(rep("E", length(swe_raup_crick[metadata_swe_K2$cluster==1]))),
                    as.matrix(rep("O", length(swe_raup_crick[metadata_swe_K2$cluster==2]))))
raup_crick_clusters_swe = data.frame(cbind(raup_crick_clusters, cluster_swe))
colnames(raup_crick_clusters_swe) = c("raup_crick", "cluster")

kruskal.test(raup_crick ~ cluster, data = raup_crick_clusters_swe)



save.image(file=file.path(figsPath, "null_model.RData"))
