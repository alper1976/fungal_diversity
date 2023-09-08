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

# load(file=file.path(figsPath, "stats_a_diversity.RData"))

# image
# save.image(file=file.path(figsPath, "stats_a_diversity.RData"))

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
if (!require("RcmdrMisc")) {
   install.packages("RcmdrMisc", dependencies = TRUE)
   library(RcmdrMisc)
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
if (!require("plsdepot")) {
   install.packages("plsdepot", dependencies = TRUE)
   library(plsdepot)
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
swe_physeq_clean <- readRDS(file = "file.path(figsPath, swe_physeq_clean.rds)")
nor_physeq_clean <- readRDS(file = "file.path(figsPath, nor_physeq_clean.rds)")

rarefied_nor = rarefy_even_depth(nor_physeq_clean, sample.size = 1000,
  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

rarefied_swe = rarefy_even_depth(swe_physeq_clean, sample.size = 1000,
  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

################################ Alpha diversity #################################################

## richness

richness_nor <- estimate_richness(rarefied_nor, measures=c("Observed", "InvSimpson", "Shannon", "ACE"))
summary(richness_nor)


richness_swe <- estimate_richness(rarefied_swe, measures=c("Observed", "InvSimpson", "Shannon", "ACE"))
summary(richness_swe)


# Pielou's Evenness
richness_nor$Pielou <- richness_nor$Shannon/log(richness_nor$Observed)
richness_swe$Pielou <- richness_swe$Shannon/log(richness_swe$Observed)


#adding metadata
metadata_nor = phyloseq::sample_data(nor_physeq_clean)
richness_nor = cbind(richness_nor, metadata_nor)
richness_nor = richness_nor[,-c(3:5,7,21)]


metadata_swe = phyloseq::sample_data(swe_physeq_clean)
richness_swe = cbind(richness_swe, metadata_swe)


###### statistics ######

corr_spearman_nor <- rcorr(as.matrix(richness_nor), type = "spearman")
corr_spearman_swe <- rcorr(as.matrix(richness_swe), type = "spearman")

# get Holm's adjusted p-value

spearman_adj_nor = print_rcorr_adjust(rcorr.adjust(richness_nor, type = "spearman",
                                              use = "pairwise.complete.obs"))

df_1 = dplyr::mutate_all(as.data.frame(spearman_adj_nor, stringsAsFactors = FALSE),
        function(x) as.numeric(as.character(x)))
df_1[is.na(df_1)] = 0.0001
df_2 = as.matrix(df_1)

corr_spearman_nor$P[is.na(corr_spearman_nor$P)] = 0.0001


postscript(file.path(figsPath, "correlation_matrix_diversity_adj_nor.eps"), width=80/25.4, height=80/25.4, pointsize = 4)
  corrplot::corrplot(corr_spearman_nor$r, type="upper", p.mat = df_2, sig.level = 0.05, insig="blank", order="hclust", addrect=2)
dev.off()

postscript(file.path(figsPath, "correlation_matrix_diversity_nor.eps"), width=80/25.4, height=80/25.4, pointsize = 4)
  corrplot::corrplot(corr_spearman_nor$r, type="upper", p.mat = corr_spearman_nor$P, sig.level = 0.05, insig="blank", order="hclust", addrect=2)
dev.off()

spearman_adj_swe = print_rcorr_adjust(rcorr.adjust(richness_swe, type = "spearman",
                                              use = "pairwise.complete.obs"))

df_1 = dplyr::mutate_all(as.data.frame(spearman_adj_swe, stringsAsFactors = FALSE),
        function(x) as.numeric(as.character(x)))
df_1[is.na(df_1)] = 0.0001
df_2 = as.matrix(df_1)

corr_spearman_swe$P[is.na(corr_spearman_swe$P)] = 0.0001


postscript(file.path(figsPath, "correlation_matrix_diversity_adj_swe.eps"), width=80/25.4, height=80/25.4, pointsize = 4)
  corrplot::corrplot(corr_spearman_swe$r, type="upper", p.mat = df_2, sig.level = 0.05, insig="blank", order="hclust", addrect=2)
dev.off()

postscript(file.path(figsPath, "correlation_matrix_diversity_swe.eps"), width=80/25.4, height=80/25.4, pointsize = 4)
  corrplot::corrplot(corr_spearman_swe$r, type="upper", p.mat = corr_spearman_swe$P, sig.level = 0.05, insig="blank", order="hclust", addrect=2)
dev.off()

##### PLS #####
# PLS/PCA
# NORWAY
independent_var_nor = scale(metadata_nor[,c(2,3, 7:9,8:13, 16:45)])
dependent_var_nor = richness_nor[,c("ACE", "Pielou")]
dependent_var_nor[which(dependent_var_nor$ACE == "NaN"),] <- richness_nor$Observed[which(dependent_var_nor$ACE == "NaN")]

pls2_reg_nor = plsreg2(independent_var_nor, dependent_var_nor, comps=5, crosval=TRUE)

summary(pls2_reg_nor)

pls2_reg_nor$expvar
pls2_reg_nor$Q2cum
pls2_reg_nor$VIP
cairo_ps(file.path(figsPath, "PLS_diversity_nor.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)

  plot(pls2_reg_nor, xlab = paste0("PC1, R2X = ",
                                   round(pls2_reg_euk$expvar[1], 2),
                                   ", R2Y = ",
                                   round(pls2_reg_euk$expvar[11],2)),
                     ylab = paste0("PC2, R2X = ",
                                   round(pls2_reg_euk$expvar[2], 2),
                                   ", R2Y = ",
                                   round(pls2_reg_euk$expvar[12],2)),
                 cex = 0.7,
                 cex.axis = 2,
                 cex.lab = 2,
                 main = NA)
dev.off()

# Sweden

independent_var_swe = scale(metadata_swe[,c(2:39)])
dependent_var_swe = richness_swe[,c("ACE", "Pielou")]
dependent_var_swe[which(dependent_var_swe$ACE == "NaN"),] <- richness_swe$Observed[which(dependent_var_swe$ACE == "NaN")]

pls2_reg_swe = plsreg2(independent_var_swe[5:nrow(independent_var_swe),], dependent_var_swe[5:nrow(dependent_var_swe),], comps=5, crosval=TRUE)

summary(pls2_reg_swe)

pls2_reg_swe$expvar
pls2_reg_swe$Q2cum
pls2_reg_swe$VIP

cairo_ps(file.path(figs_dir, "PLS_diversity_swe.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)

  plot(pls2_reg_swe, xlab = paste0("PC1, R2X = ",
                                   round(pls2_reg_bac$expvar[1], 2),
                                   ", R2Y = ",
                                   round(pls2_reg_bac$expvar[11],2)),
                     ylab = paste0("PC2, R2X = ",
                                   round(pls2_reg_bac$expvar[2], 2),
                                   ", R2Y = ",
                                   round(pls2_reg_bac$expvar[12],2)),
                 cex = 0.7,
                 cex.axis = 2,
                 cex.lab = 2,
                 main = NA)
dev.off()

library(pls)
library(tidyverse)

y_bac <- as.matrix(diversity_metadata[c(1:5,10:29),c(4)])
x_bac <- as.matrix(sample_metadata[c(1:5,10:29),c(1:3, 5, 7:9, 11,13, 15:18,21:23)])
df_plsr <- mvr(y_bac  ~ x_bac , ncomp = 2, method = "oscorespls" , scale = T)

df_coef <- as.data.frame(coef(df_plsr, ncomp =1:2))
df_coef <- df_coef %>%
  dplyr::mutate(variables = rownames(.)) %>%
  dplyr::mutate(variables = factor(variables,
                                   levels = c(
                                      "altitude", "gl_dist",  "Ccounts",  "Cond",
                                      "TOC", "DOC", "TN", "TP", "O2", "CO2", "CH4",
                                      "N2O", "bird", "CO2_sat", "CH4_sat", "N2O_sat")
                                      ))

colnames(df_coef)[1] <- "regression_coefficients"
df_coef  <- df_coef[, -c(2)]

# VIP

VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")

  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


vip_bac <- as.data.frame(VIP(df_plsr))
vip_bac_comp1 <- vip_bac[1, ]
row.names(vip_bac_comp1)<-NULL

combine_vip_coef <- vip_bac_comp1 %>%
  tidyr::gather("variables", "VIP") %>%
  dplyr::full_join(., df_coef, by = c("variables")) %>%
  dplyr::mutate(variables = factor(variables,
                                   levels = c(
                                      "altitude", "gl_dist",  "Ccounts",  "Cond",
                                      "TOC", "DOC", "TN", "TP", "O2", "CO2", "CH4",
                                      "N2O", "bird", "CO2_sat", "CH4_sat", "N2O_sat")
                                      ))


p1_vip <-
  ggplot(combine_vip_coef, aes(x = variables,  y = VIP, group =1))+
  geom_bar(stat="identity",  fill = "black") +
  geom_hline(yintercept = 1, size = 0.55, linetype = 3) +
  theme_bw()+
  theme(axis.text.x = element_text(angle=65,
                                   hjust=1,
                                   size = 6),
        axis.title.y = element_text(size = 10))+
  labs(x= "")

p2_coef <-
  ggplot(df_coef, aes(x = variables, y = regression_coefficients, group = 1))+  geom_bar(stat = "identity",  fill = "black")+
  theme(axis.text.x = element_text(angle=65,
                                   hjust=1,
                                   size = 8),
        axis.title.y = element_text(size = 2))+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  labs(x="")

cairo_ps(file.path(figs_dir, "PLS_diversity_bac_coeff.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)
p2_coef
dev.off()

cairo_ps(file.path(figs_dir, "PLS_diversity_bac_vip.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)
p1_vip
dev.off()

y_euk <- as.matrix(diversity_metadata[c(1:5,10:29),c(13)])
x_euk <- as.matrix(sample_metadata[c(1:5,10:29),c(1:3, 5, 7:9, 11,13, 15:18,21:23)])
df_plsr <- mvr(y_euk  ~ x_euk , ncomp = 2, method = "oscorespls" , scale = T)

df_coef <- as.data.frame(coef(df_plsr, ncomp =1:2))
df_coef <- df_coef %>%
  dplyr::mutate(variables = rownames(.)) %>%
  dplyr::mutate(variables = factor(variables,
                                   levels = c(
                                      "altitude", "gl_dist",  "Ccounts",  "Cond",
                                      "TOC", "DOC", "TN", "TP", "O2", "CO2", "CH4",
                                      "N2O", "bird", "CO2_sat", "CH4_sat", "N2O_sat")
                                      ))

colnames(df_coef)[1] <- "regression_coefficients"
df_coef  <- df_coef[, -c(2)]

# VIP

VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")

  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


vip_euk <- as.data.frame(VIP(df_plsr))
vip_euk_comp1 <- vip_euk[1, ]
row.names(vip_euk_comp1)<-NULL

combine_vip_coef <- vip_euk_comp1 %>%
  tidyr::gather("variables", "VIP") %>%
  dplyr::full_join(., df_coef, by = c("variables")) %>%
  dplyr::mutate(variables = factor(variables,
                                   levels = c(
                                      "altitude", "gl_dist",  "Ccounts",  "Cond",
                                      "TOC", "DOC", "TN", "TP", "O2", "CO2", "CH4",
                                      "N2O", "bird", "CO2_sat", "CH4_sat", "N2O_sat")
                                      ))


p1_vip <-
  ggplot(combine_vip_coef, aes(x = variables,  y = VIP, group =1))+
  geom_bar(stat="identity",  fill = "black") +
  geom_hline(yintercept = 1, size = 0.55, linetype = 3) +
  theme_bw()+
  theme(axis.text.x = element_text(angle=65,
                                   hjust=1,
                                   size = 6),
        axis.title.y = element_text(size = 10))+
  labs(x= "")

p2_coef <-
  ggplot(df_coef, aes(x = variables, y = regression_coefficients, group = 1))+  geom_bar(stat = "identity",  fill = "black")+
  theme(axis.text.x = element_text(angle=65,
                                   hjust=1,
                                   size = 8),
        axis.title.y = element_text(size = 2))+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  labs(x="")

cairo_ps(file.path(figs_dir, "PLS_diversity_euk_coeff.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)
p2_coef
dev.off()

cairo_ps(file.path(figs_dir, "PLS_diversity_euk_vip.eps"), width=80/25.4, height=80/25.4,
                        pointsize = 4, bg = FALSE, fallback_resolution = 300)
p1_vip
dev.off()



# NO or WEAK RELATONSHIPS between environmental properties and fungal alpha and beta diversity


###### FINAL STEPS ######
save.image(file=file.path(figsPath, "stats_a_diversity.RData"))
