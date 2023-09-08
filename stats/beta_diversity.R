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
swe_physeq_clean <- readRDS(file = "file.path(figsPath, swe_physeq_clean.rds)")
nor_physeq_clean <- readRDS(file = "file.path(figsPath, nor_physeq_clean.rds)")

################################ beta diversity ############

## Rarefaction and coverage of diversity ##
# Number of ASVs
# Norway
nor_asv_number = apply(data.frame(otu_table(nor_physeq_clean)),2,function(x) sum(x > 0))


#NORWAY
rarefied_nor = rarefy_even_depth(nor_physeq_clean, sample.size = 1000,
  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)



cairo_ps(file.path(figsPath, "rarefaction_curve_nor.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  rarecurve(data.frame(t(otu_table(nor_physeq_clean))), step = 20, xlab = "Sample Size", ylab = "ASVs", drop = TRUE, label = FALSE, cex.lab = 2, cex.axis = 1.5)
dev.off()

cairo_ps(file.path(figsPath, "barplot_taxa_nor.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
  N <- 30
  barplot(sort(taxa_sums(nor_physeq_clean), TRUE)[1:N]/nsamples(nor_physeq_clean), las=2)
dev.off()


# region-wide SAC
spec_accum_random <- specaccum(t(otu_table(nor_physeq_clean)), method = "random", permutations = 100,
          conditioned =TRUE, gamma = "jack1")

mod_random <- fitspecaccum(spec_accum_random, "arrh")

cairo_ps(file.path(figsPath, "sac_nor.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  plot(mod_random, col="hotpink", xlab = "Number of sites", ylab = "ASVs", cex.lab = 2, cex.axis = 1.5)
  boxplot(spec_accum_random, col = "#D1BBD7", border = "#F7F056", lty=1, cex=0.5, add= TRUE)
dev.off()

nor_otu = otu_table(rarefied_nor)
nor_otu = t(nor_otu[rowSums(nor_otu) != 0,])

# Standardize data
nor_otu_stand = decostand(nor_otu, method="hellinger")

# Create NMDS object
nor_mds = metaMDS(nor_otu_stand, distance = "bray", autotransform = FALSE, try = 200)

## subset NAs

nor_mdata <- data.frame(sample_data(rarefied_nor)[,-c(1,4:6, 15)])

envfit_results <- envfit(nor_mds, nor_mdata, na.rm = TRUE, permu = 999)
envfit_table_nor <- data.frame(round((envfit_results$vectors)$arrows, 3), round((envfit_results$vectors)$r, 3), round((envfit_results$vectors)$pvals, 3))
colnames(envfit_table_nor) <- c("DIM 1", "DIM 2", "R", "p")
write.csv(envfit_table_nor, file.path(figsPath, "betaddimivReg_nor.csv"))


# SWEDEN
swe_asv_number = apply(data.frame(otu_table(swe_physeq_clean)),2,function(x) sum(x > 0))

rarefied_swe = rarefy_even_depth(swe_physeq_clean, sample.size = 1000,
  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

cairo_ps(file.path(figsPath, "rarefaction_curve_swe.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  rarecurve(data.frame(t(otu_table(swe_physeq_clean))), step = 20, xlab = "Sample Size", ylab = "ASVs", label = FALSE, cex.lab = 2, cex.axis = 1.5)
dev.off()

cairo_ps(file.path(figsPath, "barplot_taxa_swe.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
  N <- 30
  barplot(sort(taxa_sums(swe_physeq_clean), TRUE)[1:N]/nsamples(swe_physeq_clean), las=2)
dev.off()

# region-wide SAC
spec_accum_random <- specaccum(t(otu_table(swe_physeq_clean)), method = "random", permutations = 100,
          conditioned =TRUE, gamma = "jack1")

mod_random <- fitspecaccum(spec_accum_random, "arrh")

cairo_ps(file.path(figsPath, "sac_swe.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  plot(mod_random, col="hotpink", xlab = "Number of sites", ylab = "ASVs", cex.lab = 2, cex.axis = 1.5)
  boxplot(spec_accum_random, col = "#D1BBD7", border = "#F7F056", lty=1, cex=0.5, add= TRUE)
dev.off()

swe_otu = otu_table(rarefied_swe)
swe_otu = t(swe_otu[rowSums(swe_otu) != 0,])


# Standardize data
swe_otu_stand = decostand(swe_otu, method="hellinger")

# Create NMDS object
swe_mds = metaMDS(swe_otu_stand, distance = "bray", autotransform = FALSE, try = 200)

## subset NAs


### NEEDS FIXING
swe_mdata <- data.frame(sample_data(rarefied_swe)[,-c(1)])

envfit_results <- envfit(swe_mds, swe_mdata, na.rm = TRUE, permu = 999)
envfit_table_swe <- data.frame(round((envfit_results$vectors)$arrows, 3), round((envfit_results$vectors)$r, 3), round((envfit_results$vectors)$pvals, 3))
colnames(envfit_table_swe) <- c("DIM 1", "DIM 2", "R", "p")
write.csv(envfit_table_swe, file.path(figsPath, "betaddimivReg_swe.csv"))

tmp <- data.frame(R_value = c(envfit_table_swe$R, envfit_table_nor$R),
                  p_value = c(envfit_table_swe$p, envfit_table_nor$p),
                  Dataset = rep(c("Sweden", "Norway"),
                               each = length(envfit_table_swe$p)),
                  env_var = rownames(envfit_table_nor))

cairo_ps(file.path(figs_dir, "envfit_plot.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
ggplot(tmp, aes(x = Domain, y = env_var)) +
  geom_point(aes(size = 1/p_value, color = R_value), alpha = 0.7) +
  scale_color_gradient(low = "#E69F00", high = "#56B4E9") +
  scale_size(breaks = c(20,100,1000)) +
  theme_bw()
dev.off()



################################ Distance decay #################################################

#Functions for calculating distances between lakes
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
   # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
   # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.

   if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
   if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
   else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
   else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
   m[tri] <- t(m)[tri]
   return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
   # Returns a matrix (M) of distances between geographic points.
   # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
   # (df.geopoints$lat[j], df.geopoints$lon[j]).
   # The row and column names are given by df.geopoints$name.

   GeoDistanceInMetres <- function(g1, g2){
      # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
      # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
      # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
      # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
      # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
      DistM <- function(g1, g2){
         require("geosphere")
         return(ifelse(g1$index > g2$index, 0, distm (c(g1$lon, g1$lat), c(g2$lon, g2$lat), fun = distHaversine)))
      }
      return(mapply(DistM, g1, g2))
   }

   n.geopoints <- nrow(df.geopoints)

   # The index column is used to ensure we only do calculations for the upper triangle of points
   df.geopoints$index <- 1:n.geopoints

   # Create a list of lists
   list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})

   # Get a matrix of distances (in metres)
   mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")

   # Set the row and column names
   rownames(mat.distances) <- df.geopoints$name
   colnames(mat.distances) <- df.geopoints$name

   return(mat.distances)
}

# Norway distance
nor_stand <- decostand(nor_otu, method="hellinger", MARGIN = 1)
nor_dist <- vegdist(nor_stand, method = "bray")

df.lakes <- data.frame(name = rownames(nor_mdata), lat  = nor_mdata[, "Latitude"], lon  = nor_mdata[, "Longitude"])
lakes_dist <- round(GeoDistanceInMetresMatrix(df.lakes) / 1000)

lakes_dist[upper.tri(lakes_dist)] <- NA

mantel_results <- mantel.rtest(nor_dist, as.dist(lakes_dist), nrepet= 9999)
mantel_results

cairo_ps(file.path(figsPath, "divDist_mantel_nor.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  plot(as.dist(lakes_dist), nor_dist,
            xlab="lake distance [km]", ylab="community distance [Bray Curtis]")
  reg <- lm(nor_dist ~ as.dist(lakes_dist))
  abline(reg)
dev.off()

#partial mantel test
meta_nor <- data.matrix(sample_data(rarefied_nor)[,-c(1:7)])
meta_dist_nor <- vegdist(meta_nor, method = "euclidean", na.rm = TRUE)

mantel.partial(nor_dist, meta_dist_nor, as.dist(lakes_dist), permutations = 999, method = "pearson")

mantel.partial(nor_dist, as.dist(lakes_dist), meta_dist_nor, permutations = 999, method = "pearson")


# Sweden distance
swe_stand <- decostand(swe_otu, method="hellinger", MARGIN = 1)
swe_dist <- vegdist(swe_stand, method = "bray")

df.lakes <- data.frame(name = rownames(swe_mdata), lat  = swe_mdata[, "Latitude"], lon  = swe_mdata[, "Longitude"])
lakes_dist <- round(GeoDistanceInMetresMatrix(df.lakes) / 1000)

lakes_dist[upper.tri(lakes_dist)] <- NA

mantel_results <- mantel.rtest(swe_dist, as.dist(lakes_dist), nrepet= 9999)
mantel_results

cairo_ps(file.path(figsPath, "divDist_mantel_swe.eps"), width=80/25.4, height=80/25.4, pointsize = 4, bg = FALSE, fallback_resolution = 300)
  plot(as.dist(lakes_dist), swe_dist,
            xlab="lake distance [km]", ylab="community distance [Bray Curtis]")
  reg <- lm(swe_dist ~ as.dist(lakes_dist))
  abline(reg)
dev.off()

#partial mantel test
meta_swe <- data.matrix(sample_data(rarefied_swe)[,-c(1:7)])
meta_dist_swe <- vegdist(meta_swe, method = "euclidean", na.rm = TRUE)

mantel.partial(swe_dist, meta_dist_swe, as.dist(lakes_dist), permutations = 999, method = "pearson")

mantel.partial(swe_dist, as.dist(lakes_dist), meta_dist_swe, permutations = 999, method = "pearson")

save.image(file=file.path(figsPath, "stats_b_diversity.RData"))
