########################## Explore metadata ############################
## Needs to be run in a separate R session

# module load R/4.0.0-fosscuda-2020a
.libPaths(c("/cluster/projects/nn9745k/Rpackages_4_0_0", .libPaths()))
norway <- file.path("/cluster/projects/nn9744k/02_results/29_euk_20220325/Pedro/final_tax_classification") # <---------- Change to your path

sweden <- file.path("/cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/final_tax_classification") # <---------- Change to your path

figsPath <- file.path("/cluster/projects/nn9744k/02_results/30_scandinavian_fungi/figs")



if (!require("xlsx")) {
   install.packages("xlsx", dependencies = TRUE)
   library(xlsx)
   }
if (!require("stringr")) {
   install.packages("stringr", dependencies = TRUE)
   library(stringr)
   }
if (!require("bestNormalize")) {
   install.packages("bestNormalize", dependencies = TRUE)
   library(bestNormalize)
   }
if (!require("mice")) {
   install.packages("mice", dependencies = TRUE)
   library(mice)
   }
library(dplyr)

# local functions
metadata_transform <- function(NotStand){
  stand_metadata = NotStand;
  for (i in 1:dim(NotStand)[2]) {
    if (is.numeric(NotStand[,i])) {
      output = bestNormalize(NotStand[,i], r = 1, k = 5);
      stand_metadata [,i] = output$x.t
    }
    i = i + 1
  }
  return(stand_metadata)
}

######################################################
### Data loading

normeta = data.frame(read.xlsx(file.path(norway, "100lakes_norway_meta_sample.xlsx"), 1, header=T, row.names = "sample", check.names=FALSE))

swemeta = data.frame(read.table(file.path("/cluster/projects/nn9745k/02_results/02_hundred_lakes/lake_metadata", "100lakes_metadata_final.csv"), sep = ";", row.names = "fungi_seq_id", header = TRUE))
swemeta = swemeta[1:103,c(2,10,11,26:73)]
swemeta_names = str_extract(rownames(swemeta), "([^_]+)")
rownames(swemeta) = swemeta_names


# norway
dim(normeta)
colnames(normeta)

# parse metadata

## remove redundant variables
variables_to_remove = c("Campaign", "ID", "Lake.name", "Date")

normeta_cleaned <- normeta[,!colnames(normeta) %in% variables_to_remove]

normeta_trans = metadata_transform(normeta_cleaned)

normeta_trans_interpolated = complete(mice(normeta_trans))

variables_to_remove = c("Latitude",
                        "Longitude")


normeta_trans_interpolated_clean <- normeta_trans_interpolated[,!colnames(normeta_trans_interpolated) %in% variables_to_remove]
normeta_trans_interpolated_clean = cbind(normeta_trans_interpolated_clean,
                                         normeta_cleaned[,c("Latitude", "Longitude")])

write.csv(normeta_trans_interpolated_clean, file.path(norway, "normeta_trans_interpolated.csv"))




# sweden

dim(swemeta)
colnames(swemeta)

# parse metadata

## remove redundant variables
variables_to_remove = c("logCLratio", "Urbankm", "Agriculturalkm", "SO4Cl",
                        "Kond_25", "Grassbarekm", "Wetlandkm", "NO2.NO3.N",
                        "AbsF254", "AbsF436")

swemeta_cleaned <- swemeta[,!colnames(swemeta) %in% variables_to_remove]

swemeta_trans = metadata_transform(swemeta_cleaned)

swemeta_trans_cleaned = swemeta_trans[1:91,]

swemeta_trans_interpolated = complete(mice(swemeta_trans_cleaned))

variables_to_remove = c("lakeID",
                        "Stationskoordinat.N.X",
                        "Stationskoordinat.E.Y",
                        "Latitude",
                        "Longitude")


swemeta_trans_interpolated_clean <- swemeta_trans_interpolated[,!colnames(swemeta_trans_interpolated) %in% variables_to_remove]
swemeta_trans_interpolated_clean = cbind(swemeta_trans_interpolated_clean,
                                         swemeta_cleaned[1:91,c("Latitude", "Longitude")])

write.csv(swemeta_trans_interpolated_clean, file.path(sweden, "swemeta_trans_interpolated.csv"))

