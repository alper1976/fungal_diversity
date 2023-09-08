# Script for merging two taxonomies



################## Work in progress ################################

## saga
#module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0

# load(file=file.path(figsPath, "stats.RData"))

library(readxl)
library(tidyverse)


.libPaths(c("/cluster/projects/nn9745k/Rpackages_4_0_0", .libPaths()))
rootPath  = file.path("/cluster/projects/nn9744k/02_results/29_euk_20220325/Pedro/fungal_pipeline_final_classification")

##load packages

tax_5_8S = read.table(file = file.path(rootPath, "5_8S", "lca_results_5_8S"), sep = "\t")
colnames(tax_5_8S) = c("Query",
                       "rank_5_8S",
                       "taxon_5_8S",
                       "kingdom_5_8S",
                       "phylum_5_8S",
                       "class_5_8S",
                       "order_5_8S",
                       "family_5_8S",
                       "genus_5_8S",
                       "species_5_8S",
                       "method_5_8S"
                        )
tax_5_8S$Query <- sub(";.*", "", tax_5_8S$Query)

tax_ITS2 = read.table(file = file.path(rootPath, "ITS2", "lca_results_ITS2"), sep = "\t")
colnames(tax_ITS2) = c("Query",
                       "rank_ITS2",
                       "taxon_ITS2",
                       "kingdom_ITS2",
                       "phylum_ITS2",
                       "class_ITS2",
                       "order_ITS2",
                       "family_ITS2",
                       "genus_ITS2",
                       "species_ITS2",
                       "method_ITS2"
                        )
tax_ITS2$Query <- sub(";.*", "", tax_ITS2$Query)

otus = read_excel(file.path(rootPath, "OTU_table_lulu_curated.xlsx"))


merged_tax = full_join(tax_5_8S, tax_ITS2, by="Query")

final_tax = merged_tax
final_tax[is.na(final_tax)] = "no identification"

for(i in 1:nrow(final_tax)) {
  line = final_tax[i,]
  if ( (line$rank_ITS2 == "no identification") &&
    (line$kingdom_5_8S != "no identification") ) {
    final_tax[i, 12:21] = line[2:11]
  }
  else if ( (line$kingdom_ITS2 == "no identification") &&
    (line$kingdom_5_8S != "no identification") ) {
    final_tax[i, 12:21] = line[2:11]
  }
  else if ( (line$phylum_ITS2 == "no identification") &&
    (line$phylum_5_8S !="no identification") ) {
    final_tax[i, 12:21] = line[2:11]
  }
  else if ( (line$class_ITS2 == "no identification") &&
    (line$class_5_8S != "no identification") ) {
    if(line$phylum_ITS2 == line$phylum_5_8S) {
    final_tax[i, 12:21] = line[2:11]
    }
    else {
      final_tax[i, "Match"] = "No class"
    }
  }
  else if ( (line$order_ITS2 == "no identification") &&
    (line$order_5_8S != "no identification") ) {
    if(line$class_ITS2 == line$class_5_8S) {
    final_tax[i, 12:21] = line[2:11]
    }
    else {
      final_tax[i, "Match"] = "No order"
    }
  }
  else if ( (line$family_ITS2 == "no identification") &&
    (line$family_5_8S != "no identification") ) {
    if(line$order_ITS2 == line$order_5_8S) {
    final_tax[i, 12:21] = line[2:11]
    }
    else {
      final_tax[i, "Match"] = "No family"
    }
  }
  else if ( (line$genus_ITS2 == "no identification") &&
    (line$genus_5_8S != "no identification") ) {
    if(line$family_ITS2 == line$family_5_8S) {
    final_tax[i, 12:21] = line[2:11]
    }
    else {
      final_tax[i, "Match"] = "No genus"
    }
  }
  else if ( (line$species_ITS2 == "no identification") &&
    (line$species_5_8S != "no identification") ) {
    if(line$genus_ITS2 == line$genus_5_8S) {
    final_tax[i, 12:21] = line[2:11]
    }
    else {
      final_tax[i, "Match"] = "No species"
    }
  }
  else {
   final_tax[i, "Match"] = "No identification"
  }
}

final_tax = final_tax[c(1, 12:22)]

colnames(final_tax) = c("Query",
                       "rank",
                       "taxon",
                       "kingdom",
                       "phylum",
                       "class",
                       "order",
                       "family",
                       "genus",
                       "species",
                       "method",
                       "Match"
                        )

taxonomy = final_tax[,c("kingdom",
                       "phylum",
                       "class",
                       "order",
                       "family",
                       "genus",
                       "species")]
row.names(taxonomy) = final_tax$Query


