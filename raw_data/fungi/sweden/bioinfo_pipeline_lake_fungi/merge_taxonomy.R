####### Script for merging two taxonomies and preparing the needed files

## In Saga
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
R
.libPaths(c("/cluster/projects/nn9745k/Rpackages_4_0_0", .libPaths()))

##load packages
library(readxl)
library(tidyverse)
library(openxlsx)
library(phyloseq)

# Path
rootPath  = file.path("/cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/final_tax_classification")


tax_5_8S = read.table(file = file.path(rootPath, "blast_results_lca_8_80_85_only_lca_flh_5_8S"), sep = "\t")
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
tax_5_8S$Query <- sub(";.*", "", tax_5_8S$Query) # Maybe this (delete the size of the OTU names/Query is problematic)


tax_ITS2 = read.table(file = file.path(rootPath, "blast_results_lca_8_80_85_only_lca_flh_ITS2"), sep = "\t")
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
tax_ITS2$Query <- sub(";.*", "", tax_ITS2$Query)# delete the size of the OTU names/Query


merged_tax = full_join(tax_5_8S, tax_ITS2, by="Query")


final_tax = merged_tax

#convert all factors to characters to avoid warning related to NA
final_tax <- data.frame(lapply(final_tax, as.character), stringsAsFactors=FALSE)

final_tax[is.na(final_tax)] = "no identification" # transform NA in "no identification"# but not all NA were transformed initially... 



for(i in 1:nrow(final_tax)) {
  line = final_tax[i,]
  if ( (line$kingdom_ITS2 == "no identification") &&
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



#load the OTU table
otus = read.xlsx(file.path(rootPath, "OTU_table_lulu_curated.xlsx"))

# Merge OTU table (lulu curated) with taxonomy ## Note there are blanks->NAslacking the taxonomic assignments
taxonomy2 = taxonomy %>% rownames_to_column("OTUid")
otus_tax = left_join(otus, taxonomy2)

# Transform NAs to the string "no identification" in the taxonomy columns
otus_tax2 = otus_tax %>% replace_na(list(kingdom='no identification',
                                         phylum='no identification',
                                         class='no identification',
                                         order='no identification',
                                         family='no identification',
                                         genus='no identification',
                                         species='no identification'))

## save LCA merged taxonomy
write.xlsx(taxonomy2, "/cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/final_tax_classification/lca_merged_taxonomy.xlsx")


## save final OTU table with taxonomy 
write.xlsx(otus_tax2, "/cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/final_tax_classification/OTU_table_lulu_curated_lca.xlsx")


######### Create the phyloseq objects

#Adapt OTU table
otus2 = otus %>% remove_rownames()%>% column_to_rownames('OTUid')

otus3 <- otu_table(otus2, taxa_are_rows = TRUE)
dim(otus3)

#Adapt the taxonomy - using that merged to the OTU table.  
taxonomy3 = otus_tax2 %>% select(OTUid, c(77:83))
taxonomy4 =taxonomy3 %>% column_to_rownames('OTUid')

taxonomy5 <- tax_table(taxonomy4)

taxa_names(taxonomy5)
taxa_names(taxonomy5) <- taxonomy3$OTUid

# CREATE PHYLOSEQ OBJECT (OTUTABLE WITH TAXONOMY)
sweden_physeq <- phyloseq(otus3, taxonomy5) 


#save 
saveRDS(sweden_physeq,"/cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/final_tax_classification/phyloseq_sweden.rds") ####### phyloseq object to use for further analysis on local R.
