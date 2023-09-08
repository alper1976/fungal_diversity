##### Adapt_blast_lca.R ################
### IN SAGA

module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
R

.libPaths(c("/cluster/projects/nn9745k/RpackagCD its2es_4_0_0", .libPaths()))


# Load libraries
library(tidyverse)


# Set the working directory
#setwd("/cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/DADA2_extracted_samples_no_chim/ITS2/blast_lca_test")
setwd("/cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/DADA2_extracted_samples_no_chim/5_8S/blast_lca_test/")
path <- getwd()

# Load blast results
#blast_results_db <- read.table(file.path(path, "blast_results_ITS2")) # blast results from above
blast_results_db <- read.table(file.path(path, "blast_results_5_8S")) # blast results from above

colnames(blast_results_db) <- c( "Query ID",  "Subject", "Identity percentage", 
"Coverage", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
"S.start", "S.end", "Evalue", "Bitscore" ) 


summary(blast_results_db) #check if everything works
dim(blast_results_db)
length(unique(blast_results_db$Subject))
length(unique(blast_results_db$"Query ID"))

#tax = as.matrix(stringr::str_split_fixed(blast_results_clean$SubjectID, "[|;]", 9)) # to seperate 
tax_db = stringr::str_replace_all(blast_results_db$Subject, ";", " / ") # replaces ";" with "/", compatible with downstream analysis
tax_db = stringr::str_replace_all(tax_db, "k__", "")
tax_db = stringr::str_replace_all(tax_db, "p__", "")
tax_db = stringr::str_replace_all(tax_db, "c__", "")
tax_db = stringr::str_replace_all(tax_db, "o__", "")
tax_db = stringr::str_replace_all(tax_db, "f__", "")
tax_db = stringr::str_replace_all(tax_db, "g__", "")
tax_db = stringr::str_replace_all(tax_db, "s__", "")

taxonomy = gsub(".*[|]", "", tax_db) # removes "|" from header, replaces with nothing, to make compatible with downstream analysis

blast_results_db$"Taxonomy" = taxonomy

tax_acc = unique(blast_results_db$Subject)

tax_accno1 = str_split_fixed (tax_acc, "[|]", 3)[,2]  #extracts the accession number from subject id

####tax_accno = gsub(".*[|]([^|])[|].*", "\\1", tax_acc) # previous script from Eivind
####unlist(lapply(strsplit(as.character(tax_acc), split ='[|]'), "[[", 2)) ## from Franck

tax_accno2 = data.frame(cbind(as.data.frame(tax_acc), tax_accno1))

colnames(tax_accno2) = c("Subject", "Subject accession") #accession still in subject

#fix_ query and subject need accession in them.

blast_results_db = left_join(blast_results_db, tax_accno2, by = "Subject")
blast_results_db$"Subject Taxonomy ID" = blast_results_db$"Subject"
blast_results_db$"Source" = "UNITE_all"


blast_results_db_final = blast_results_db %>% select("Query ID", "Subject", "Subject accession","Subject Taxonomy ID", "Identity percentage", "Coverage","Evalue","Bitscore","Source","Taxonomy") 

                                        
colnames(blast_results_db_final) = c("#Query ID",
                                     "#Subject",
                                     "#Subject accession",
                                     "#Subject Taxonomy ID",
                                     "#Identity percentage",
                                     "#Coverage",
                                     "#Evalue",
                                     "#Bitscore",
                                     "#Source",
                                     "#Taxonomy"
                                          )

length(unique(blast_results_db_final$"#Subject"))
length(unique(blast_results_db_final$"#Query ID"))


#write.table(blast_results_db_final,file.path(path, "blast_results_for_lca_ITS2"), row.names = FALSE, quote = FALSE, sep = "\t") # should be used for python script underneath.
write.table(blast_results_db_final,file.path(path, "blast_results_for_lca_5_8S"), row.names = FALSE, quote = FALSE, sep = "\t") # should be used for python script underneath.


################# end ###########