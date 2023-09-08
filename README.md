# fungal_diversity

Here, we use high-resolution 5.8s gene and ITS2 sequencing to explore fungal diversity in two Scandinavian lake datasets and resolve fungal community assembly mechanisms across large-scale latitudinal, climatic (nemoral to arctic), and nutrient (oligotrophic to eutrophic) gradients. This repository includes scripts to run cutadapt (https://github.com/marcelm/cutadapt), dada2 (https://github.com/benjjneb/dada2) and statistical analyses for studying community assembly rules in aquatic fungi.

## pull repository

```
cd path to repositories
git clone github.com/alper1976/fungal_diversity.git
```

## Authors and acknowledgment
Scripts were written by Pedro Maria Martin Sanchez, Even Werner; Laurent Fontaine and Alexander Eiler.

## License
This Code is subject to the terms of the MIT License. 

## Project status
Results from this project have been submitted to a peer-reviewed scientifc journal.

## Folders and code
The analyses code is divided into two folders "raw_data" and "stats" representing the code to analyze raw sequencing data and perform statistical analysis, respectively. The "Data" folder contains the data input and output used in statistical analyses and data visualization.

### metadata
In this folder you can find the metadata for both the Norwegian and Swedish datasets.
'100lakes_norway_meta_sample.xlsx' represents the Norwegian metadata. '100lakes_metadata_final.xlsx' represents the Swedish metadata. There are also '.csv' files. The latter represents data from the Trend Lake national monitoring program that were used in this manuscript. For more detail see https://www.havochvatten.se/overvakning-och-uppfoljning/miljoovervakning/organisation-och-programomraden/miljoovervakningens-programomrade-sotvatten/delprogram-trendstationer-sjoar.html#h-Datavardarochutforare

### data
Here you can find all code for running statistical analyses as well as the Amplicon sequence variants (ASV) tables used for statistical analyses as well as the taxonomic annotations files. In addition subfolders including code, input and output files from various analyses are provided.

### raw_data
Here you can find the code to obtain the Amplicon sequence variants (ASV) tables from the raw sequencing data.

## released version 1.0.0

