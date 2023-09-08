# stats
These are scripts (R and phyton scripts) used to run the statistical analyses

## individual scripts
### danube_indicator_stats.R
This contains the code to produce the boxplots for the environmental variables and the graphs for the microbial community analysis. 

### Spatiotemporal_ASV_screening.py
The file 'Spatiotemporal_ASV_screening.py' was run without issues in a jupyter notebook on an Anaconda distribution of python 3.6.
It takes as input files 'j3s_L.df.export.tsv', 'j3s_M.df.export.tsv', 'j3s_R.df.export.tsv' provided in the 'Data' folder. These files are the ICDPR data used for this manuscript.
The script outputs the files 'Danube_WQ_predictors.txt', 'Danube_chl_a_predictors.txt', 'Danube_WQ_dictionaries.txt' and 'Danube_chl_a_dictionaries.txt'.
The latter outputs are used as inputs for ASV informational redundancy analyses. 

### ASV_screening_output_analysis.R
The file 'ASV_screening_output_analysis.R' was run on R 3.6.2 without issues.
It processes the outputs of 'Spatiotemporal_ASV_screening.py' to yield the co-occurence and co-exclusion networks used to assess ASV information content and redundancy.
