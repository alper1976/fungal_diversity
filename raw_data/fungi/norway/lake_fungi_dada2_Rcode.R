
# Set the library path 
.libPaths(c("/cluster/projects/nn9745k/Rpackages_4_0_0", .libPaths()))

##load packages
library(optparse)
library(dada2)
library(phyloseq)
library(vegan)

library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")


## set parametersR

option_list = list(

   # make_option(c("-f", "--truncate_length_fwd"), # can be altered
    #    type="integer",
     #   default=120,
      #  help="Truncate length of forward reads [default %default]",
       # metavar="number"),

   # make_option(c("-r", "--truncate_length_rev"), # can be altered
    #    type="integer",
     #   default=80,
      #  help="Truncate length of reverse reads [default %default]",
       # metavar="number"),

    make_option(c("-q", "--truncate_quality"), #try with Phred quality 30, suggested Miya et al 2015 
        type="integer",
        default=20,
        help="Truncates based on quality score [default %default]",
        metavar="character"),

    make_option(c("-p", "--path"), 
        type="character",
        default="/cluster/projects/nn9744k/02_results/29_euk_20220325/Pedro/AdaptersRemoved",   #path to adapters_removed
        help="path to input",
        metavar="character"),

    make_option(c("-o", "--output"),
        type="character",
        default="/cluster/projects/nn9744k/02_results/29_euk_20220325/Pedro/",
        help="path to output",
        metavar="character"),

  #  make_option(c("-s", "--path_to_database"), #path to scandifish
   #     type="character",
    #    default="/cluster/projects/nn9745k/03_databases/fish/ScandiFish_12s_v1.1/ScandiFish12s_v1.1_db",                      ### location of database
     #   help="path to database",
      #  metavar="character"),

    make_option(c("-n", "--max_n"),
        type="integer",
        default=0,
        help="maximum of ambiguous bases (Ns) [default %default]",
        metavar="number")


    )

opt = parse_args(OptionParser(option_list=option_list))



##################################


path_to_input = opt$path
path_to_output = opt$output

path_to_fwd = list.files(path_to_input, pattern = "R1_001.fastq.gz", all.files = TRUE, full.names = TRUE, recursive = FALSE)
path_to_rev = list.files(path_to_input, pattern = "R2_001.fastq.gz", all.files = TRUE, full.names = TRUE, recursive = FALSE)



##################################################################
# Identify primers - for FUNGI's workflow, but in our case the primers were removed before...
#####################################################################

##### look at the primers in the reads ###########

FWD <-"CAWCGATGAAGAACGCAG"    #INSERT_FWD_PRIMER
REV <- "TCCTCCGCTTAYTRATATGC"    #INSERT_REV_PRIMER

allOrients <- function(primer) { # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients


### I have not removed the ambiguous bases (Ns)


#We are now ready to count the number of times the primers appear in the forward and reverse read, 
#while considering all possible primer orientations. 
#Identifying and counting the primers on one set of paired end FASTQ files is sufficient, 
# assuming all the files were created using the same library preparation, so we’ll just process the first sample.

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = path_to_fwd[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = path_to_rev[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = path_to_fwd[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = path_to_rev[[1]]))



    #####################################################################
    # clean and match
    #####################################################################

    ## Removes all sequences containing any ambiguous bases (N).
    ## Tries to match up all pairs of forward and reverse reads by ID,
    ## Removes all sequences that have been orphaned by filtering.
    ## Creates new directory for clean and matched files

    fastq_fwd <- basename(path_to_fwd)
    fastq_rev <- basename(path_to_rev)
    filtered_path <- file.path(path_to_output, "filtered")
    filtered_fwd <- file.path(filtered_path, fastq_fwd)
    filtered_rev <- file.path(filtered_path, fastq_rev)

    ##Plot quality profiles


    plot_fwd_qual <- plotQualityProfile(file.path(path_to_fwd)[1:10])
    pdf(file.path(path_to_output,"plot_fwd_qual.pdf"))
    plot_fwd_qual
    dev.off()

    plot_rev_qual <- plotQualityProfile(file.path(path_to_rev)[1:10])
    pdf(file.path(path_to_output, "plot_rev_qual.pdf"))
    plot_rev_qual
    dev.off()


  if(length(path_to_fwd) != length(path_to_fwd))
        # DevNote: Returns error to stderr
        write("The number of input files looks dodgy", stderr())

   write(file.path(path_to_fwd), stdout())
   write(file.path(filtered_fwd), stdout())
   write(file.path(path_to_rev), stdout())
   write(file.path(filtered_rev), stdout())


####### Filtering: these parameters need to be set by user
    out <- filterAndTrim(fwd = file.path(path_to_fwd), filt = file.path(filtered_fwd),
                       rev = file.path(path_to_rev), filt.rev = file.path(filtered_rev),
                       minLen = 50, ##### For ITS we delete that: "truncLen = c(opt$truncate_length_fwd, opt$truncate_length_rev)"; and adding "minLen" athough is not needed after cutadapt (min lenght = 200)
                       maxEE = 2,
                       truncQ = opt$truncate_quality, 
                       maxN = opt$max_n,
                       rm.phix=TRUE,
                       compress=TRUE,
                       verbose=TRUE,
                       multithread=FALSE)

    if (length(path_to_fwd)*2 != length(list.files(filtered_path))) {
        write(paste0("Filtering did not result in any sequences from some files.\n",
                     "The number of filtered sequence files is ",
                     length(list.files(filtered_path)),
                     " \n",
                     "while the number of raw files is ",
                     length(path_to_fwd) + length(path_to_rev),
                     " \n"),
                     stdout())
    }

  # Infer Sequence variance

    filtered_fwd = list.files(filtered_path, pattern = "R1_001.fastq.gz", all.files = TRUE, full.names = TRUE, recursive = TRUE)
    filtered_rev = list.files(filtered_path, pattern = "R2_001.fastq.gz", all.files = TRUE, full.names = TRUE, recursive = TRUE)

    sample.names <- basename(filtered_fwd)  # Assumes filename = samplename_XXX.fastq.gz
    sample.names.rev <- basename(filtered_rev) # Assumes filename = samplename_XXX.fastq.gz
    names(filtered_fwd) <- sample.names
    names(filtered_rev) <- sample.names
    set.seed(100)

    # Learn forward error rates
    error_fwd <- learnErrors(filtered_fwd, nbases=1e9, randomize = TRUE, multithread=FALSE)
    # Learn reverse error rates
    error_rev <- learnErrors(filtered_rev, nbases=1e9, randomize = TRUE, multithread=FALSE)

    pdf(file.path(path_to_output, "plot_fwd_errors.pdf"))
    plotErrors(error_fwd, nominalQ = T)
    dev.off()

    pdf(file.path(path_to_output, "plot_rev_errors.pdf"))
    plotErrors(error_rev, nominalQ = T) #plotErrors(error_fwd)???
    dev.off()


    # Sample inference and merger of paired-end reads
    mergers <- vector("list", length(sample.names))
    names(mergers) <- sample.names
    dada_fwd <- vector("list", length(sample.names))
    dada_rev <- vector("list", length(sample.names))

        for(sam in sample.names) {
            cat("Processing:", sam, "\n")
            derepF <- derepFastq(filtered_fwd[[sam]])
            ddF <- dada(derepF, err = error_fwd, multithread=FALSE)
            dada_fwd[[sam]] <- ddF
            derepR <- derepFastq(filtered_rev[[sam]])
            ddR <- dada(derepR, err = error_rev, multithread=FALSE)
            dada_rev[[sam]] <- ddR
            concatenate <- mergePairs(ddF, derepF, ddR, derepR, justConcatenate=FALSE)
            mergers[[sam]] <- concatenate
        }
        seqtab <- makeSequenceTable(mergers)
        write(paste0("Merging the sequences resulted in ", ncol(seqtab), " sequences."))

        seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)  ### should be used for plotting taxonomy are done in file seq_data

    saveRDS(seqtab, file.path(path_to_output, "seqtab.rds"))
    seqtab.final <- seqtab.nochim

    saveRDS(seqtab.final, file.path(path_to_output, "seqtab_nochim.rds"))


    #####################################################################
    # track sequences through analysis
    #####################################################################
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out,
                   sapply(mergers, getN),
                   rowSums(seqtab.final))
    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    colnames(track) <- c("input", "filtered", "merged", "final")
    rownames(track) <- sample.names

    # write to file
    write.csv(track, file.path(path_to_output, "track_sequence_stats.csv"))

    write(paste0("Write tracking file to outp
	ut."), stdout())

    uniquesToFasta(seqtab.final, fout = file.path(path_to_output, "rep-seqs.fna"), ids=colnames(seqtab.final))


    


    save(file=file.path(path_to_output, "rrun_final.RData"))

    write(paste0("Analysis run to the end."), stdout())



######################### ADAPTATION OF AQUA-OMG PIPELINES #############

## Load libraries
library(tidyverse)

setwd("/cluster/projects/nn9744k/02_results/29_euk_20220325/Pedro/")
path <- getwd()


##### ASV TABLES FROM DADA2

# load the seqtab (652 ASVs)
seqtab = readRDS("seqtab.rds")

# load the seqtab_nochim (226 ASVs)
seqtab_nochim = readRDS("seqtab_nochim.rds")


# Edit/clean sample names (just those starting by L) 
seqtab_names = str_extract(rownames(seqtab_nochim), "(L)([^_]+)")
rownames(seqtab_nochim) = seqtab_names


###### From the OMG script:

#Transpose table, assign names, extract sequences and saving table, for further processing:
trasBoth_nochim_sumtable <- as.data.frame(t(seqtab_nochim))

#Get DNA sequences
sequences <- row.names(trasBoth_nochim_sumtable)

#Assign new rownames
row.names(trasBoth_nochim_sumtable) <- paste0("seq",seq.int(nrow(trasBoth_nochim_sumtable)))

tbname <- file.path(path,"DADA2_noChime.table")
{write.table(trasBoth_nochim_sumtable,tbname,sep="\t",col.names = NA, quote=FALSE)}

#Extract OTUs (sequences)
sinkname <- file.path(path,"DADA2_noChime_raw.otus")
{
  sink(sinkname)
  for (seqX in seq.int(nrow(trasBoth_nochim_sumtable))) {
    header <- paste0(">","seq",seqX,"\n")
    cat(header)
    seqq <- paste0(sequences[seqX],"\n")
    cat(seqq)
  }
  sink()
}

#Define function to extract sequences sample-wise
extrSamDADA2 <- function(my_table) {
  out_path <- file.path(path, "DADA2_extracted_samples_no_chim")
  if(!file_test("-d", out_path)) dir.create(out_path)
  for (sampleX in seq(1:dim(my_table)[1])){
    sinkname <- file.path(out_path, paste0(rownames(my_table)[sampleX],".fas"))
    {
      sink(sinkname)
      for (seqX in seq(1:dim(my_table)[2])) {
        if (my_table[sampleX,seqX] > 0) {
          header <- paste0(">",rownames(my_table)[sampleX],";size=",my_table[sampleX,seqX],";","\n")
          cat(header)
          seqq <- paste0(colnames(my_table)[seqX],"\n")
          cat(seqq)
        }
      }
      sink()
    }
  }
}

#Extract samplewise sequences from the non-chimera table using the above function:
extrSamDADA2(seqtab_nochim)


############## Generated OUTPUTS:

## folder with extracted reads in fasta files (one per sample): "DADA2_extracted_samples_no_chim" 
## "DADA2_noChime.table" where the rows are ASVs (Seq1-Seq226) and columns are samples with namees like: 25-L41_S25_L001.fastq.gz) 
## "DADA2_noChime.raw.otus": sequences in fasta format, like: >Seq1.... 


#### get and save read length  (not sure if the function getSequences is the loaded dad2 version?....)

reads.per.seqlen.no.chim <- tapply(colSums(seqtab_nochim), factor(nchar(getSequences(seqtab_nochim))), sum)
df.without.chim <- data.frame(length=as.numeric(names(reads.per.seqlen.no.chim)), count=reads.per.seqlen.no.chim)


write.xlsx(df.without.chim, "No_chimeras_seq_length_table.xlsx")

df.without.chim$count<-as.numeric(df.without.chim$count)
ggplot(data=df.without.chim, aes(x=length, y=count)) + geom_col()
ggsave("ASV_no_chime_length.pdf", width = 18, height = 14, units = c("in"))

################################



