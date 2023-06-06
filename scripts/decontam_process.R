library(stringr)
suppressPackageStartupMessages(library(furrr))
library(phyloseq)
library(tictoc)
library(purrr)
library(ggplot2)
library(decontam)
library(tibble)
library(tidyr)
suppressPackageStartupMessages(library(dplyr))
library(optparse)
suppressPackageStartupMessages(library(argparser))
options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)

# Function to parse tuple-like argument
parseTupleArgument <- function(arg) {
  values <- strsplit(arg, ",")[[1]]
  values <- as.numeric(values)
  return(values)
}

path_dir <- args[1]
path_otu <- args[2]
path_meta <- args[3]
dtype <- args[4]
thresh_interval <- parseTupleArgument(args[5])
workers <- args[6]


#options <- list(
#    make_option(c("-d", "--path_dir"), dest = "path_dir", default = getwd(),help = "Output Directory", type = "character"),
#    make_option(c("-o", "--path_otu"), dest = "path_otu", default = NULL,help = "Path to OTU table", type = "character"),
#    make_option(c("-m", "--path_meta"), dest = "path_meta", default = NULL, help = "Path to Meta table", type = "character"),
#    make_option(c("-t", "--dtype"), dest = "path_otu", default = NULL,help = "OTU file type (options: kraken, xtree)", type = "character"),
#    make_option(c("-i", "--thresh_interval"), dest = "thresh_interval", default = .1,.5,help = "Threshold Interval, seperated by ,", type = "character"),
#    make_option(c("-w", "--workers"), dest = "workers", default = 5, help = "Number of Cores", type = "integer")
#)

#parser <- add_option(parser, c("-h", "--help"), action = "help", type = "character",
#                    help = "Refer to https://github.com/qiutip/decontam_processing")

# Parse the command-line arguments
#opt <- parse_args(OptionParser(option_list = options))

#path_dir <- opt$path_dir
#path_otu <- opt$path_otu
#path_meta <- opt$path_meta
#dtype <- opt$dtype
#thresh_interval <- parseTupleArgument(opt$thresh_interval)
#workers <- opt$workers


setwd(path_dir)
dir.create("output")

#### FUNCTIONS ####

process_dataframes <- function(df.otu_path, df.meta_path, dtype){
    #read in data
    df.otu <- read.table(file = df.otu_path, sep = '\t', header = TRUE)
    df.meta <- read.table(file = df.meta_path, sep = '\t', header = TRUE)
    
    #tax_ranks
    ranks = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
    ranks_complete <- list()
    
    if (dtype == "kraken" | dtype == "kraken2" | dtype == "bracken"){
        kraken_rank = c("D", "K","P","C","O", "F", "G","S", "T")
        
        for (x in seq(1,length(ranks))){
            if (any(grepl(kraken_rank[x], df.otu$taxonomy_lvl)) == TRUE){
                ranks_complete[x] <- ranks[x]
            }
            }
        ranks_complete <- compact(ranks_complete)

        rownames(df.otu) <- df.otu$name
        df.otu <- select(df.otu, c(seq(5, ncol(df.otu), by = 2)))
        df.otu$sample_id <- rownames(df.otu)

        df.taxa<- df.otu %>% select(sample_id) %>% rename_with(~unlist(ranks_complete), sample_id)
        df.taxa<- df.taxa %>% mutate_if(is.character, as.factor)
    } else if (dtype == "xtree" | dtype =="Xtree") {
        #process df.otu
        df.otu$sample_id <- rownames(df.otu)

        xtree_rank = c("d__", "k__","p__","c__","o__", "f__", "g__","s__", "t__")
        
        ## Checking taxanomic ranks within otu
        for (x in seq(1,length(ranks))){
            if (any(grepl(xtree_rank[x], df.otu$sample_id)) == TRUE){
                ranks_complete[x] <- ranks[x]
            }
            }
        ranks_complete <- compact(ranks_complete)
        
        df.otu$sample_id <- rownames(df.otu)
        df.otu <- df.otu %>% mutate(sample_id = str_remove_all(sample_id, "d__|k__|p__|c__|o__|f__|g__|s__|t__"))
        rownames(df.otu) <- df.otu$sample_id
        

        df.taxa<- df.otu %>% select(sample_id) %>% separate(sample_id, unlist(ranks_complete),";")
        df.taxa<- df.taxa %>% mutate_if(is.character, as.factor)
        df.taxa["Unknown", ] <- "Unknown"
        } else if (dtype == "metaphlan" | dtype == "metaplhan"){
        #process df.otu
        
        rownames(df.otu) <- df.otu$clade_name
        df.otu$sample_id <- rownames(df.otu)
        
        metaphlan_rank = c("d__", "k__","p__","c__","o__", "f__", "g__","s__", "t__")

        ## Checking taxanomic ranks within otu
        
        for (x in seq(1,length(ranks))){
            if (any(grepl(metaphlan_rank[x], df.otu$clade_name)) == TRUE){
                ranks_complete[x] <- ranks[x]
            }
            }
        df.otu <- df.otu[grepl(metaphlan_rank[length(ranks_complete)],df.otu$clade_name), ]
        ranks_complete <- compact(ranks_complete)
        
        df.otu <- df.otu %>% mutate(sample_id = str_remove_all(sample_id, "d__|k__|p__|c__|o__|f__|g__|s__|t__"))
        rownames(df.otu) <- df.otu$sample_id

        df.taxa<- df.otu %>% select(sample_id) %>% separate(sample_id, unlist(ranks_complete),"\\|")
        df.taxa<- df.taxa %>% mutate_if(is.character, as.factor)
        
        df.otu <- select(df.otu, c(-clade_name))
        
        } else {
            print("Not the correct dtype, type again")
            q()
        }
    

    #process df.meta
    df.meta<- df.meta %>% column_to_rownames(var = "sample_id")

    #process misc
    
    df.otu <- select(df.otu, c(-sample_id))
    otu_mat<- as.matrix(df.otu)
    tax_mat<- as.matrix(df.taxa)

    phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
    phylo_TAX<- tax_table(tax_mat)
    phylo_samples<- sample_data(df.meta)
    
    ps<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples)
    sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
    
    return(ps)
}

contam_space <- function(ps, thresh) {
  # Perform some operations based on the threshold
  contamdf <- isContaminant(ps, method="prevalence", neg="is.neg", threshold = thresh)
  return(contamdf)
}

process <- function(ps, threshold_values){
  complete_ranks <- colnames(tax_table(ps))
  print(complete_ranks)
  
  results <- list() 
  contam_count <- list()
  for (x in seq(1,length(threshold_values))) {
    print(paste0("THRESHOLD VALUE = ", threshold_values[x]))
    
    ## Calling prevalence function
    result <- contam_space(ps, threshold_values[x])
    results[[x]] <- result 
    
    ## Counting number of contanimated taxa
    contam_count[[x]] <- sum(table(result$contaminant)[2])
    
    ## Creating TSV file with Filtered out Taxa based on thresholds
    rows_with_true <- row.names(result)[which(result$contaminant == TRUE)]
    df <- data.frame("taxonomy" = rows_with_true)
    #df$species <- lookup$name[lookup$taxonomy_id %in% df$tax_id]
    file_path <- paste0("output/taxonomy_id_filter_list", paste0(".thresh-", as.character(threshold_values[x])),".tsv")
    write.table(df, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    ## Updating Phyloseq Object with Removed Taxa
    print("Removed Contaminant Lowest Taxanomic Unit")
    update_ps = prune_taxa(rows_with_false <- row.names(result)[which(result$contaminant == FALSE)], ps)
    
    print("Processing Abundance Matrices")
    #Re-formating datamatrix with absolute abundance
    
    for (y in complete_ranks){
      
      df_relative <- ps %>% tax_glom(taxrank = y) %>%
        transform_sample_counts(function(x) {x/sum(x)})%>% psmelt() %>% 
        select(y, Sample, Abundance) %>% spread(Sample, Abundance)
      
      #updates threshold/taxa_rank
      file_path_r <- paste0("output/", paste0("relative_abundance.", y), 
                            paste0(".thresh-", as.character(threshold_values[x])),".tsv")
      
      ##writes dataframes
      write.table(df_relative, file = file_path_r, sep = "\t", quote = F, row.names = F, col.names = T)
    }
    file_path_all <- paste0("output/relative_abundance.all", 
                            paste0(".thresh-", as.character(threshold_values[x])),".tsv")
    write.table(ps %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% arrange(OTU) %>%
                  select(OTU, unlist(complete_ranks), Sample, Abundance) %>%
                  spread(Sample, Abundance), file = file_path_all, 
                sep = "\t", quote = F, row.names = F, col.names = T)
  }
  
  print("Contanimated Count vs Threshold")
  contam_count <- lapply(contam_count, function(x) replace(x, is.na(x), 0))
  df_count <- data.frame(thresh = unlist(threshold_values), count_contaminants = unlist(contam_count))
  write.table(df_count, 
              file = "output/thresh_vs_count.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
  ggplot(data=df_count, aes(x=threshold_values, y=count_contaminants)) + geom_point()
  ggsave("output/thresh_vs_count.png")
}

process_parallel <- function(ps, threshold_values, workers = 5){
    
    plan("multisession", workers=workers)

    complete_ranks <- colnames(tax_table(ps))
    print(complete_ranks)
    
    results <- list() 
    contam_count <- list()
    contam_count <- future_map(1:length(threshold_values), ~ {
        
        print(paste0("THRESHOLD VALUE = ", threshold_values[.x]))
        
        ## Calling prevalence function
        result <- contam_space(ps, threshold_values[.x])
        results[[.x]] <- result 

        ## Counting number of contanimated taxa
        contam_count[[.x]] <<- sum(table(result$contaminant)[2])

        ## Creating TSV file with Filtered out Taxa based on thresholds
        rows_with_true <- row.names(result)[which(result$contaminant == TRUE)]
        df <- data.frame("taxonomy" = rows_with_true)
        file_path <- paste0("output/taxonomy_id_filter_list", paste0(".thresh-", as.character(threshold_values[.x])),".tsv")
        write.table(df, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE)

        ## Updating Phyloseq Object with Removed Taxa
        print("Removed Contaminant Lowest Taxanomic Unit")
        update_ps = prune_taxa(rows_with_false <- row.names(result)[which(result$contaminant == FALSE)], ps)

        print("Processing Abundance Matrices")
        
        #Re-formating datamatrix with absolute abundance

        for (y in complete_ranks){

            df_relative <- ps %>% tax_glom(taxrank = y) %>%
            transform_sample_counts(function(x) {x/sum(x)})%>% psmelt() %>% 
            select(all_of(y), Sample, Abundance) %>% spread(Sample, Abundance)
            
            #updates threshold/taxa_rank
            file_path_r <- paste0("output/", paste0("relative_abundance.", y), 
                paste0(".thresh-", as.character(threshold_values[.x])),".tsv")
            
            ##writes dataframes
            write.table(df_relative, file = file_path_r, sep = "\t", quote = F, row.names = F, col.names = T)
        }
        file_path_all <- paste0("output/relative_abundance.all", 
                              paste0(".thresh-", as.character(threshold_values[.x])),".tsv")
        write.table(ps %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% arrange(OTU) %>%
            select(OTU, unlist(complete_ranks), Sample, Abundance) %>%
            spread(Sample, Abundance), file = file_path_all, 
            sep = "\t", quote = F, row.names = F, col.names = T)
        
        contam_count[[.x]]
    })
    plan("sequential")
    
    print("Contanimated Count vs Threshold")
    contam_count <- lapply(contam_count, function(x) replace(x, is.na(x), 0))
    df_count <- data.frame(thresh = unlist(threshold_values), count_contaminants = unlist(contam_count))
    write.table(df_count, 
                file = "output/thresh_vs_count.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
    ggplot(data=df_count, aes(x=threshold_values, y=count_contaminants)) + geom_point()
    ggsave("output/thresh_vs_count.png")
    }


ps <- process_dataframes(path_otu, path_meta, dtype)
threshold_values<- seq(from = thresh_interval[1], to = thresh_interval[2], length.out = 5)
process_parallel(ps, threshold_values, workers)

