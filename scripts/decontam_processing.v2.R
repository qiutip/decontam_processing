library(phyloseq)
library(ggplot2)
library(decontam)
library(tibble)
library(tidyr)
library(dplyr)

path_output <- args[1]
path_input <- args[2]
path_meta <- args[3]

setwd(path_output)
dir.create("output")

### FUNCTIONS ###
contam_space <- function(phylo_object, thresh) {
  # Perform some operations based on the threshold
  contamdf <- isContaminant(phylo_object, method="prevalence", neg="is.neg", threshold = thresh)
  return(contamdf)
}

#### READING IN DATA ###
df <- read.table(file = path_input, sep = '\t', header = TRUE)
df.otu <- select(df.otu, c(2,seq(4, ncol(df.otu), by = 2))) %>% 
  column_to_rownames(var = "taxonomy_id") %>% 
  as.matrix()

df.taxa <- df %>% select(taxonomy_id) %>% rename(Species = taxonomy_id) %>% 
  mutate_if(is.integer, as.factor) %>% mutate(taxonomy_id = df_bs$taxonomy_id) %>% 
  column_to_rownames(var = "taxonomy_id") %>% as.matrix()

df.meta <- read.table(file = path_meta, sep = '\t', header = TRUE)
df.meta<- df.meta %>% column_to_rownames(var = "sample_id")

### CREATING PHYLOSEQ OBJECT ###
phylo_otu<- otu_table(df.otu, taxa_are_rows = TRUE)
phylo_tax<- tax_table(df.taxa)
phylo_samples<- sample_data(df.meta)

phylo_object<- phyloseq(phylo_otu, phylo_tax, phylo_samples)
sample_data(phylo_object)$is.neg <- sample_data(phylo_object)$Sample_or_Control == "Control"


# Define a sequence of threshold values
#threshold_values <- seq(.01, .5, by=.05)
threshold_values <- c(.1, .2,. .3, .4 .5)


### IDENTIFYING AND REMOVING CONTANIMATED TAXA ###

results <- list() 
contam_count <- list()
for (x in seq(1,length(threshold_values))) {
  ## Calling prevalence function
  result <- contam_space(phylo_object, threshold_values[x])
  results[[x]] <- result 
  
  ## Counting number of contanimated taxa
  contam_count[[x]] <- sum(table(result$contaminant)[2])
  
  ## Creating TSV file with Filtered out Taxa based on thresholds
  rows_with_true <- row.names(result)[which(result$contaminant == TRUE)]
  df <- data.frame("tax_id" = rows_with_true)
  df$species <- lookup$name[lookup$taxonomy_id %in% df$tax_id]
  file_path <- paste0("output/taxonomy_id_filter_list", paste0(".thresh-", as.character(threshold_values[x])),".txt")
  write.table(df, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  ## Updating Phyloseq Object with Removed Taxa
  print("Removed Contaminant Taxa")
  update_ps = prune_taxa(rows_with_false <- row.names(result)[which(result$contaminant == FALSE)], phylo_object)
  
  file_path_r <- paste0("output/relative_abundance.species", paste0(".thresh-", as.character(threshold_values[x])),".tsv")
  
  print("Processing Abundance Matrices")
  
  #Reformating Relative Abundance Maxtrix across Samples
  df_relative <- update_ps %>% tax_glom(taxrank = "Species") %>%
    transform_sample_counts(function(x) {x/sum(x)})%>% psmelt() %>% 
    select(Species, Sample, Abundance) %>% spread(Sample, Abundance) %>%  mutate(Species = case_when(
      Species %in% lookup$taxonomy_id ~ lookup$name[match(Species, lookup$taxonomy_id)],
      TRUE ~ Species))
  
  write.table(df_relative, 
              file = file_path_r, sep = "\t", quote = F, row.names = F, col.names = T)
}

### CREATING THRESHOLD VS CONTANIMATED TAXAS ###
contam_count <- lapply(contam_count, function(x) replace(x, is.na(x), 0))
df_count <- data.frame(thresh = unlist(threshold_values), count_contaminants = unlist(contam_count))
write.table(df_count, 
            file = "output/thresh_vs_count.species.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
ggplot(data=df_count, aes(x=threshold_values, y=count_contaminants)) + geom_point()
ggsave("output/thresh_vs_count.png")
