library(phyloseq)
library(ggplot2)
library(decontam)
library(tibble)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

path_output <- args[1]
path_input <- args[2]
path_meta <- args[3]

setwd(path_output)
dir.create("output")

df <- read.table(file = path_input, sep = '\t', header = TRUE)
df.otu <- select(df, c(2,seq(4, ncol(df), by = 2))) %>% 
  column_to_rownames(var = "taxonomy_id") %>% 
  as.matrix()

df.taxa <- df %>% select(taxonomy_id) %>% rename(Species = taxonomy_id) %>% 
  mutate_if(is.integer, as.factor) %>% mutate(taxonomy_id = df$taxonomy_id) %>% 
  column_to_rownames(var = "taxonomy_id") %>% as.matrix()

df.meta <- read.table(file = path_meta, sep = '\t', header = TRUE)
df.meta<- df.meta %>% column_to_rownames(var = "sample_id")

phylo_otu<- otu_table(df.otu, taxa_are_rows = TRUE)
phylo_tax<- tax_table(df.taxa)
phylo_samples<- sample_data(df.meta)

phylo_object<- phyloseq(phylo_otu, phylo_tax, phylo_samples)

df <- as.data.frame(sample_data(phylo_object)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phylo_object)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
ggsave("output/library_size.png")

sample_data(phylo_object)$is.neg <- sample_data(phylo_object)$Sample_or_Control == "Control"
contamdf <- isContaminant(phylo_object, method="prevalence", neg="is.neg")
table(contamdf$contaminant)

rows_with_true <- row.names(contamdf)[which(contamdf$contaminant == TRUE)]
# Specify the file path and name
file_path <- "output/taxonomy_id_filter_list.txt"
# Save the row names to the text file
writeLines(rows_with_true, file_path)
