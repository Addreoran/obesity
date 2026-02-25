

library(phyloseq)


read_and_parse_data <- function(metadata_path, otu_path, tax_path) {
  
  otu <- read.csv(otu_path, row.names = 1) 
  otu <- otu_table(as.matrix(otu), taxa_are_rows = TRUE) 
  colnames(otu)<-sub("^X", "", colnames(otu))
  
  metadata <- read.csv(metadata_path, row.names = 1)
  metadata <- sample_data(metadata) 
  rownames(metadata)<-metadata$id
  
  tax <- read.csv(tax_path, row.names = 1) 
  tax <- tax_table(as.matrix(tax)) 
  
  ps <- phyloseq(otu, metadata, tax) 
  return(ps) 
}
