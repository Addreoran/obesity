

library(phyloseq)


read_and_parse_data<-funtion(metadata_path, otu_path, tax_path){
  metadata <- read.csv(metadata_path)
  otu <- read.csv(otu_path)
  tax <- read.csv(tax_path)
  ps<-phyloseq(otu, metadatam tax)
  return(ps)
}
