source("./src/read_data.R")
source("./src/differential_bacteria.R")
source("./src/PCoA.R")


run <- function(metadata_path, otu_path, tax_path, suffix){
   ps <- read_and_parse_data(metadata_path, otu_path, tax_path)
   ps_genus <- tax_glom(ps, "genus")

  ##
  deseq_res <- deseq2_result(ps_genus)
  ancomb_res <- ancombc2_result(ps_genus)
  
  tax_tab <- as.data.frame(unclass(tax_table(ps_genus)))
  otu_tab <- as.data.frame(unclass(otu_table(ps_genus)))
  
  differential_save_path <- paste0("./result/differential_ancomb_deseq_",suffix,".csv")
  deseq2_ancombc2_result(deseq_res, ancomb_res, tax_tab, otu_tab, differential_save_path)
  
  ##
  save_path_PCoAAitch <- paste0("./result/PCoA_Aitch", suffix, ".svg")
  PCoAAitch(ps, save_path_PCoAAitch)
  save_path_PCoABray <- paste0("./result/PCoA_Bray", suffix, ".svg")
  PCoABray(ps, save_path_PCoABray)
  
  ##
  Anosim
  Permanova
  
  ##
  alfa_diversity
  
  ##
  tax_differences
}

# analyse by CAP with threshold 250
metadata_path <- "./data/metadata_by_metric.csv"
otu_path <- "./data/otu_by_metric.csv"
tax_path <- "./data/tax_data_by_metric.csv"

run(metadata_path, otu_path, tax_path, "by_metric")
