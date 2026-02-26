library(vegan)
library(permute)
#install.packages('devtools')
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

Permanova <-function(ps, method){
  data_prop_bray <- phyloseq::distance(ps, method = method)
  sampledf <- data.frame(sample_data(ps))
  ht_well <- adonis2(data_prop_bray ~ research, data = sampledf)
  ps_dist_matrix <- phyloseq::distance(ps, method=method)
  permanova_pairwise <- vegan::adonis2(ps_dist_matrix ~ phyloseq::sample_data(ps)$research, permutations = 9999)
  betadisper(ps_dist_matrix, sample_data(ps)$research)
  return(c(permanova_pairwise, permanova_pairwise$R2[1]))
}


Anosim <- function(ps, method){
  metadata <- data.frame(sample_data(ps))
  anosim_test<-anosim(phyloseq::distance(ps, method=method), 
                             metadata$research, permutations = 9999)
  return(c(anosim_test$signif, anosim_test$statistic))
}
