library(vegan)
library(permute)
#install.packages('devtools')
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

Permanova <-function(ps, method){
  ps_dist_matrix <- phyloseq::distance(ps, method=method)
  permanova_pairwise <- vegan::adonis2(ps_dist_matrix ~ phyloseq::sample_data(ps)$research, permutations = 9999)
  betadisper(ps_dist_matrix, sample_data(ps)$research)
  return(c(permanova_pairwise$`Pr(>F)`[1], permanova_pairwise$R2[1]))
}


Anosim <- function(ps, method){
  metadata <- data.frame(sample_data(ps))
  anosim_test<-anosim(phyloseq::distance(ps, method=method), 
                             metadata$research, permutations = 9999)
  return(c(anosim_test$signif, anosim_test$statistic))
}


PermanovaCLR <-function(ps, method){
  physeq_clr <- microbiome::transform(ps, "clr")
  otu_clr <- as.data.frame(otu_table(physeq_clr))
  dist_matrix <- vegan::vegdist(otu_clr, method = method)

  
  permanova_pairwise <- vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps)$research, permutations = 9999)
  betadisper(ps_dist_matrix, sample_data(ps)$research)
  return(c(permanova_pairwise$`Pr(>F)`[1], permanova_pairwise$R2[1]))
}


AnosimCLR <- function(ps, method){
  metadata <- data.frame(sample_data(ps))
  physeq_clr <- microbiome::transform(ps, "clr")
  otu_clr <- as.data.frame(otu_table(physeq_clr))
  dist_matrix <- vegan::vegdist(otu_clr, method = "euclidean")
  anosim_test<-anosim(dist_matrix, 
                             metadata$research, permutations = 9999)
  return(c(anosim_test$signif, anosim_test$statistic))
}
