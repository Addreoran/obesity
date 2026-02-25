library(vegan)
library(permute)



ps_merged_glom_ark5<-tax_glom(ps, "genus")


###
metadata_sub <- data.frame(sample_data(ps_merged_glom_ark5))
permanova_pairwise <- adonis2(phyloseq::distance(ps_merged_glom_ark5, method = "euclidean") ~ research, 
                              data = metadata_sub)
p <- c(p, permanova_pairwise$`Pr(>F)`[1])


metadata_sub <- data.frame(sample_data(ps_merged_glom_ark5))
permanova_anosim<- anosim(phyloseq::distance(ps_merged_glom_ark5, method="euclidean", binary = TRUE), 
                             metadata_sub$research)
p <- c(p, permanova_pairwise$signif)


###
metadata_sub <- data.frame(sample_data(ps_merged_glom_ark5))
permanova_pairwise <- adonis2(phyloseq::distance(ps_merged_glom_ark5, method = "euclidean") ~ research, 
                              data = metadata_sub)
p <- c(p, permanova_pairwise$`Pr(>F)`[1])


metadata_sub <- data.frame(sample_data(ps_merged_glom_ark5))
permanova_anosim<- anosim(phyloseq::distance(ps_merged_glom_ark5, method="euclidean", binary = TRUE), 
                             metadata_sub$research)
p <- c(p, permanova_pairwise$signif)
