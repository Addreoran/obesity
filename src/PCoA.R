#library(BiocManager)
#BiocManager::install("microbiome")
library(microbiome)

PCoAAitch <-function(ps, save_path, width=10, height=8){
  #ps <- transform_sample_counts(ps, function(x) x + 1)
  ps.clr <- microbiome::transform(ps, "clr")
  pcoa_euclidean <- ordinate(ps.clr, "PCoA", "euclidean")
  image <- phyloseq::plot_ordination(ps.clr, pcoa_euclidean, color = "research") +
    geom_point(size = 3) +
    geom_text(aes(label = sample_data(ps.clr)$CAP), vjust = -0.5) +
    geom_polygon(stat = "ellipse", aes(fill = research), alpha = 0.3)
  ggsave(file=save_path, plot=image, width=10, height=8)
}

PCoABray <- function(ps, save_path, width=10, height=8){
  ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
  ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
  
  image<-plot_ordination(ps.prop, ord.nmds.bray, color = "research") +
    geom_point(size = 3)+
    geom_text(aes(label = sample_data(ps.prop)$CAP), vjust = -0.5)  +
    geom_polygon(stat = "ellipse", aes(fill = research), alpha = 0.3)
  ggsave(file=save_path, plot=image, width=10, height=8)
}
