# alfa-diversity
library(dplyr)
library(tibble)

AlfaDiversity <- function(ps, folder, suffix){

  p <- plot_richness(ps_genus, "research", measures = c("Chao1",  "Observed"))
  p <- p + theme_bw() + geom_boxplot(aes(fill = research))+ theme(axis.text.x = element_text(angle=90, hjust=1), legend.position="none") 
  
  ggsave(
    paste0(folder, suffix, "entropy_1.svg"),
    plot = p,
    width=8, height=8
  )
    p <- plot_richness(ps_genus, "research", measures = c("Shannon", "Simpson", "Fisher"))
  p <- p + theme_bw() + geom_boxplot(aes(fill = research))+ theme(axis.text.x = element_text(angle=90, hjust=1), legend.position="none") 
  
  ggsave(
    paste0(folder, suffix, "entropy_2.svg"),
    plot = p,
    width=12, height=8
  )
  
}


TestDiversity <- function(ps, folder, suffix){
  richness<-estimate_richness(ps)
  richness$id<-sub("X", "", rownames(richness))
  
  df_sample<-as.data.frame(as.matrix(sample_data(ps)))
  df_sample$id<-sub(" ", "", df_sample$id)
  df_sample$id<-sub(" ", "", df_sample$id)
  tmp <- merge(df_sample, richness, by="id")

  res_chao<-kruskal.test(Chao1 ~ research, data = tmp)$p.value
  res_shan<-kruskal.test(Shannon ~ research, data = tmp)$p.value
  res_simpson<-kruskal.test(Simpson ~ research, data = tmp)$p.value
  res_fisher<-kruskal.test(Fisher ~ research, data = tmp)$p.value
  res_observed<-kruskal.test(Observed ~ research, data = tmp)$p.value


  res_chao_wilcoxon <- pairwise.wilcox.test(tmp$Chao1, tmp$research, 
                               p.adjust.method = "BH")$p.value[1]
  res_shan_wilcoxon <- pairwise.wilcox.test(tmp$Shannon, tmp$research, 
                               p.adjust.method = "BH")$p.value[1]
  res_simpson_wilcoxon <- pairwise.wilcox.test(tmp$Simpson, tmp$research, 
                                  p.adjust.method = "BH")$p.value[1]
  res_fisher_wilcoxon <- pairwise.wilcox.test(tmp$Fisher, tmp$research, 
                                 p.adjust.method = "BH")$p.value[1]
  res_observed_wilcoxon <- pairwise.wilcox.test(tmp$Observed, tmp$research, 
                                   p.adjust.method = "BH")$p.value[1]

  result <- data.frame(wilcoxon=c(res_chao_wilcoxon, res_shan_wilcoxon, res_simpson_wilcoxon, res_fisher_wilcoxon, res_observed_wilcoxon), 
                     kruskal=c(res_chao, res_shan, res_simpson, res_fisher, res_observed), row.names=c("Chao1", "Shannon", "Simpson", "Fisher", "Observed"))
  write.csv(result, paste0(folder, suffix, "_wilcoxon_kruskal.csv"))

  
}
