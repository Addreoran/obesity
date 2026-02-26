# alfa-diversity



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




