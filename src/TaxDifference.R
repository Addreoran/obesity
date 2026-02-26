

tax_difference <- function(ps, folder, suffix, tax_lvl="genus"){
  ps_differential<-tax_glom(ps, tax_lvl)
  ps.rel = transform_sample_counts(ps_differential, function(x) x/sum(x)*100)
  glom <- tax_glom(ps.rel, taxrank = tax_lvl, NArm = FALSE)
  ps.melt <- psmelt(glom)

  ps.melt$Tax <- as.character(ps.melt[,tax_lvl])
  ps.melt$Tax <- gsub(".*__", "", ps.melt$Tax)
  ps.melt <- ps.melt %>%
    group_by(research, Tax) %>%
    mutate(median=median(Abundance))
  keep <- unique(ps.melt$Tax[ps.melt$median > 1])
  ps.melt$Tax[!(ps.melt$Tax %in% keep)] <- "< 1%"
  
  ps.melt_sum <- ps.melt %>%
    group_by(research,Tax) %>%
    summarise(Abundance=sum(Abundance))
  
  g <- ggplot(ps.melt_sum, aes(x = research, y = Abundance, fill = Tax)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x="", y="Relative abundance (%)") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90))
  ggsave(
    paste0(folder, suffix, "_", tax_lvl, "_tax_diff.png"),
    plot = g,
  )
}
