

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
  tax_order <- ps.melt_sum %>%
    group_by(Tax) %>%
    summarise(total = sum(Abundance)) %>%
    arrange(desc(total)) %>%
    pull(Tax)

  ps.melt_sum$Tax <- factor(ps.melt_sum$Tax,
                          levels = c("< 1%", setdiff(tax_order, "< 1%")))
  g <- ggplot(ps.melt_sum, aes(x = research, y = Abundance, fill = Tax)) + 
    geom_bar(stat = "identity", position = "fill", colour="white", linewidth=0.2) + 
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


tax_stats <- function(ps, tax_lvl){
  sample_sums(ps)
  ps_phylum <- tax_glom(ps, taxrank = tax_lvl)
  ps_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x) * 100)
  df <- psmelt(ps_rel)
  df_summary <- df %>%
    group_by(research, phylum) %>%
    summarise(
        mean_abundance = mean(Abundance),
        median_abundance = median(Abundance)
      )
  write.csv(df_summary, paste0(folder, suffix,"_", tax_lvl,"_percentage.csv"))
  write.csv(sample_sums(ps), paste0(folder, suffix,"_", tax_lvl,"_reads_per_probe.csv"))                                   
}
