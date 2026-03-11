

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


tax_stats <- function(ps, tax_lvl, folder){
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


stats_to_publication <- function(ps, tax_lvl, control, researched, folder){
    ps_genus <- tax_glom(ps, "genus")
  
    sum(otu_table(ps_genus))/nrow(sample_data(ps_genus))
    
    ps_phylum <- tax_glom(ps_genus, "phylum")
    ps_prop <- transform_sample_counts(physeq = ps_phylum, fun = function(x){x/sum(x)})
    keep <- filter_taxa(ps_prop, function(x){mean(x)  > 0.01})
    # there is 1%
    ps_filtered_phylum <- prune_taxa(taxa = keep, x = ps_phylum)
    
    ps_filtered_phylum
    ps_prop <- transform_sample_counts(physeq = ps_genus, fun = function(x){x/sum(x)})
    keep <- filter_taxa(ps_prop, function(x){mean(x) > 0.0001})
    # there is 0.01%
    ps_filtered <- prune_taxa(taxa = keep, x = ps_genus)
                        
    
    subset_samples(ps_genus, research == control) -> ps_kontrola
    keep <- filter_taxa(ps_kontrola, function(x){mean(x) > 0})
    ps_kontrola_filtered <- prune_taxa(taxa = keep, x = ps_kontrola)
    
    
                                        
    ps_kontrola_filtered_phylum<-tax_glom(ps_kontrola, taxrank = 'phylum')
    
    ps_kontrola_filtered_phylum <- phyloseq::transform_sample_counts(ps_kontrola_filtered_phylum, function(x) { x/sum(x) } )
    phyla_rel_bact <- otu_table(subset_taxa(ps_kontrola_filtered_phylum, phylum == "p__Bacteroidota"))
    phyla_rel_firm <- suppressWarnings(otu_table(subset_taxa(ps_kontrola_filtered_phylum, phylum == "p__Firmicutes")))
    bf_rati_kontrola <- phyla_rel_bact /  phyla_rel_firm
    mean(bf_rati_kontrola)
    
    subset_samples(ps_genus, research == researched') -> ps_other
    keep <- filter_taxa(ps_other, function(x){mean(x) > 0})
    ps_other_filtered <- prune_taxa(taxa = keep, x = ps_other)
    sum(tax_table(ps_other_filtered)[,'genus'] %in% tax_table(ps_kontrola_filtered)[,'genus'])
    ps_other_filtered_phylum<- tax_glom(ps_other, taxrank = 'phylum')
                                        
    ps_other_filtered_phylum <- phyloseq::transform_sample_counts(ps_other_filtered_phylum, function(x) { x/sum(x) } )
    phyla_rel_bact <- otu_table(subset_taxa(ps_other_filtered_phylum, phylum == "p__Bacteroidota"))
    phyla_rel_firm <- suppressWarnings(otu_table(subset_taxa(ps_other_filtered_phylum, phylum == "p__Firmicutes")))
    bf_ratio_other <- phyla_rel_bact /  phyla_rel_firm
    mean(bf_ratio_other)
    wilcox.test(bf_ratio_other[1,], bf_rati_kontrola[1,], 
                                       p.adjust.method = "BH")
  }
        
