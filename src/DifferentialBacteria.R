#BiocManager::install('DESeq2')
#BiocManager::install('ANCOMBC')

library(DESeq2)
library(ANCOMBC)
library(ggplot2)
library(tidyr)

deseq2_result<-function(ps){
  no_tries <- rowSums(otu_table(ps)>0)
  number_of_probe <- nrow(sample_data(ps))
  tries <- names(no_tries)[no_tries>0.2*number_of_probe] #58
  otu_table(ps) <- otu_table(ps)[tries,]
  
  diagdds = phyloseq_to_deseq2(ps, ~research)
  diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
  res = results(diagdds, cooksCutoff = FALSE)
  res<-as(res, "data.frame")
  return(res)
}


ancombc2_result <- function(ps){
  output_ancomb2_ark5 <- ancombc2(data = ps, 
                                p_adj_method = "BH",
                               fix_formula = "research",
                               prv_cut=0.2,
                               pseudo_sens =FALSE)
  result <- output_ancomb2_ark5$res
  return(result)
}

 # tax_tab<-as.data.frame(unclass(tax_table(ps)))
 # otu_tab<-as.data.frame(unclass(otu_table(ps)))
deseq2_ancombc2_result <- function(deseq_res, ancomb_res, tax_table, otu_table, save_path){
  deseq2_result<-deseq_res
  ancombc2_result<-ancomb_res
  deseq2_result$taxon <- rownames(deseq2_result)  
  df_to_save <- merge(deseq2_result,ancombc2_result, by.x=0,by.y="taxon",all = T )
  df_additional_info<-merge(tax_table, otu_table, by=0, all = T )
  write.csv2(merge(df_to_save, df_additional_info, by='Row.names', all = T ), save_path)
}


plot_ancomb2_significant_results<-function(ancombc2, tax_tab, save_path){

df_ancomb <- merge(ancombc2, tax_tab,by.x="taxon" ,by.y=0, all = T)
q_name <- colnames(df_ancomb)[grepl("q_research", colnames(df_ancomb))]
lfc_name <- colnames(df_ancomb)[grepl("lfc_research", colnames(df_ancomb))]

genus_name <- colnames(df_ancomb)[grepl("genus", colnames(df_ancomb))]
df_ancomb <- df_ancomb[df_ancomb[,q_name] <= 0.05,]
df <- data.frame(
  group = paste(df_ancomb[,"family"], df_ancomb[,genus_name], sep=";"),
  left  =  ifelse(df_ancomb[,lfc_name] < 0, df_ancomb[,lfc_name], 0),
  right =  ifelse(df_ancomb[,lfc_name] >= 0, df_ancomb[,lfc_name], 0)
)
  df <- df[!is.na(df$left) & !is.na(df$right), ]
df_long <- pivot_longer(df, cols = c(left, right),
                        names_to = "side",
                        values_to = "value")
  

df_long$group <- reorder(df_long$group, df_long$value)
df_long <- df_long |>
  dplyr::group_by(side) |>
  dplyr::mutate(group = reorder(group, value)) |>
  dplyr::ungroup()
p <- ggplot(df_long, aes(x = value, y = group)) +
  geom_col(fill = "lightblue") +
  
  facet_grid(. ~ side, scales = "free_x", space = "free_x") +
  
    theme_minimal() +
    theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(3, "cm")  # więcej miejsca na etykiety
  ) + geom_text(aes(label=round(value,2))) +
  
  geom_text(data = subset(df_long, side == "left"),
            aes(x = 0, label = gsub(";", "\n", group)),
            hjust = 0.5,
            nudge_x = 1,
            
            size = 3)
  ggsave(
          paste0(save_path),
          plot = p,
  )
}
