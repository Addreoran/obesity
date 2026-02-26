#BiocManager::install('DESeq2')
#BiocManager::install('ANCOMBC')

library(DESeq2)
library(ANCOMBC)

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
