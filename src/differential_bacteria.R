#BiocManager::install('DESeq2')
#BiocManager::install('ANCOMBC')

library(DESeq2)
library(ANCOMBC)

deseq2_result<-function(ps, ){
  no_tries=rowSums(otu_table(ps)>0)
  names(no_tries) # sprawdzić czy dobrze liczone
  tries=names(no_tries)[no_tries>0.2*no_tries] #58
  otu_table(ps)<-otu_table(ps)[tries,]
  
  
  diagdds = phyloseq_to_deseq2(ps, ~research)
  diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
  res = results(diagdds, cooksCutoff = FALSE)
  alpha = 1
  sigtab = cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  return(sigtab)
}

ancombc2_result <- function(ps, ps_tax_level){
  output_ancomb2_ark5 <- ancombc2(data = ps, tax_level = ps_tax_level,
                                p_adj_method = "BH",
                               fix_formula = "research",
                               prv_cut=0.2,
                               pseudo_sens =FALSE)
  result <- output_ancomb2_ark5$res
  tmp_tax <- data.frame(tax_table(ps))
  tmp_tax$otu_name <- rownames(tax_table(ps))
  tmp_otu <- data.frame(otu_table(ps))
  tmp_otu$otu_name <- rownames(otu_table(ps))
  tmp_otu_tax <- merge(tmp_tax, tmp_otu, by="otu_name")
  ancom_new_res <- merge(tmp_otu_tax,result, by.x="genus", by.y='taxon',all = T )
  return(ancom_new_res)
}


deseq2_ancombc2_result <- function(deseq2_result, ancombc2_result, save_path){
    write.csv2(merge(deseq2_result,ancombc2_result, by="genus",all = T ), save_path)

}
