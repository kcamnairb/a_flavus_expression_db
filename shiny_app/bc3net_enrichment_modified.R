enrichment_test = function(gene_list, functional_annotation, columns_list){
  columns_list = columns_list %>% purrr::set_names()
  named_group_split = function(...) {
    data <- group_by(...)
    names <- group_keys(data) %>% 
      map(as.character) %>% 
      reduce(paste, sep = "~~")
    group_split(data) %>% 
      set_names(names)
  }
  enrichment_bc3net = function (genes, reference, genesets){
    tab = lapply(1:length(genesets), function(i) {
      reference = reference[!reference %in% genes]
      RinSet = sum(reference %in% genesets[[i]])
      RninSet = length(reference) - RinSet
      GinSet = sum(genes %in% genesets[[i]])
      GninSet = length(genes) - GinSet
      fmat = matrix(c(GinSet, RinSet, GninSet, RninSet), nrow = 2, 
                    ncol = 2, byrow = F)
      colnames(fmat) = c("inSet", "ninSet")
      rownames(fmat) = c("genes", "reference")
      fish = fisher.test(fmat, alternative = "greater")
      pval = fish$p.value
      inSet = RinSet + GinSet
      genes_in_cat = genes[genes %in% genesets[[i]]] %>% str_flatten(collapse=';')
      #genes_in_cat = genes_in_cat[genes_in_cat != ''] %>% str_flatten(collapse=';')
      res = c(GinSet, inSet, pval, genes_in_cat)
      res
    })
    rtab = do.call("rbind", tab)
    rtab = data.frame(as.vector(names(genesets)), rtab)
    rtab = rtab[order(rtab[, 4]), ]
    #colnames(rtab) = c("TermID", "genes", "all", "pval")
    colnames(rtab) = c("TermID", "genes", "all", "pval", "gene_ids") 
    tab.out = data.frame(rtab)
    return(tab.out)
  }
  enrichment_one_col = function(column, genes, annotation){
    cat2gene_lists = annotation %>% 
      select(gene_id, column_to_select = {{ column }}) %>%
      separate_rows(column_to_select, sep=";") %>%
      drop_na(column_to_select) %>%
      named_group_split(column_to_select) %>%
      map(~ .x$gene_id) %>% keep(~length(.x) > 3)
    enrichment_bc3net(genes, annotation$gene_id, genesets=cat2gene_lists) %>%
      mutate(annotation_category = {{column}})
  }
  
  run_gene_list_with_columns = function(gene_list, gene_list_name, column_name){
    enrichment_one_col(column_name, gene_list, functional_annotation) %>%
      mutate(gene_list_name = {{gene_list_name}})
  }
  combine_gene_list_with_columns = function(gene_list, gene_list_name){
    imap(columns_list, ~ run_gene_list_with_columns(gene_list, gene_list_name, .y))
  }
    imap(gene_list, ~combine_gene_list_with_columns(.x, .y)) %>%
    flatten_dfr()
}
