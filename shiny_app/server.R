#save.image('data/data.RData')
#load('data/data.RData')
server = function(input, output, session) {
  dataset_input = reactive({
    switch(paste(input$dataset, input$normalization_method),
           'GCA_000006275.2 (JCVI-afl1-v2.0) TPM (transcripts per million)' = tpm_jcvi,
           'GCA_009017415.1 (chromosome level) TPM (transcripts per million)' = tpm_chrom_level,
           'GCA_000006275.2 (JCVI-afl1-v2.0) VST (variance stabilizing transformation)' = vst_jcvi,
           'GCA_009017415.1 (chromosome level) VST (variance stabilizing transformation)' = vst_chrom_level)
  })
  rv <- reactiveValues()
  single_gene_barplot = function(df, gene_of_interest){
    dataset_input() %>%
      filter(gene_id == gene_of_interest) %>%
      pivot_longer(-gene_id, names_to='sra_run', values_to='TPM') %>%
      left_join(metadata, by=c('sra_run' = 'run')) %>%
      arrange(bioproject, sra_run) %>%
      mutate(bioproject = fct_inorder(bioproject),
             sra_run = fct_inorder(sra_run)) %>%
      ggplot(aes(sra_run, TPM, fill=bioproject, label=sample_description, key=sra_run)) +
      geom_col() +
      scale_fill_manual(values = mpn65) +
      labs(title = gene_of_interest, x='Sample') +
      ggeasy::easy_remove_legend() +
      ggeasy::easy_remove_x_axis(what=c('tics', 'text', 'line')) +
      ggeasy::easy_center_title()
  }
  
  output$barplot = renderPlotly({
    single_gene_barplot(dataset_input(), input$gene_id)
  }) 
  ## Change default gene_id value if the dataset changes
  observe({
    updateTextInput(inputId = 'gene_id', 
                    value = ifelse(input$dataset == 'GCA_000006275.2 (JCVI-afl1-v2.0)', 
                                   'AFLA_139360', 'F9C07_7811'))
    rv$normalization_method_short = str_remove(input$normalization_method, ' .*')
  })
  ## Create metadata table when bar in plot is clicked
  output$sample_metadata_table = renderDataTable({
    click_data = event_data('plotly_click')
    if (is.null(click_data)) 'Click on a bar to get further sample data' else {
      metadata %>% filter(run == click_data$key) %>%
        select(-sample_description) %>%
        mutate(across(everything(), as.character)) %>%
        pivot_longer(everything(), names_to='category', values_to='value') %>%
        filter(!is.na(value))  
    }
  }, escape = FALSE, options = list(paginate=FALSE, info = FALSE, sort=FALSE))
  ## Add a download button to download a csv of the expression data
  output$download_single_gene_data = downloadHandler(
    filename = function() {
      paste0(input$gene_id, '_data.csv')
    },
    content = function(file) {
      df = dataset_input()
      df = df %>% filter(gene_id == input$gene_id) %>%
        pivot_longer(-gene_id, names_to='sra_run', values_to='TPM') %>%
        left_join(metadata, by=c('sra_run' = 'run')) %>%
        arrange(bioproject, sra_run)
      write.csv(df, file, row.names = FALSE)
    })
  ### Multi-gene heatmap 
  output$download_multi_gene_data = downloadHandler(
    filename = function(){
      paste0('A_flavus_', rv$normalization_method_short, '.csv')
      },
    content = function(file){
      df = dataset_input()
      sra_runs_to_include = metadata %>%
        filter(bioproject %in% input$bioprojects_to_include) %>%
        pull(run)
      gene_ids_to_include =  annotation_list[[input$annotation_category]] %>% 
        filter(display_text %in% input$gene_categories) %>% 
        pull(gene_id) %>% unlist()
      df = df %>% 
        filter(gene_id %in% gene_ids_to_include) %>%
        select(all_of(c('gene_id', sra_runs_to_include)))
      write.csv(df, file, row.names = FALSE)
    }
  )
  multi_gene_heatmap = function(df, genes_of_interest, bioprojects_to_include){
    #genes_of_interest = str_split_1(genes_of_interest, ' +|,+')
    df = df %>%
      filter(gene_id %in% genes_of_interest) %>%
      pivot_longer(-gene_id, names_to='sra_run', values_to='expression_value') %>%
      left_join(metadata, by=c('sra_run' = 'run')) %>%
      filter(bioproject %in% bioprojects_to_include) 
    if (rv$normalization_method_short == 'TPM'){
      df = df %>% mutate(expression_value = log1p(expression_value))
    } 
    if (length(input$gene_categories) > 1) {
      df = df %>% 
        left_join(
          annotation_list[[input$annotation_category]] %>% 
            filter(display_text %in% input$gene_categories) %>%
            unnest(gene_id),
          by='gene_id') %>%
        data.table::setnames('column_to_select', input$annotation_category)
    }
    if (length(bioprojects_to_include) > 1) {
      ht = df %>%
        group_by(bioproject) %>%
        heatmap(gene_id, sra_run, expression_value, column_names_gp = gpar(fontsize = 8),
                row_names_gp = gpar(fontsize = 8),
                heatmap_legend_param = list(title = str_replace(rv$normalization_method_short, 'TPM', 'log_TPM')))
    } else{
      ht = df %>%
        heatmap(gene_id, sra_run, expression_value, column_names_gp = gpar(fontsize = 8),
                row_names_gp = gpar(fontsize = 8),
                heatmap_legend_param = list(title = str_replace(rv$normalization_method_short, 'TPM', 'log_TPM')))
    }
    if (length(input$gene_categories) > 1 & input$annotation_category != 'Gene list (Comma separated)') {
      column = sym(input$annotation_category)
      ht = ht %>% annotation_tile(!!column, palette = mpn65)
    }
    ht = ht %>% as_ComplexHeatmap()
    rv$m = ht@matrix
    ht = draw(ht)
    return(ht)
  }
  observeEvent(input$generate_heatmap,{
    makeInteractiveComplexHeatmap(
      input, output, session,
      ht_list = multi_gene_heatmap(dataset_input(), 
                                   annotation_list[[input$annotation_category]] %>% 
                                     filter(display_text %in% input$gene_categories) %>% 
                                     pull(gene_id) %>% unlist() %>%
                                     str_split(',') %>% unlist(), 
                                   input$bioprojects_to_include),
      click_action = click_action_heatmap
    )
  })
  observe({
    updateSelectInput(session, "gene_categories",
                      choices = annotation_list[[input$annotation_category]] %>% 
                        pull(display_text))
  })
  click_action_heatmap = function(df, output) {
    output[["gene_info"]] = renderUI({
      if(!is.null(df)) {
        row_index = unique(unlist(df$row_index))
        column_index = unique(unlist(df$column_index))
        gene = rownames(rv$m)[row_index]
        sra_run = colnames(rv$m)[column_index]
        log_TPM = rv$m[row_index, column_index, drop = TRUE]
        run_data = metadata %>% 
          filter(run == sra_run) %>%
          select(-sample_description) %>%
          mutate(across(everything(), as.character)) %>%
          pivot_longer(everything(), names_to='category', values_to='value') %>%
          filter(!is.na(value))  
        gene_data = functional_annotation_jcvi %>%
          filter(gene_id == gene) %>%
          select(all_of(c('gene_id', annotation_categories[1:5]))) %>%
          pivot_longer(everything(), names_to='category', values_to='value') %>%
          filter(!is.na(value)) %>%
          add_row(category = 'log_TPM', value = as.character(round(log_TPM, 2)))
        kable(bind_rows(gene_data, run_data), 'html', col.names = NULL) %>%
          kable_styling(bootstrap_options = c("striped")) %>%
          pack_rows('gene data',  1, nrow(gene_data)) %>%
          pack_rows('SRA run data',  nrow(gene_data) + 1, nrow(run_data)) %>%
          HTML()
      }
    })
  }
}