#save.image('data/data.RData')
#load('data/data.RData')
server = function(input, output, session) {
  ### Barplot ###
  dataset_input = reactive({
    switch(paste(input$dataset, input$normalization_method),
           'GCA_000006275.2 (JCVI-afl1-v2.0) TPM (transcripts per million)' = tpm_jcvi,
           'GCA_009017415.1 (chromosome level) TPM (transcripts per million)' = tpm_chrom_level,
           'GCA_000006275.2 (JCVI-afl1-v2.0) VST (variance stabilizing transformation)' = vst_jcvi,
           'GCA_009017415.1 (chromosome level) VST (variance stabilizing transformation)' = vst_chrom_level)
  })
  dataset_pca = reactive({
    switch(input$dataset_pca,
           'GCA_000006275.2 (JCVI-afl1-v2.0)' = vst_jcvi,
           'GCA_009017415.1 (chromosome level)' = vst_chrom_level)
  })
  rv = reactiveValues()
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
  ### Multi-gene heatmap ###
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
  ## Create metadata table when a sample in the PCA plot is clicked
  output$sample_metadata_table_pca = renderDataTable({
    click_data = event_data('plotly_click')
    if (is.null(click_data)) 'Click on a sample dot to get further sample data' else {
      metadata %>% filter(run == click_data$key) %>%
        select(-sample_description) %>%
        mutate(across(everything(), as.character)) %>%
        pivot_longer(everything(), names_to='category', values_to='value') %>%
        filter(!is.na(value))  
    }
  }, escape = FALSE, options = list(paginate=FALSE, info = FALSE, sort=FALSE))
  output$pca = renderPlotly({
    create_pca()
  }) 
  ### Co-expression network ###
  create_network = function(gene_ids, network){
    #gene_ids = str_split_1(gene_ids, ',|\\s')
    gene_ids = gene_ids[gene_ids %in% V(network)$name]
    if (input$show_neighbors){
      node_neighbors = ego(network, order = 1, nodes = gene_ids)
      g = induced_subgraph(network, unlist(node_neighbors))
      V(g)$selected = V(g)$name  %in% gene_ids
      V(g)$group = ifelse(V(g)$name %in% gene_ids, 'query_genes', 'neighbors')
      data = toVisNetworkData(g)
      visNetwork(nodes = data$nodes, edges = data$edges) %>%
        visIgraphLayout(randomSeed=1234, type = "full") %>%
        visEdges(width=1, color='lightgray', smooth = FALSE) %>%
        visNodes(opacity=0.5) %>%
        visGroups(groupname = "neighbors", color = "lightgrey") %>%
        visGroups(groupname = "query_genes", color = "lightblue", size=30) %>%
        visLegend(position = 'right')
    } else { # Since not coloring which genes are original vs neighbor, can color by annotation
      g = induced_subgraph(network, gene_ids)
      data = toVisNetworkData(g)
      if (input$annotation_category_network != 'Gene list (Comma separated)'){
        annotations = input$gene_categories_network %>% str_remove(' \\(\\d+? genes\\)')
        temp = data$nodes %>% select(id, group = !!sym(input$annotation_category_network)) %>% 
          separate_longer_delim(group, ';') %>% 
          filter(group %in% annotations) %>% 
          group_by(id) %>%
          summarize(group = str_c(sort(group), collapse=';'))
        data$nodes = data$nodes %>% left_join(temp, by='id')
        visnet = visNetwork(nodes = data$nodes, edges = data$edges)
        #for (idx in seq(1:length(annotations))){
        #  visnet = visGroups(visnet, groupname = annotations[idx], color = mpn65[idx])
        #}
        visnet %>%
          visIgraphLayout(randomSeed=1234, type = "full") %>%
          visEdges(width=1, color='lightgray', smooth = FALSE) %>%
          visNodes(opacity=0.5) %>%
          visLegend(position = 'right')
      } else { # Gene list provided
      visNetwork(nodes = data$nodes, edges = data$edges) %>%
        visIgraphLayout(randomSeed=1234, type = "full") %>%
        visEdges(width=1, color='lightgray', smooth = FALSE) %>%
        visNodes(opacity=0.5) %>%
        visGroups(groupname = "neighbors", color = "lightgrey") %>%
        visGroups(groupname = "query_genes", color = "lightblue", size=30) %>%
        visLegend(position = 'right')
      }
    }
  }
  observeEvent(input$generate_network,{
      output$network_vis = renderVisNetwork({
        create_network(annotation_list[[input$annotation_category_network]] %>% 
                         filter(display_text %in% input$gene_categories_network) %>% 
                         pull(gene_id) %>% unlist() %>%
                         str_split(',') %>% unlist(), 
                       network) %>%
          visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
        })
      })
  observe({
    updateSelectInput(session, "gene_categories_network",
                      choices = annotation_list[[input$annotation_category_network]] %>% 
                        pull(display_text))
  })
  myNode = reactiveValues(selected = '')
  observeEvent(input$current_node_id, {
    myNode$selected <<- input$current_node_id
  })
  output$nodes_data_from_shiny  = function(){
    gene_data = functional_annotation_jcvi %>%
      filter(gene_id == myNode$selected) %>%
      select(all_of(c('gene_id', annotation_categories[1:5]))) %>%
      pivot_longer(everything(), names_to='category', values_to='value') %>%
      filter(!is.na(value)) 
    kable(gene_data, 'html', col.names = NULL) %>%
      kable_styling(bootstrap_options = c("striped")) %>%
      pack_rows('gene data',  1, nrow(gene_data)) %>%
      HTML()
  }
  
  output$dt_UI = renderUI({
    if (length(myNode$selected) > 0) {
    #if (TRUE) {
      tableOutput('nodes_data_from_shiny')
    } else{
    }
  })

  ### PCA ###
  create_pca = function(){
    vsd = dataset_pca()
    gene_variance = matrixStats::rowVars(vsd %>% column_to_rownames('gene_id') %>% as.matrix())
    # select the top 4,000 genes by variance
    most_variable_genes = order(gene_variance, decreasing=TRUE)[1:4000]
    pca = prcomp(t((vsd %>% column_to_rownames('gene_id') %>% as.matrix())))
    # the contribution to the total variance for each component
    percent_var = pca$sdev^2 / sum(pca$sdev^2) 
    percent_var = round(100 * percent_var, 2)
    # assemble data for biplot
    pca_scores = pca$x[,1:10] %>% as_tibble() %>% 
      bind_cols(run = names(vsd)[2:length(vsd)]) %>%
      left_join(metadata, by='run')
    pca_scores %>%
      mutate(across(c(strain, sample_type,  source_name, genotype), 
                    ~str_remove(.x, regex('^A.* flavus ', ignore_case = TRUE)))) %>%
      mutate(across(c(strain, sample_type,  source_name, genotype, treatment, isolate), 
                    ~na_if(.x, 'missing'))) %>%
      mutate(strain = str_replace(strain, 'NRRL3357|NR3357', 'NRRL 3357')) %>%
      unite('sample_description', strain, sample_type,  source_name, genotype, treatment, isolate,
            sep = ';', na.rm =TRUE, remove = FALSE) %>% 
      ggplot(aes(.data[[paste0('PC',input$pc_x)]], .data[[paste0('PC',as.numeric(input$pc_x) + 1)]], 
                 color=.data[[input$category_to_color_pca]], label=sample_description)) +
      geom_point(alpha=0.5) +
      xlab(paste0("PC", input$pc_x, ": ",percent_var[as.numeric(input$pc_x)],"% variance")) +
      ylab(paste0("PC", as.numeric(input$pc_x) + 1, ": ",percent_var[as.numeric(input$pc_x) + 1],"% variance")) +
      theme_bw() +
      scale_fill_manual(values = mpn65) +
      ggeasy::easy_remove_legend()
  }
}