server = function(input, output, session) {
  rv = reactiveValues()
  ### Barplot ###
  dataset_input_barplot = reactive({
    switch(paste(input$normalization_method_barplot, input$dataset_barplot),
           'TPM JCVI' = read_fst('data/A_flavus_jcvi_tpm.fst') %>% as_tibble(),
           'TPM chrom_level' = read_fst('data/A_flavus_chrom_level_tpm.fst') %>% as_tibble(),
           'VST JCVI' = read_fst('data/vst_jcvi.fst') %>% as_tibble(),
           'VST chrom_level' = read_fst('data/vst_chrom_level.fst') %>% as_tibble())
  })
  single_gene_barplot = function(df, gene_of_interest){
    dataset_input_barplot() %>%
      filter(gene_id == gene_of_interest) %>%
      pivot_longer(-gene_id, names_to='sra_run', values_to='expr_value') %>%
      left_join(metadata, by=c('sra_run' = 'run')) %>%
      arrange(bioproject, sra_run) %>%
      mutate(bioproject = fct_inorder(bioproject),
             sra_run = fct_inorder(sra_run)) %>%
      ggplot(aes(sra_run, expr_value, fill=bioproject, label=sample_description, key=sra_run)) +
      geom_col() +
      scale_fill_manual(values = mpn65) +
      labs(title = gene_of_interest, x='Sample', y=input$normalization_method_barplot) +
      ggeasy::easy_remove_legend() +
      ggeasy::easy_remove_x_axis(what=c('tics', 'text', 'line')) +
      ggeasy::easy_center_title()
  }
  output$barplot = renderPlotly({
    single_gene_barplot(dataset_input_barplot(), input$gene_id)
  }) 
  # Change default gene_id value if the dataset changes
  observe({
    updateTextInput(inputId = 'gene_id', 
                    value = ifelse(input$dataset_barplot == 'JCVI', 'AFLA_139360', 'F9C07_7811'))
  })
  # Create metadata table when bar in plot is clicked
  output$sample_metadata_table = DT::renderDT({
    click_data = event_data('plotly_click')
    validate(need(!is.null(click_data), 'Click on a bar to get further sample data'))
      metadata %>% filter(run == click_data$key) %>%
        select(-sample_description) %>%
        mutate(across(everything(), as.character)) %>%
        pivot_longer(everything(), names_to='category', values_to='value') %>%
        filter(!is.na(value))  
  }, escape = FALSE, options = list(paginate=FALSE, info = FALSE, sort=FALSE, dom = 't'), rownames = FALSE)
  ## Add a download button to download a table of the expression data
  output$download_single_gene_data = downloadHandler(
    filename = function() {
      paste0(input$gene_id, '_data.xlsx')
    },
    content = function(file) {
      df = dataset_input_barplot()
      df = df %>% filter(gene_id == input$gene_id) %>%
        pivot_longer(-gene_id, names_to='sra_run', values_to='TPM') %>%
        left_join(metadata, by=c('sra_run' = 'run')) %>%
        arrange(bioproject, sra_run)
      openxlsx::write.xlsx(df, file, firstRow = TRUE)
    })
  ### Multi-gene heatmap ###
  dataset_input_heatmap = reactive({
    switch(input$dataset_heatmap,
           'JCVI' = read_fst('data/vst_jcvi.fst') %>% as_tibble(),
           'chrom_level' = read_fst('data/vst_chrom_level.fst') %>% as_tibble())
  })
  output$download_multi_gene_data = downloadHandler(
    filename = 'A_flavus_VST.xlsx',
    content = function(file){
      df = dataset_input_heatmap()
      sra_runs_to_include = metadata %>%
        filter(bioproject %in% input$bioprojects_to_include) %>%
        pull(run)
      gene_ids_to_include =  annotation_list[[input$annotation_category]] %>% 
        filter(display_text %in% input$gene_categories) %>% 
        pull(gene_id) %>% unlist()
      df = df %>% 
        filter(gene_id %in% gene_ids_to_include) %>%
        select(all_of(c('gene_id', sra_runs_to_include)))
      annotation_data = functional_annotation %>% 
        filter(gene_id %in% gene_ids_to_include)
      openxlsx::write.xlsx(list(df, annotation_data) %>% 
                             setNames(c('VST', 'functional_annotation')), 
                           file, firstRow = TRUE)
    }
  )
  output$download_entire_expression_dataset = downloadHandler(
    filename = 'A_flavus_VST_entire_dataset.xlsx',
    content = function(file){
      df = dataset_input_heatmap()
      openxlsx::write.xlsx(list(df, functional_annotation %>% filter(genome == input$dataset_heatmap)) %>% 
                             setNames(c('VST', 'functional_annotation')),
                           file = file,
                           firstRow = TRUE)
    }
  )
  multi_gene_heatmap = function(df, genes_of_interest, bioprojects_to_include){
    df = df %>%
      filter(gene_id %in% genes_of_interest) %>%
      pivot_longer(-gene_id, names_to='sra_run', values_to='VST') %>%
      left_join(metadata, by=c('sra_run' = 'run')) %>%
      filter(bioproject %in% bioprojects_to_include) 
    if (length(input$gene_categories) > 1) {
      df = df %>% 
        left_join(
          annotation_list[[input$annotation_category]] %>% 
            filter(display_text %in% input$gene_categories) %>%
            unnest(gene_id),
          by='gene_id') %>%
        group_by(gene_id, sra_run) %>%
        mutate(column_to_select = str_c(column_to_select, collapse = ';')) %>%
        ungroup() %>%
        data.table::setnames('column_to_select', input$annotation_category)
    }
    if (length(bioprojects_to_include) > 1) {
      ht = df %>%
        group_by(bioproject) %>%
        heatmap(gene_id, sra_run, VST, column_names_gp = gpar(fontsize = 8),
                row_names_gp = gpar(fontsize = 8))
    } else{
      ht = df %>%
        heatmap(gene_id, sra_run, VST, column_names_gp = gpar(fontsize = 8),
                row_names_gp = gpar(fontsize = 8))
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
      ht_list = multi_gene_heatmap(dataset_input_heatmap(), 
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
                        filter(genome == input$dataset_heatmap) %>%
                        pull(display_text))
  })
  click_action_heatmap = function(df, output) {
    output[["gene_info"]] = renderUI({
      if(!is.null(df)) {
        row_index = unique(unlist(df$row_index))
        column_index = unique(unlist(df$column_index))
        gene = rownames(rv$m)[row_index]
        sra_run = colnames(rv$m)[column_index]
        VST = rv$m[row_index, column_index, drop = TRUE]
        run_data = metadata %>% 
          filter(run == sra_run) %>%
          select(-sample_description) %>%
          mutate(across(everything(), as.character)) %>%
          pivot_longer(everything(), names_to='category', values_to='value') %>%
          filter(!is.na(value))  
        gene_data = functional_annotation %>%
          filter(gene_id == gene) %>%
          mutate(`Ensemble Fungi` = paste0('https://fungi.ensembl.org/Aspergillus_flavus/Gene/Summary?g=', gene_id),
                 `Ensemble Fungi` = cell_spec(gene_id, 'html', link = `Ensemble Fungi`, new_tab=TRUE),
                 `Ensemble Fungi` = if_else(genome == 'chrom_level', NA_character_, `Ensemble Fungi`)) %>%
          select(all_of(c('gene_id', 'Ensemble Fungi', 'protein name', annotation_categories[1:5]))) %>%
          pivot_longer(everything(), names_to='category', values_to='value') %>%
          filter(!is.na(value)) %>%
          add_row(category = 'VST', value = as.character(round(VST, 2)))
        kable(bind_rows(gene_data, run_data), 'html', col.names = NULL, escape = FALSE) %>%
          kable_styling(bootstrap_options = c("striped")) %>%
          pack_rows('gene data',  1, nrow(gene_data)) %>%
          pack_rows('SRA run data',  nrow(gene_data) + 1, nrow(run_data)) %>%
          HTML()
      }
    })
  }
 
  ### Co-expression network ###
  add_network_annotation = function(network){
    vertex_attr(network) = functional_annotation %>%
      filter(gene_id %in% V(network)$name) %>%
      mutate(gene_id = factor(gene_id, levels=V(network)$name)) %>%
      arrange(gene_id) %>% 
      rename(name=gene_id) %>%
      full_join(tibble(name = V(network)$name), by='name')
    edge_attr(network, 'abs_weight') = E(network)$weight %>% abs()
    return(network)
  }
  network = reactive({
    switch(input$dataset_network,
           'JCVI' = read_rds('data/get_nondownsampled_data-no_correction-upper_quartile_normalize-no_scaling-pearson.rds') %>%
             add_network_annotation(),
           'chrom_level' = read_rds('data/get_nondownsampled_data_chrom_level-no_correction-upper_quartile_normalize-no_scaling-pearson.rds') %>%
             add_network_annotation()
           )
  })
  create_network = function(gene_ids, network){
    ledges = data.frame(color = c('lightgray', 'hotpink'), 
                        label = c('positive correlation', 'negative correlation'), 
                        font.align = "top") 
    gene_list_provided = input$annotation_category_network == 'Gene list (Comma separated)'
    gene_ids = gene_ids[gene_ids %in% V(network)$name]
    if (input$show_neighbors & gene_list_provided){
      
      ## Showing first order neighboring nodes of provided gene ids
      #gene_ids = 'AFLA_006380'
      node_neighbors = ego(network, order = 1, nodes = gene_ids)
      g = induced_subgraph(network, unlist(node_neighbors))
      ## Limit the number of edges to the top 50000.
      edge_weight_cutoff = E(g)$abs_weight %>% sort(decreasing = TRUE) %>% .[50000]
      ## Creating a new network with the minimum edge weight was a straightforward way to remove isolated nodes.
      network2 = delete.edges(network, which(E(network)$abs_weight < edge_weight_cutoff))
      node_neighbors = ego(network2, order = 1, nodes = gene_ids)
      g = induced_subgraph(network2, unlist(node_neighbors))
      rv$displayed_nodes = unlist(node_neighbors) %>% names() %>% unique()
      V(g)$group = ifelse(V(g)$name %in% gene_ids, 'query_genes', 'neighbors')
      data = toVisNetworkData(g) 
      print(table(data$edges$weight > 0))
      data$edges['color'] = ifelse(data$edges$weight > 0, 'lightgray', 'hotpink')
      rv$network_data = data
      visNetwork(nodes = data$nodes, edges = data$edges) %>%
        visIgraphLayout(randomSeed=1234, type = "full", layout = 'layout_with_fr') %>%
        visEdges(width=3, smooth = FALSE) %>%
        visNodes(opacity=0.5) %>%
        visLegend(position = 'right', main = "Legend", addEdges = ledges)
    } else if (!input$show_neighbors & !gene_list_provided) { 
      
      ## No gene list provided and neighbors not shown, coloring by annotation
      g = induced_subgraph(network, gene_ids)
      rv$displayed_nodes = gene_ids
      data = toVisNetworkData(g)
      data$edges['color'] = ifelse(data$edges$weight > 0, 'lightgray', 'hotpink')
      annotations = input$gene_categories_network %>% str_remove(' \\(\\d+? genes\\)')
      temp = data$nodes %>% select(id, group = !!sym(input$annotation_category_network)) %>% 
        separate_longer_delim(group, ';') %>% 
        filter(group %in% annotations) %>% 
        group_by(id) %>%
        summarize(group = str_c(sort(group), collapse=';'))
      data$nodes = data$nodes %>% left_join(temp, by='id')
      rv$network_data = data
      visNetwork(nodes = data$nodes, edges = data$edges) %>%
        visIgraphLayout(randomSeed=1234, type = "full") %>%
        visEdges(width=3, smooth = FALSE) %>%
        visNodes(opacity=0.5) %>%
        visLegend(position = 'right', main = "Legend", addEdges = ledges)
      } else if (!input$show_neighbors & gene_list_provided) { 
        
      # Gene list provided but don't need neighbors
        validate(need(length(gene_ids) > 1, 'Select at least two genes if "Show neighbors" is not selected'))
        g = induced_subgraph(network, gene_ids)
        rv$displayed_nodes = gene_ids
        data = toVisNetworkData(g)
        validate(need(length(data$edges) > 0, 'No edges between selected genes'))
        rv$network_data = data
      visNetwork(nodes = data$nodes, edges = data$edges) %>%
        visIgraphLayout(randomSeed=1234, type = "full") %>%
        visEdges(width=3, smooth = FALSE) %>%
        visNodes(opacity=0.5) %>%
        visLegend(position = 'right', main = "Legend", useGroups = FALSE, addEdges = ledges)
      } else if (input$show_neighbors & !gene_list_provided) { 
        
        # No gene list provided and show neighbors
        node_neighbors = ego(network, order = 1, nodes = gene_ids)
        g = induced_subgraph(network, unlist(node_neighbors))
        edge_weight_cutoff = E(g)$abs_weight %>% sort(decreasing = TRUE) %>% .[50000]
        network2 = delete.edges(network, which(E(network)$abs_weight < edge_weight_cutoff))
        gene_ids = gene_ids[gene_ids %in% V(network2)$name]
        node_neighbors = ego(network2, order = 1, nodes = gene_ids)
        g = induced_subgraph(network2, unlist(node_neighbors))
        rv$displayed_nodes = unlist(node_neighbors) %>% names() %>% unique()
        data = toVisNetworkData(g)
        data$edges['color'] = ifelse(data$edges$weight > 0, 'lightgray', 'hotpink')
        annotations = input$gene_categories_network %>% str_remove(' \\(\\d+? genes\\)')
        temp = data$nodes %>% select(id, group = !!sym(input$annotation_category_network)) %>% 
          separate_longer_delim(group, ';') %>% 
          filter(group %in% annotations) %>% 
          group_by(id) %>%
          summarize(group = str_c(sort(group), collapse=';'))
        data$nodes = data$nodes %>% left_join(temp, by='id') %>%
          replace_na(list(group = "neighbor"))
        rv$network_data = data
        visNetwork(nodes = data$nodes, edges = data$edges) %>%
          visIgraphLayout(randomSeed=1234, type = "full") %>%
          visEdges(width=3, smooth = FALSE) %>%
          visNodes(opacity=0.5) %>%
          visLegend(position = 'right', main = "Legend", addEdges = ledges)
      }
    }
  observeEvent(input$generate_network,{
      output$network_vis = renderVisNetwork({
        validate(need(!is.null(input$gene_categories_network), message = FALSE))
        create_network(annotation_list_network[[input$annotation_category_network]] %>% 
                         filter(display_text %in% input$gene_categories_network) %>% 
                         pull(gene_id) %>% unlist() %>%
                         str_split(',') %>% unlist(), 
                       network()) %>%
          visEvents(click = "function(nodes){
                Shiny.onInputChange('click', nodes.nodes[0]);
                Shiny.onInputChange('node_selected', nodes.nodes.length);
                ;}"
          )
        })
      })
  observe({
      updateSelectInput(session, "gene_categories_network",
                        choices = annotation_list_network[[input$annotation_category_network]] %>% 
                          filter(genome == input$dataset_network) %>%
                          pull(display_text)
                        )
  })
  output$node_data_from_network  = function(){
    ## https://stackoverflow.com/questions/49913752/undo-click-event-in-visnetwork-in-r-shiny
    if (!is.null(input$node_selected) && (input$node_selected == 1)){
      gene_data = functional_annotation %>%
        filter(gene_id == input$click) %>%
        mutate(across(everything(), as.character)) %>%
        mutate(`Ensemble Fungi` = paste0('https://fungi.ensembl.org/Aspergillus_flavus/Gene/Summary?g=', gene_id),
               `Ensemble Fungi` = cell_spec(gene_id, 'html', link = `Ensemble Fungi`, new_tab=TRUE),
               `Ensemble Fungi` = if_else(genome == 'chrom_level', NA_character_, `Ensemble Fungi`)) %>%
        select(any_of(c('gene_id','protein name',  'KEGG pathways','Gene Ontology','Interpro domains',
                        'Subcellular localization','deeploc_score', 'Ensemble Fungi',
                        'biosynthetic gene clusters','mibig_accession','mibig_compound'))) %>%
        pivot_longer(everything(), names_to='category', values_to='value') %>%
        filter(!is.na(value)) 
      kable(gene_data, 'html', escape = FALSE, col.names = NULL) %>%
        kable_styling(bootstrap_options = c("striped")) %>%
        pack_rows('gene data',  1, nrow(gene_data)) %>%
        HTML()
    } else {
      invisible()
    }
  }
  observeEvent(input$perform_enrichment,{
    output$enrichment_table = DT::renderDT({
      validate(need(!is.null(input$gene_categories_network), message = FALSE))
      enrich_res = enrichment_test(gene_list = list('network' = rv$displayed_nodes),
                                   columns_list = c('Gene Ontology', 'KEGG pathways', 'biosynthetic gene clusters', 
                                                    'Subcellular localization', 'Interpro domains'),
                                   functional_annotation %>% 
                                     filter(genome == input$dataset_network) %>%
                                     semi_join(network_genes, by='gene_id')) %>% 
        mutate(padjust = p.adjust(pval, method = 'fdr')) %>%
        filter(padjust < 0.05) %>% 
        select(-pval, -gene_list_name) %>%
        relocate(c(padjust, annotation_category), .after = all) %>% 
        rename(n_genes_subnetwork = genes, n_genes_all_network = all, pvalue_adjusted = padjust) %>% 
        arrange(pvalue_adjusted)
      enrich_res
    }, escape = FALSE, options = list(paginate=FALSE, info = FALSE, sort=FALSE, dom = 't',
                                      caption = 'Enrichment results:'), rownames = FALSE)
  })
  output$download_network_data = downloadHandler(
    filename = 'network_data.xlsx',
    content = function(file){
      df = rv$network_data
      df$edges = df$edges %>% select(from, to, weight)
      df$nodes = df$nodes %>% select(any_of(c('id','protein name','KEGG pathways','Gene Ontology','Interpro domains',
                                            'Subcellular localization','deeploc_score',
                                            'biosynthetic gene clusters','mibig_accession','mibig_compound')))
      openxlsx::write.xlsx(df, file, firstRow = TRUE)
    }
  )

  ### PCA ###
  dataset_pca = reactive({
    switch(input$dataset_pca,
           'JCVI' = read_rds('data/pca_jcvi.rds'),
           'chrom_level' = read_rds('data/pca_chrom_level.rds'))
  })
  create_pca = function(){
    pca = dataset_pca()
    percent_var = pca$sdev^2 / sum(pca$sdev^2) 
    percent_var = round(100 * percent_var, 2)
    # assemble data for biplot
    pca_scores = pca$x[,1:10] %>% as.data.frame() %>% 
      rownames_to_column('run') %>%
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
                 color=.data[[input$category_to_color_pca]], label=sample_description, key=run)) +
      geom_point(alpha=0.5) +
      xlab(paste0("PC", input$pc_x, ": ",percent_var[as.numeric(input$pc_x)],"% variance")) +
      ylab(paste0("PC", as.numeric(input$pc_x) + 1, ": ",percent_var[as.numeric(input$pc_x) + 1],"% variance")) +
      theme_bw() +
      scale_fill_manual(values = mpn65) +
      ggeasy::easy_remove_legend()
  }
  output$pca = renderPlotly({
    create_pca()
  })
  ## Create metadata table when a sample in the PCA plot is clicked
  output$sample_metadata_table_pca = DT::renderDT({
    click_data_pca = event_data('plotly_click')
    validate(need(!is.null(click_data_pca), 'Click on a sample dot to get further sample data'))
    metadata %>% filter(run == click_data_pca$key) %>%
      select(-sample_description) %>%
      mutate(across(everything(), as.character)) %>%
      pivot_longer(everything(), names_to='category', values_to='value') %>%
      filter(!is.na(value))  
  }, escape = FALSE, options = list(paginate=FALSE, info = FALSE, sort=FALSE), rownames= FALSE)
  
  ### JBrowse ###
  assembly_jcvi = assembly(paste0(jbrowse_resources_base, "jcvi.fa.gz"), bgzip = TRUE)
  annot_track_jcvi = track_feature(paste0(jbrowse_resources_base, "jcvi.gff.gz"),
                                   assembly_jcvi) %>% 
    str_replace('name\": \"jcvi\"', 'name\": \"genes\"')
  assembly_chrom_level = assembly(paste0(jbrowse_resources_base, "chrom_level.fa.gz"), bgzip = TRUE)
  annot_track_chrom_level = track_feature(paste0(jbrowse_resources_base, "chrom_level.gff.gz"),
                                          assembly_chrom_level) %>% 
    str_replace('name\": \"chrom_level\"', 'name\": \"genes\"')
  bw_track_jcvi = track_wiggle(paste0(jbrowse_resources_base, "jcvi_mean_unstranded.bw"), assembly_jcvi) %>% 
    str_replace('name\": \"jcvi_mean_unstranded\"', 'name\": \"mean read coverage\"')
  bw_track_chrom_level = track_wiggle(paste0(jbrowse_resources_base, "chrom_level_mean_unstranded.bw"), 
                                      assembly_chrom_level) %>% 
    str_replace('name\": \"chrom_level_mean_unstranded\"', 'name\": \"mean read coverage\"')
  assembly = reactive({
    switch(input$dataset_jbrowse,
           'JCVI' = assembly_jcvi,
           'chrom_level' = assembly_chrom_level
    )
  })
  track_set = reactive({
    switch(input$dataset_jbrowse,
           'JCVI' = tracks(annot_track_jcvi, bw_track_jcvi),
           'chrom_level' = tracks(annot_track_chrom_level, bw_track_chrom_level)
    )
  })
  default_session_jcvi = default_session(
    assembly_jcvi,
    c(annot_track_jcvi, bw_track_jcvi)
  )
  default_session_chrom_level = default_session(
    assembly_chrom_level,
    c(annot_track_chrom_level, bw_track_chrom_level)
  )
  default_session = reactive({
    switch(input$dataset_jbrowse,
           'JCVI' = default_session_jcvi,
           'chrom_level' = default_session_chrom_level
    )
  })
  output$jbrowse_output = renderJBrowseR(
    JBrowseR(
      "View",
      assembly = assembly(),
      tracks = track_set(),
      theme = theme('#2596be', '#2596be', '#2596be', '#2596be'),
      defaultSession = default_session(),
      location = ifelse(input$dataset_jbrowse == "JCVI", 
                        "EQ963478:2233141..2251880", "CP044620.1:5,082,369..5,101,183"),
       text_index = ifelse(input$dataset_jbrowse == "JCVI",
                           text_index(paste0(jbrowse_resources_base, "trix/jcvi.ix"),
                                      paste0(jbrowse_resources_base, "trix/jcvi.ixx"),
                                      paste0(jbrowse_resources_base, "trix/jcvi_meta.json"),
                                      "jcvi"),
                           text_index(paste0(jbrowse_resources_base, "trix/chrom_level.ix"),
                                      paste0(jbrowse_resources_base, "trix/chrom_level.ixx"),
                                      paste0(jbrowse_resources_base, "trix/chrom_level_meta.json"),
                                      "chrom_level"))
    )
  )
}