ui = navbarPage(
  # tags$style("#fieldsPanel {font-size:10px;height:10px;}
  #            #group input {height:10px;}"
  # ),
  useShinyjs(),
  id='tabs',
  title = 'Aspergillus flavus expression database',
  tabPanel('Single gene barplot ',
           sidebarLayout(
             sidebarPanel(
               selectInput(inputId = 'dataset_barplot',
                           label = 'Choose an assembly/annotation:',
                           choices = c('GCA_000006275.2 (JCVI-afl1-v2.0)' = 'JCVI', 
                                       'GCA_009017415.1 (chromosome level)' = 'chrom_level')),
               selectInput(inputId = 'normalization_method_barplot',
                           label = 'Choose normalization method:',
                           choices = c('TPM (transcripts per million)' = 'TPM', 
                                       'VST (variance stabilizing transformation)' = 'VST')),
               textInput(inputId = 'gene_id',
                         label = 'Gene ID',
                         value = 'AFLA_139360'),
               downloadButton('download_single_gene_data', label = 'Download expression data for this gene')
             ),
             mainPanel(
               plotlyOutput('barplot'),
               DT::DTOutput('sample_metadata_table')
             )
           )
  ),
  tabPanel('Multiple gene heatmap ',
           sidebarLayout(
             sidebarPanel(
               selectInput(inputId = 'dataset_heatmap',
                           label = 'Choose an assembly/annotation:',
                           choices = c('GCA_000006275.2 (JCVI-afl1-v2.0)' = 'JCVI', 
                                       'GCA_009017415.1 (chromosome level)' = 'chrom_level')),
               selectInput(inputId = 'annotation_category', 
                           label = 'Choose type of annotation', 
                           choices = annotation_categories, 
                           selected = 'Gene Ontology'),
               selectizeInput(inputId = 'gene_categories', 
                           label = 'Choose genes to include', 
                           choices = '', 
                           multiple=TRUE, 
                           options = list(
                             delimiter = ',',
                             create = I("function(input, callback){
                                return {
                                  value: input,
                                  text: input
                                 };
                              }"))),
               selectInput(inputId = 'bioprojects_to_include', 
                           label = 'Choose bioprojects to include', 
                           choices = bioproject_sizes %>% filter(n > 10) %>% 
                             mutate(bioproject_and_title = str_c(bioproject, study_title, sep=': ')) %>%
                             pull(bioproject, name = bioproject_and_title), 
                           selected = bioproject_sizes$bioproject[1:2],
                           multiple=TRUE),
               actionButton("generate_heatmap", "Generate heatmap"),
               br(),
               br(),
               downloadButton('download_multi_gene_data', label = 'Download expression data for these genes'),
               br(),
               br(),
               downloadButton('download_entire_expression_dataset', label = 'Download entire expression dataset')
             ),
             mainPanel(
               conditionalPanel(
                 condition = "input.generate_heatmap > 0",
                 InteractiveComplexHeatmapOutput(title1 = "Full heatmap", layout = "1|(2-3)", width1=900, height1=450, 
                                               action = 'click', output_ui = htmlOutput("gene_info"))
             ))
           )
  ),
  
  tabPanel(title = 'Co-expression network ',
           sidebarLayout(
             sidebarPanel(
               selectInput(inputId = 'dataset_network',
                           label = 'Choose an assembly/annotation:',
                           choices = c('GCA_000006275.2 (JCVI-afl1-v2.0)' = 'JCVI', 
                                       'GCA_009017415.1 (chromosome level)' = 'chrom_level')),
               selectInput(inputId = 'annotation_category_network', 
                           label = 'Choose type of annotation', 
                           choices = annotation_categories, 
                           selected = 'Gene Ontology'),
               selectizeInput(inputId = 'gene_categories_network', 
                              label = 'Choose genes to include', 
                              choices = '', 
                              multiple=TRUE, 
                              options = list(
                                delimiter = ',',
                                create = I("function(input, callback){
                                return {
                                  value: input,
                                  text: input
                                 };
                              }"))),
               checkboxInput(inputId = 'show_neighbors', 
                             label = 'Show neighbors', 
                             value = TRUE),
               actionButton("generate_network", "Generate network"),
               br(),
               br(),
               actionButton("perform_enrichment", "Enrichment analysis"),
               br(),
               br(),
               downloadButton('download_network_data', label = 'Download network data for shown genes')
             ),
             mainPanel(
               visNetworkOutput("network_vis", height = '600px'),
               htmlOutput("node_data_from_network"),
               DT::DTOutput('enrichment_table')
             )
           )
  ),
  tabPanel('PCA ',
           sidebarLayout(
             sidebarPanel(
               selectInput(inputId = 'dataset_pca',
                           label = 'Choose an assembly/annotation:',
                           choices = c('GCA_000006275.2 (JCVI-afl1-v2.0)' = 'JCVI', 
                                       'GCA_009017415.1 (chromosome level)' = 'chrom_level')),
               selectInput(inputId = 'pc_x',
                           label = 'principal component to display',
                           choices = seq(1, 10),
                           selected = 1),
               selectizeInput(inputId = 'category_to_color_pca', 
                              label = 'Category to color samples by', 
                              choices = c('bioproject', 'strain', 'growth medium'))
             ),
             mainPanel(
               plotlyOutput('pca'),
               DT::DTOutput('sample_metadata_table_pca')
             )
           )
  ),
  tabPanel('Genome browser ',
           sidebarLayout(
            sidebarPanel(
              selectInput(inputId = 'dataset_jbrowse',
                          label = 'Choose an assembly/annotation:',
                          choices = c('GCA_000006275.2 (JCVI-afl1-v2.0)' = 'JCVI', 
                                      'GCA_009017415.1 (chromosome level)' = 'chrom_level'))
            ),
            mainPanel(
              JBrowseROutput("jbrowse_output")
            )
  )),
  tabPanel('About',
           tags$p('This expression database was made using the Aspergillus flavus RNA-Seq datasets available on NCBI. 
                  Details of the methods used to create the database and interface can be found at ', 
                  tags$a('https://github.com/kcamnairb/flavus_sra_rnaseq', href = 'https://github.com/kcamnairb/flavus_sra_rnaseq', target='_blank'),
                  '. Each tab has the option to see the data produced using reads mapped to two different genome assemblies of Aspergillus flavus NRRL 3357. 
                  GCA_000006275.2 (JCVI-afl1-v2.0) is an older assembly submitted in 2009 that has been used the most frequently and most gene ids in the literature refer to this assembly. 
                  GCA_009017415.1 is a chromosome-level assembly that was submitted in 2019.'),
           tags$h5('Single gene barplot'),
           tags$p('This tab allows you to query any single gene and see the expression values
                  plotted using either TPM or VST normalization for all samples.
                  The samples are colored by their bioproject and further details about the sample can be accessed by clicking on the bar for each sample.'),
           tags$h5('Multiple gene heatmap'),
           tags$p('This tab allows you create a heatmap from either a user provided list of gene ids that are comma separated, 
                  for example "AFLA_023020,AFLA_023020...", or you can choose groups of genes by functional annotation using the dropdown menus.
                  When the type of annotation is changed, the list of terms in the dropdown menu below it will be updated. 
                  When multiple functional categories are selected, the category each gene belongs to is indicated by a different color on the left of the heatmap.
                  To provide a comma separated list of gene ids, first select "Gene list (Comma separated)" in the type of annotation menu.
                  The samples included in the heatmap can be selected by bioproject in the final dropdown menu. 
                  When multiple bioprojects are selected, each bioproject will be shown in a separate group since there can be technical biases specific to each project.
                  Once the genes and bioprojects are selected click "Generate heatmap" to produce the visualization.
                  Clicking on a cell in the heatmap will display the functional annotation for the gene and metadata for the sample.
                  If the heatmap is too dense to see the text on the axes properly, drawing a rectangle on a portion of the heatmap will display a sub-heatmap below.
                  The expression values and functional data for the selected dataset or the entire dataset can be downloaded to an Excel file using the two buttons at the bottom left.
                  '),
           tags$h5('Co-expression network'),
           tags$p('This tab allows you to see genes that are positively and negatively correlated with your genes of interest using the same dropdown menus as the heatmap tab to select the genes.
                  After the subnetwork containing your genes of interest is displayed you can click on the "Enrichment analysis" button to look for functional
                  terms that are enriched in the subnetwork using a one-sided Fishers exact test. Clicking on individual nodes in the network will show information about the corresponding gene.
                  The edge weights for the subnetwork and functional annotation of the genes can be downloaded by clicking the button on the bottom left.'),
           tags$h5('PCA'),
           tags$p('This tab shows a principal component analysis of the samples using VST counts as input. 
                  Additional principal components can be shown using the second drop down menu.
                  The color of the samples defaults to the bioproject, but can also be changed to the strain or growth medium.'),
           tags$h5('Genome browser'),
           tags$p('This tab displays a JBrowse genome browser showing an RNA-Seq track that is a mean of the coverage of all the samples.
                  Entering your gene of interest in the search bar will navigate the browser to the locus of that gene.
                  Looking at the coverage of a gene can allow you to evaluate the accuracy of a gene prediction.')
  ),
  theme = bs_theme(bootswatch = 'cerulean')
)
