ui = navbarPage( 
  tags$style("#fieldsPanel {font-size:10px;height:10px;}
             #group input {height:10px;}"
  ),
  title = 'Aspergillus flavus expression database',
  #tabPanel('Single gene barplot',
  #         sidebarLayout(
  #           sidebarPanel(
  #             selectInput(inputId = 'dataset',
  #                         label = 'Choose an assembly/annotation:',
  #                         choices = c('GCA_000006275.2 (JCVI-afl1-v2.0)', 'GCA_009017415.1 (chromosome level)')),
  #             textInput(inputId = 'gene_id',
  #                       label = 'Gene ID',
  #                       value = 'AFLA_139360'),
  #             downloadButton('download_single_gene_data', label = 'Download expression data for this gene')
  #           ),
  #           mainPanel(
  #             plotlyOutput('barplot'),
  #             dataTableOutput('sample_metadata_table')
  #           )
  #         )
  #),
  tabPanel('Multiple gene heatmap',
           sidebarLayout(
             sidebarPanel(
               selectInput(inputId = 'dataset',
                           label = 'Choose an assembly/annotation:',
                           choices = c('GCA_000006275.2 (JCVI-afl1-v2.0)', 'GCA_009017415.1 (chromosome level)')),
               selectInput(inputId = 'normalization_method',
                           label = 'Choose an normalization method:',
                           choices = c('TPM (transcripts per million)', 'VST (variance stabilizing transformation)')),
               selectInput(inputId = 'annotation_category', 
                           label = 'Choose type of annotation', 
                           choices = annotation_categories, 
                           selected = 'Gene Ontology'),
               selectizeInput(inputId = 'gene_categories', 
                           label = 'Choose types of genes to include', 
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
               downloadButton('download_multi_gene_data', label = 'Download expression data for these genes')
             ),
             mainPanel(
               conditionalPanel(
                 condition = "input.generate_heatmap > 0",
                 InteractiveComplexHeatmapOutput(title1 = "Full heatmap", layout = "1|(2-3)", width1=1000, height1=450, 
                                               action = 'click', output_ui = htmlOutput("gene_info"))
             ))
           )
  ),
  theme = bs_theme(bootswatch = 'cerulean')
)
