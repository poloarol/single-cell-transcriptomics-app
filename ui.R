

if (!require('pacman')) install.packages("pacman")

pacman::p_load(pacman, tidyverse, shiny, DT)

ui <- fluidPage(
  
  titlePanel("A Shiny app to facilate analysis of single-cell transcriptomics data"),
  
  navbarPage(
    "single-cell transcriptomics",
    tabPanel("Single-Cell Transcriptomics",
             mainPanel(width=12,
               
               h2("1. Quality Control and Cell Selection"),
               
               fluidRow(align="center",
                 column(4,
                        align = "center",
                        h4("Project Name"),
                        textInput("proj_name", NULL, value = "pbmc3k", width = NULL, placeholder = NULL)
                 ),
                 column(4,
                        align = "center",
                        h4("Min. Number of Cells"),
                        textInput("min.cells", NULL, value = 10, width = NULL, placeholder = NULL)
                 ),
                 column(4,
                        align = "center",
                        h4("Min. Number of features"),
                        textInput("min.feats", NULL, value = 200, width = NULL, placeholder = NULL)
                 )
               ),
               
               fluidRow(
                 column(6, align = "center", "QC Metrics",plotOutput("metrics")),
                 column(6, align = "center", "Feature-Feature Relationships", plotOutput("features"))
               ),
               
               h2("2. Filtering and Doublet removal"),
              
              fluidRow(
                column(3,
                       align = "center",
                       h4("Min. Number of Genes Expressed"),
                       textInput("min_cells", NULL, value = 200, width = NULL, placeholder = NULL)
                ),
                column(3,
                       align = "center",
                       h4("Max. number of Genes Expressed"),
                       textInput("max_cells", NULL, value = 2500, width = NULL, placeholder = NULL)
                ),
                column(3,
                       align = "center",
                       h4("Percentage Mitochondrial Genes"),
                       textInput("mt", NULL, value = 5, width = NULL, placeholder = NULL)
                ),
                column(3,
                       h4("Dealing with Doublets"),
                       checkboxInput("dblt", "Doublet Removal"))
              ),
              
              h2("3. Identification of highly variable genes"),
              
              fluidRow(
                column(4,
                       h4("Parameters"),
                       div(class="option-group",
                           radioButtons("normalization",
                                        "Normalization Strategy",
                                        choices = c("LogNormalize", "CLR", "RC"),
                                        inline = TRUE),
                           radioButtons("ftselection",
                                        "Feature Selection",
                                        choices = c("vst", "mvp", "disp"),
                                        inline = TRUE
                                        )
                           # checkboxGroupInput("scaling-data",
                           #              "Scaling Data Options",
                           #              choices = c("Speed Up", "Remove unwanted sources of variation"))
                           )),
                       # actionButton('process-data-scaling', "Process Data")),
                column(8,
                       align = "center",
                       h4("Top 10 most variable genes"),
                       plotOutput("topvariable"))
              ),
              
              h2("4. Linear Dimensional Reduction (PCA)"),
              
              fluidRow(
                column(4,
                       align = "center",
                       h4("Visualise PCA"),
                       plotOutput("pca")),
                column(4,
                       align = "center",
                       h4("Identify the number of significant PCs"),
                       plotOutput("jack")),
                column(4,
                       align = "center",
                       h4("Elbow Plot: Ranking of PCs by % of variance explained"),
                       plotOutput("elbow"))
              ),
              
              h2("5. Find Clusters and Run non-linear dimensionality reduction"),
              
              fluidRow(
                column(2),
                column(4,
                  align = "center",
                  h4("UMAP Resolution"),
                  sliderInput("range", "Resolution:",min = 0, max = 100, value = 50)
                ),
                column(4,
                       align = "center",
                       h4("UMAP Clusters"),
                       plotOutput("umap")
                ),
                column(2)
              ),
              
              h2("6. Identification of Differential Expressed features"),
              
              fluidRow(
                column(2,
                       align = "center",
                       h4("Top Expressed biomarkers"),
                       DT::dataTableOutput("biomarkers")),
                column(4,
                       align = "center",
                       h4("Violin Plot of Top 4 Expressed biomarkers across clusters"),
                       plotOutput("deplot")),
                column(6,
                       align = "center",
                       h4("Expression HeatMap of top 20 biomarkers"),
                       plotOutput("heatmap"))
              ),
              
              
              )),
    tabPanel("Single-Cell Intergration",
             mainPanel()),
    tabPanel("Spatial Transcriptomics Analysis",
             mainPanel()),
    navbarMenu("More",
               "Repository")
  )
)