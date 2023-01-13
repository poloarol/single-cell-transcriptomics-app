

if (!require('pacman')) install.packages("pacman")

pacman::p_load(pacman, tidyverse, shiny, DT, shinycssloaders, shinyFiles, shinyWidgets)

options(shiny.maxRequestSize = 30*1024^2)

ui <- fluidPage(
  
  titlePanel("A Shiny app to facilate analysis of single-cell transcriptomics data"),
  
  navbarPage(
    "single-cell app",
    tabPanel("Single-Cell RNAseq clustering",
             br(),
             mainPanel(width=12,
                       tabsetPanel(type = "pills",
                                   br(),
                                   tabPanel("Quality Control & Filtering",
                                            sidebarLayout(
                                              sidebarPanel(
                                                width = 3,
                                                fileInput("rna1", "Load 10X Genomics data",
                                                          multiple = TRUE,
                                                          accept = c(".csv", ".tsv", ".mtx", 
                                                                     ".RDS", ".HDF5", ".loom", "h5ad")),
                                                h4("Adjust parameters for MTX dataset"),
                                                textInput("proj.name", "Project name",value = "abc13"),
                                                numericInput("min.cells", "Min. number of cells", 5),
                                                numericInput("min.feats", "Min. number of features", 200),
                                                br(),
                                                h4("Adjust parameters to remove low quality cells and empty droplets"),
                                                numericInput("min.genes", "Min. number of genes", 200),
                                                numericInput("max.genes", "Max. number of genes", 2500),
                                                numericInput("mito.pcts", "Percentage of Mitochonrial genes", 5),
                                                actionButton("subset", "Process")
                                              ),
                                              mainPanel(
                                                br(), br(), br(),
                                                 fluidRow(
                                                   column(6, align = "center", h4("QC Metrics"), plotOutput("metrics") %>% withSpinner(color="#0dc5c1")),
                                                   column(6, align = "center", h4("Feature-Feature Relationships"), plotOutput("features") %>% withSpinner(color="#0dc5c1"))
                                                 )))),
                                   tabPanel("Identification of highly variable genes",
                                            sidebarLayout(
                                              sidebarPanel(
                                                width = 3,
                                                h4("Normalization Parameters"),
                                                br(), br(),
                                                div(class="option-group",
                                                    numericInput("nfeatures", "Number of features", 2000),
                                                    radioButtons("normalization",
                                                                 "Normalization Strategy",
                                                                 choices = c("LogNormalize", "CLR", "RC"), inline = TRUE),
                                                    radioButtons("ftselection",
                                                                 "Feature Selection",
                                                                 choices = c("vst", "mvp", "disp"), inline = TRUE)),
                                                actionButton("run.norm", "Run Normalization")),
                                              mainPanel(
                                                br(), br(),
                                                column(12, plotOutput("topvariable") %>% withSpinner(color="#0dc5c1"))))),
                                   tabPanel("Linear dimensionality reduction",
                                            fluidRow(
                                              column(4,
                                                     align = "center",
                                                     h4("Viisualize PCA"),
                                                     plotOutput("pca") %>% withSpinner(color="#0dc5c1")),
                                              column(4,
                                                     align = "center",
                                                     h4("Number of significant PCs"),
                                                     plotOutput("jack") %>% withSpinner(color="#0dc5c1")),
                                              column(4,
                                                     align = "center",
                                                     h4("Elbow Plot: Ranking of PCs by % of variance explained"),
                                                     plotOutput("elbow") %>% withSpinner(color="#0dc5c1")))),
                                   tabPanel("Clustering and non-linear dimensionaliry reduction",
                                            sidebarLayout(
                                              sidebarPanel(
                                                width = 3,
                                                h4("Dimensionality Reduction Parameters"),
                                                numericInput("num.dim", "Number of dimensions", 10),
                                                sliderInput("range", "UMAP Resolution:",min = 0, max = 100, value = 50),
                                                selectInput("ref", "Reference organism", 
                                                            choices = c("", "HumanPrimaryCellAtlasData", "BlueprintEncodeData",
                                                                        "MouseRNAseqData", "ImmGenData", "DatabaseImmuneCellExpressionData",
                                                                        "NovershternHematopoieticData", "MonacoImmuneData")),
                                                selectInput("algo", "Clustering Algorithm", 
                                                            choices = c("Louvain algorithm",
                                                                        "Louvain algorithm with Multivelel Refinement",
                                                                        "SLM Algorithm", "Leiden Algorithm")),
                                              ),
                                            mainPanel(
                                              fluidRow(
                                                column(3),
                                                column(6,
                                                       align = "center",
                                                       h4("Clustering of cells using UMAP"),
                                                       plotOutput("umap") %>% withSpinner(color="#0dc5c1")))))),
                                   tabPanel("Differentially Expressed Genes",
                                            sidebarLayout(
                                              sidebarPanel(
                                                width = 3,
                                                textInput("gene.list", "Enter gene name", value = "")
                                            ),
                                            fluidRow(
                                              column(4,
                                                     align = "center",
                                                     h4("Top Expressed biomarkers"),
                                                     DT::dataTableOutput("biomarkers") %>% withSpinner(color="#0dc5c1")),
                                              column(4,
                                                     align = "center",
                                                     plotOutput("deplot") %>% withSpinner(color="#0dc5c1"))))))),
  ),
  tabPanel("Single-Cell ATACseq workflow",
           mainPanel()),
  tabPanel("Single-Cell RNA-seq Integration",
           mainPanel()),
  tabPanel("Single-Cell ATAC-seq Integration",
           mainPanel()),
  )
)