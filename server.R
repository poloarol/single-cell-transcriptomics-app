

if (!require('pacman')) install.packages("pacman")

# Load contributed packages with pacman
pacman::p_load(pacman, Seurat, tidyverse, shiny)


load_data <- function(project, cells = 3, features){
    # Load the PBMC dataset
    reads <- Read10X(data.dir = "data/pbmc3k/filtered_gene_bc_matrices/hg19/")
    # Initialize the Seurat object with the raw (non-normalized data).
    data <- CreateSeuratObject(
    counts = reads, 
    project = project, 
    min.cells = as.numeric(cells), 
    min.features = as.numeric(features))
    
    print(data)
    
    
}


metricsplot <- function(data){
  
  # seurat_object <- load_data(project, cells, features)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  plt <- VlnPlot(data, 
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 ncol = 3)
  
  return(plt)
}

featureplot <- function(data){
  
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  return (plot1 + plot2)
  
}


server <- function(input, output, session) {
  
  seurat_obj <- reactive({
    load_data(input$proj_name, input$min.cells, input$min.feats)
    })
  
  output$metrics <- renderPlot(metricsplot(seurat_obj))
  output$features <- renderPlot(featureplot(seurat_obj))
  
}