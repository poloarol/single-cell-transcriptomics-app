

if (!require('pacman')) install.packages("pacman")

# Load contributed packages with pacman
pacman::p_load(pacman, Seurat, SeuratObject, tidyverse, shiny, DT)


load_data <- function(project_name, num_cells = 3, num_features = 200){
    # Load the PBMC dataset
    reads <- Read10X(data.dir = "data/pbmc3k/filtered_gene_bc_matrices/hg19/")
    # Initialize the Seurat object with the raw (non-normalized data).
    datum <- CreateSeuratObject(
            counts = reads,
            project = project_name,
            min.cells = as.numeric(num_cells),
            min.features = as.numeric(num_features))
    
    datum[["percent.mt"]] <- PercentageFeatureSet(datum, pattern = "^MT-")
    
    return(datum)
}


metricsplot <- function(data){
  
  # seurat_object <- load_data(project, cells, features)
  plt <- VlnPlot(data, 
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 ncol = 3)
  
  return(plt)
}

featureplot <- function(data){
  plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  return (plot1 + plot2)
  
}

subset_dataset <- function(data, min_rna = 200, max_rna = 2500, mito = 5){
  results <- subset(data, 
                    subset = nCount_RNA > min_rna & nFeature_RNA < max_rna &
                      percent.mt < 5)
  return(results)
}

normalize_dt <- function(data, strategy = "LogNormalize"){
  
  results <-  NormalizeData(data,  normalization.method = strategy)
  return(results)
}

variable_features <- function(data, strategy = "vst", features = 2000){
  results <- FindVariableFeatures(data, 
                                  selection.method = strategy, 
                                  nfeatures = features)
  return(results)
}

variable_features_plot <- function(data){
  top10 <- head(VariableFeatures(data), 10)
  plt1 <- VariableFeaturePlot(data)
  plt2 <- LabelPoints(plot = plt1, points = top10, repel = TRUE)
  plt <- plt1 + plt2
  
  return(plt)
}

server <- function(input, output, session) {
  
  seurat_obj <- eventReactive(
    eventExpr = 
      {
      input$proj_name
      input$min.cells
      input$min.feats
      }, {
        load_data(
          project_name = input$proj_name,
          num_cells = input$min.cells,
          num_features = input$min.feats)
  })
  
  output$metrics <- renderPlot(metricsplot(seurat_obj()))
  output$features <- renderPlot(featureplot(seurat_obj()))

  subset_seurat <- eventReactive(
    eventExpr =
      {
        input$min_cells
        input$max_cells
        input$mt
        input$dblt
      }, {
        subset_dataset(
          data = seurat_obj(),
          min_rna = as.numeric(input$min_cells),
          max_rna = as.numeric(input$max_cells),
          mito = as.numeric(input$mt))
      }
  )
  
  nm_seurat <- eventReactive(
    eventExpr = {
      input$normalization
    },{
      normalize_dt(subset_seurat(), input$normalization)
    })
  
  hv_seurat <- eventReactive(
    eventExpr = {
      input$ftselection
    }, {
      variable_features(nm_seurat(), input$ftselection)
    }
  )

  plt <- reactive({
    variable_features_plot(hv_seurat())
  })
    
  output$topvariable <- renderPlot(plt())
  
  
  scaled_seurat <- reactive({
    all.genes <- rownames(hv_seurat())
    ScaleData(hv_seurat(), all.genes)
  })
  
  pca_seurat <- reactive({
    RunPCA(scaled_seurat(), features = VariableFeatures(object = scaled_seurat()))
  })
  
  jack_seurat <- reactive({
    data <- JackStraw(pca_seurat(), num.replicate = 100)
    data <- ScoreJackStraw(data, dims = 1:20)
  })
  
  output$pca <- renderPlot(DimPlot(pca_seurat(), reduction = "pca"))
  
  output$jack <- renderPlot(JackStrawPlot(jack_seurat(), dims = 1:15))
  
  elbow_plot <- reactive({ElbowPlot(jack_seurat(), ndims=20, reduction = "pca")})
  output$elbow <- renderPlot(elbow_plot())
  
  cluster_seurat <- reactive({
    FindClusters(
      FindNeighbors(jack_seurat(), resolution = input$range/100),
      dims = 1:10
    )
  })
  
  
  umap_seurat <- reactive({RunUMAP(cluster_seurat(), dims = 1:10)})
  umap_plot <- reactive({DimPlot(umap_seurat(), reduction = "umap")})
  output$umap <- renderPlot(umap_plot())
  
  biomarkers <- reactive({
    FindAllMarkers(umap_seurat(), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
      group_by(clusters) %>%
      slice_max(n=2, order_by = avg_log2FC)
  })

  output$biomarkers <- DT::renderDataTable({biomarkers()})
  output$deplot <- renderPlot(VlnPlot(umap_seurat(), features = c("MS4A1", "CD79A")))
  
  output$heatmap <- renderPlot({
    biomarkers() %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
    DoHeatmap(umap_seurat(), features = top10$gene) + NoLegend()
  })

}






