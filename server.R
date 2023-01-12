

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
            min.cells = num_cells,
            min.features = num_features)
    
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

normalize_data <- function(data, strategy = "LogNormalize", feats = 200){
  results <-  NormalizeData(data,  normalization.method = strategy, nfeatures = feats)
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
  
  observed <- reactiveValues(norm = NULL, hv = NULL)
  
  seurat_obj <- eventReactive(
    eventExpr = 
      {
      input$proj.name
      input$min.cells
      input$min.feats
      }, {
        load_data(
          project_name = input$proj.name,
          num_cells = input$min.cells,
          num_features = input$min.feats)
  })
  
  output$metrics <- renderPlot(metricsplot(seurat_obj()))
  output$features <- renderPlot(featureplot(seurat_obj()))

  subset_seurat <- eventReactive(
    eventExpr =
      {
        input$subset
      }, {
        subset_dataset(
          data = seurat_obj(),
          min_rna = as.numeric(input$min.genes),
          max_rna = as.numeric(input$max.genes),
          mito = as.numeric(input$mito.pcts))
      }
  )
  
  observe({
    if(input$subset){
      print("Doublet removal")
    }
  })
  
  observe({
    if(input$run.norm){
      data <- normalize_data(subset_seurat(), input$normalization, input$nfeatures)
      observed$hv <- variable_features(data, input$ftselection)
    }
  })
  
  hv.plot <- reactive({variable_features_plot(observed$hv)})
  output$topvariable <- renderPlot(hv.plot())

  
  scaled_seurat <- reactive({
    all.genes <- rownames(observed$hv)
    ScaleData(observed$hv, all.genes)
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
      dims = 1:input$num.dim
    )
  })


  umap_seurat <- reactive({RunUMAP(cluster_seurat(), dims = 1:input$num.dim)})
  umap_plot <- reactive({DimPlot(umap_seurat(), reduction = "umap")})
  output$umap <- renderPlot(umap_plot())
  
  tsne_seurat <- reactive({RunTSNE(cluster_seurat(), dims = 1:input$num.dim)})
  tsne_plot <- reactive({DimPlot(tsne_seurat(), reduction = "tsne")})
  output$tsne <- renderPlot(tsne_plot())
  

  biomarkers <- reactive({
    markers <- FindAllMarkers(umap_seurat(), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    df_grouped <- split(markers, markers$cluster)
    top2_rows <- lapply(df_grouped, function(x) x[order(x$avg_log2FC, decreasing = TRUE)[1:2],])
    df_result <- do.call(rbind, top2_rows)
    df <- df_result[c("gene","avg_log2FC","p_val_adj")]
    df <- df[order(-df$avg_log2FC),]
    colnames(df) <- c("Gene", "Average Log2 Fold Change","Adjusted P-Value")
    df
  })

  output$biomarkers <- DT::renderDataTable({biomarkers()}, rownames = FALSE)
  
  gene.list <- reactive({
    if(is.null(input$gene.list)){
      df <- biomarkers()
      df$Gene[1:4]
    }else{
      unlist(strsplit(input$gene.list, ","))
    }
  })
  
  output$deplot <- renderPlot(VlnPlot(umap_seurat(), features = gene.list()))

  # output$heatmap <- renderPlot({
  #   df_grouped <- split(biomarkers(), biomarkers$cluster)
  #   top10_rows <- lapply(df_grouped, function(x) x[order(x$avg_log2FC, decreasing = TRUE)[1:10],])
  #   top10 <- do.call(rbind, top10_rows)
  #   DoHeatmap(umap_seurat(), features = top10$gene) + NoLegend()
  # })

}






