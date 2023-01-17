



if (!require('pacman')) install.packages("pacman")

# Load contributed packages with pacman
pacman::p_load(tools, pacman, Seurat, SeuratObject, tidyverse, shiny, DT, shinyFiles, shinyWidgets)

# source("src/doublet_removal.R")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("SingleR")
# BiocManager::install("celldex")
# BiocManager::install("SingleCellExperiment")

require(SingleR)
require(celldex)
require(SingleCellExperiment)


# Load the Seurat object
load_data <- function(reads, project_name, num_cells = 3, num_features = 200){
  datum <- CreateSeuratObject(
    counts = reads,
    project = project_name,
    names.delim = "-",
    min.cells = num_cells,
    min.features = num_features)
  
  datum[["percent.mt"]] <- PercentageFeatureSet(datum, pattern = "^MT-")
  return(datum)
}

# Provide metrics plot for RNA seq
metricsplot <- function(data){
  plt <- VlnPlot(data, 
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 ncol = 3)
  
  return(plt)
}

# provide feature plot for RNAseq
featureplot <- function(data){
  plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  return (plot1 + plot2)
  
}

# Subset the data based on min. number of genes, max. number of genes,
# percentage of mRNA and number of features
subset_dataset <- function(data, min_rna = 200, max_rna = 2500, mito = 5){
  results <- subset(data, 
                    subset = nCount_RNA > min_rna & nFeature_RNA < max_rna &
                      percent.mt < 5)
  return(results)
}

# Normalize the dataset
normalize_data <- function(data, strategy = "LogNormalize", feats = 200){
  results <-  NormalizeData(data,  normalization.method = strategy, nfeatures = feats)
  return(results)
}

# Find the variable features or differentially expressed genes
variable_features <- function(data, strategy = "vst", features = 2000){
  results <- FindVariableFeatures(data, 
                                  selection.method = strategy, 
                                  nfeatures = features)
  return(results)
}

# Produce a variable feature plot of top10 most expressed genes
variable_features_plot <- function(data){
  top10 <- head(VariableFeatures(data), 10)
  plt1 <- VariableFeaturePlot(data)
  plt2 <- LabelPoints(plot = plt1, points = top10, repel = TRUE)
  plt <- plt1 + plt2
  
  return(plt)
}

# Annonate cellular data
add_cell_data <- function(data, reference){
  result <- SingleR(test = as.SingleCellExperiment(data), ref = reference, labels = reference$label.main)
  data$singlr_labels <- result$labels
  return(data)
}


# Run RNAseq clustering pipeline
rnaseq_analysis <- function(input, output, session){
  
  observed <- reactiveValues(rna1 = NULL, norm = NULL, hv = NULL)
  
  seurat_obj <- eventReactive(
    eventExpr = {
      input$rna1
    }, {
      if(length(input$rna1$datapath) > 1){
        dt_read <- ReadMtx(
          mtx = input$rna1$datapath[3],
          cells = input$rna1$datapath[1],
          features = input$rna1$datapath[2]
        )
        data <- load_data(dt_read, input$proj.name, input$min.cells, input$min.feats)      
      } else {
        extension = tolower(file_ext(input$rna1$datapath[1]))
        if(extension == "loom"){
          dt_read <-as.Seurat(Connect(filename = input$rna1$datapath[1], mode = "r"))
          data <- load_data(dt_read, input$proj.name, input$min.cells, input$min.feats) 
        }else if(extension == "h5ad"){
          obj <- Connect(filename = input$rna1$datapath[1], des="h5seurat", overwrite = TRUE)
          dt_read <- LoadH5Seurat(obj)
          data <- load_data(dt_read, input$proj.name, input$min.cells, input$min.feats)
        }else if(extension == "hdf5"){
          obj <- Read10X_h5(input$rna1$datapath[1])
          load_data(obj, input$proj.name, input$min.cells, input$min.feats)
        }else if(extension == "rds"){
          readRDS(input$rna1$datapath[1])
          data <- load_data(dt_read, input$proj.name, input$min.cells, input$min.feats)
        }else {
          print("Enter an appropriate format")
        }
      }
      
      
    }
  )
  
  
  output$metrics <- renderPlot(metricsplot(seurat_obj()))
  output$features <- renderPlot(featureplot(seurat_obj()))
  
  subset_seurat <- eventReactive(
    eventExpr =
      {
        input$min.genes
        input$max.genes
        input$mito.pcts
      }, {
        subset_dataset(
          data = seurat_obj(),
          min_rna = as.numeric(input$min.genes),
          max_rna = as.numeric(input$max.genes),
          mito = as.numeric(input$mito.pcts))
      }
  )
  
  datasetInput <- reactive({
    switch(input$ref,
           "HumanPrimaryCellAtlasData" = celldex::HumanPrimaryCellAtlasData(),
           "BlueprintEncodeData" = celldex::BlueprintEncodeData(),
           "MouseRNAseqData" = celldex::MouseRNAseqData(),
           "ImmGenData" = celldex::ImmGenData(),
           "DatabaseImmuneCellExpressionData" = celldex::DatabaseImmuneCellExpressionData(),
           "NovershternHematopoieticData" = celldex::NovershternHematopoieticData(),
           "MonacoImmuneData" = celldex::MonacoImmuneData()
    )
  })
  
  
  algorithm <- reactive({
    switch(input$algo,
           "Louvain Algorithm" = 1,
           "Louvain algorithm with Multivelel Refinement" = 2,
           "SLM Algorithm" = 3,
           "Leiden Algorithm" = 4
           )
  })
  
  observe({
    if(input$subset){
      print(seurat_obj())
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
      dims = 1:input$num.dim,
      algorithm = algorithm()
    )
  })
  
  
  umap_seurat <- reactive({RunUMAP(cluster_seurat(), dims = 1:input$num.dim)})
  tsne_seurat <- reactive({RunTSNE(cluster_seurat(), dims = 1:input$num.dim)})
  
  umap_plot <- reactive({
    if(is.null(datasetInput())){
      DimPlot(umap_seurat(), reduction = "umap")
    }else{
      dt <- add_cell_data(umap_seurat(), datasetInput())
      DimPlot(dt, group.by = "singlr_labels", reduction = "umap", label = TRUE)
    }
  })
  
  
  tsne_plot <- reactive({
    if(is.null(datasetInput())){
      DimPlot(tsne_seurat(), reduction = "tsne")
    }else{
      dt <- add_cell_data(tsne_seurat(), datasetInput())
      DimPlot(dt,group.by = "singlr_labels", reduction = "tsne")
    }
    
  })
  
  output$umap <- renderPlot(umap_plot())
  output$tsne <- renderPlot(tsne_plot())
  
  
  
  biomarkers <- reactive({
    markers <- FindAllMarkers(umap_seurat(), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    df_grouped <- split(markers, markers$cluster)
    top2_rows <- lapply(df_grouped, function(x){x[order(x$avg_log2FC, decreasing = TRUE)[1:2],]})
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


