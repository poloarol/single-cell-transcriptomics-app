


if (!require('pacman')) install.packages("pacman")

# Load contributed packages with pacman
pacman::p_load(pacman, Seurat, SeuratObject, Signac, 
               tidyverse, shiny, DT, shinyFiles, shinyWidgets,
               Rsamtools, patchwork, biovizBase, GenomeInfoDb,
               EnsDb.Hsapiens.v75, EnsDb.Hsapiens.v79, hdf5r, 
               EnsDb.Mmusculus.v79, BiocManager)

load_atac <- function(matrix, scfile, fragments, min.cells = 5, min.features = 200){
  
  metadata <- read.csv(
    file = scfile,
    header = TRUE,
    row.names = 1
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = matrix,
    sep = c(":", "-"),
    genome = 'hg19',
    fragments = fragments,
    min.cells = min.cells,
    min.features = min.features
  )
  
  data <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  
  return(data)
}

load_rnaseq <- function(path){
  data <- readRDS(path)
  return(data)
}

add_annotations <- function(data, db = EnsDb.Hsapiens.v75){
  annotations <- GetGRangesFromEnsDb(ensdb = db)
  seqlevelsStyle(annotations) <- 'UCSC'
  Annotation(data) <- annotations
  return(data)
}

filter_dataset <- function(data, min.peak = 3000, max.peak = 20000, pct.reads = 15, 
                           pct.blacklist = 5, nuc.signal = 4, tss.enrich = 2){
  obj <- subset(
    x = data,
    subset = peak_region_fragments > min.peak &
      peak_region_fragments < max.peak &
      pct_reads_in_peaks > pct.reads &
      blacklist_ratio < pct.blacklist /100 &
      nucleosome_signal < nuc.signal &
      TSS.enrichment > tss.enrich
  )
  return(obj)
}

normalise_dataset <- function(data, enrichment = 2, signal = 4){
  
  # compute nucleosome signal score per cell
  data <- NucleosomeSignal(object = data)
  
  # compute TSS enrichment score per cell
  data <- TSSEnrichment(object = data, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100
  data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments
  
  data$high.tss <- ifelse(data$TSS.enrichment > enrichment, 'High', 'Low')
  
  minimum <- paste("NS >", as.character(signal))
  maximum <- paste("NS < ", as.character(signal))
  data$nucleosome_group <- ifelse(data$nucleosome_signal > signal, minimum, maximum)
  
  return(data)
}

linear_dim_reduction <- function(data, min.cutoff = "q0"){
  obj <- RunTFIDF(data)
  obj <- FindTopFeatures(obj, min.cutoff = min.cutoff)
  obj <- RunSVD(obj)
  return(obj)
}

nonlinear_dim_reduction <- function(data, algorithm = 3, ndims = 30){
  obj <- RunUMAP(object = data, reduction = 'lsi', dims = 2:ndims)
  obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:ndims)
  obj <- FindClusters(object = obj, verbose = FALSE, algorithm = algorithm)
  
  return(obj)
}

get_gene_activities <- function(data, strategy = "LogNormalize"){
  gene.activities <- GeneActivity(data)
  
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  data[['RNA']] <- CreateAssayObject(counts = gene.activities)
  obj <- NormalizeData(
    object = data,
    assay = 'RNA',
    normalization.method = strategy,
    scale.factor = median(data$nCount_RNA)
  )
  
  DefaultAssay(obj) <- 'RNA'
  
  return(obj)
  
}

integrate_rnaseq <- function(seurat_atac, seurat_rna, ndims = 30){
  
  transfer.anchors <- FindTransferAnchors(
    reference = seurat_rna,
    query = seurat_atac,
    reduction = 'cca'
  )
  
  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = seurat_rna$celltype,
    weight.reduction = seurat_atac[['lsi']],
    dims = 2:ndims
  )
  
  data <- AddMetaData(object = seurat_atac, metadata = predicted.labels)
  
  return(data)
  
}

diff_accessibility <- function(data){
  # change back to working with peaks instead of gene activities
  DefaultAssay(data) <- 'peaks'
  
  da_peaks <- FindMarkers(
    object = data,
    ident.1 = "CD4 Naive",
    ident.2 = "CD14 Mono",
    test.use = 'LR',
    latent.vars = 'peak_region_fragments'
  )
  
  return(da_peaks)
  
}


atacseq_analysis <- function(input, output, session){

  peaks <- Read10X_h5(filename = "data/scatac-seq/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
  data <- load_atac(peaks, scfile = "data/scatac-seq/atac_v1_pbmc_10k_singlecell.csv", fragments = "data/scatac-seq/atac_v1_pbmc_10k_fragments.tsv.gz")
  
  subset_seurat <- add_annotations(data)
  subset_seurat <- normalise_dataset(subset_seurat)
  
  # subset_seurat <- readRDS("subset_seurat.rds")

  output$tss <- renderPlot(TSSPlot(subset_seurat, group.by = 'high.tss') + NoLegend())

  output$nucleosome <- renderPlot(FragmentHistogram(object = subset_seurat, group.by = 'nucleosome_group'))

  output$qcmetrics <- renderPlot(VlnPlot(
    object = subset_seurat,
    features = c('pct_reads_in_peaks', 'peak_region_fragments',
                 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
    pt.size = 0.1,
    ncol = 5
  ))


  # # observe({
  # #   if(input$run.atac){
  # #      seurat$filtered<- filter_dataset(
  # #       data, min.peak = input$min.peaks.frag, max.peak = input$max.peaks.frag,
  # #       pct.blacklist = input$blacklist, pct.reads = input$pct.reads,
  # #       nuc.signal = input$min.nuc.signal, tss.enrich = input$tss.enrichment)
  # #     
  # #   }
  # # })


  subset_seurat <- filter_dataset(subset_seurat)
  subset_seurat <- linear_dim_reduction(subset_seurat)
  output$svd <- renderPlot(DepthCor(subset_seurat))

  subset_seurat <- nonlinear_dim_reduction(subset_seurat)
  output$atac.umap <- renderPlot(DimPlot(object = subset_seurat, label = TRUE) + NoLegend())

  subset_seurat <- get_gene_activities(subset_seurat)

  seurat_rna <- load_rnaseq("data/scatac-seq/pbmc_10k_v3.rds")

  subset_seurat <- integrate_rnaseq(subset_seurat, seurat_rna)

  output$atac.rna.seq <- renderPlot({
    plot1 <- DimPlot(
      object = seurat_rna,
      group.by = 'celltype',
      label = TRUE,
      repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

    plot2 <- DimPlot(
      object = subset_seurat,
      group.by = 'predicted.id',
      label = TRUE,
      repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

    plot1 | plot2
  })
# 
# 
# 
#   da_peaks <- reactive({diff_accessibility(subset_seurat())})
}





























