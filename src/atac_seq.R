
if (!require('pacman')) install.packages("pacman")

# Load contributed packages with pacman
pacman::p_load(pacman, Seurat, SeuratObject, Signac, 
               tidyverse, shiny, DT, shinyFiles, shinyWidgets,
               Rsamtools, patchwork)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19', 'EnsDb.Hsapiens.v75'))
# BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
# BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))

require(GenomeInfoDb)
require(EnsDb.Hsapiens.v75)
# require(EnsDb.Hsapiens.v79)
# require(EnsDb.Mmusculus.v79)

load_data <- function(matrix, scfile, fragments, reference, min.cells, min.features){
  
  metadata <- read.csv(
    file = scfile,
    header = TRUE,
    row.names = 1
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = reference,
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

add_annotations <- function(data){
  # if(ensdb == "BSgenome.Hsapiens.UCSC.hg19" || ensdb == "EnsDb.Hsapiens.v75"){
  #   annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
  # }else if(ensdb == "BSgenome.Hsapiens.UCSC.hg38" || ensdb == "EnsDb.Hsapiens.v86"){
  #   annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  # }else{
  #   annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v79)
  # }
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
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
  obj
}


atacseq_analysis <- function(input, output, session){
  seurat <- reactiveValues(filtered = NULL)
  
  seurat_obj <- eventReactive(
    eventExpr = {
      input$atac.peaks
    }, {
      counts <- Read10X_h5(filename = input$atac.peaks$datapath)
      load_data(couts, 
                input$atac.cellranger$datapath, 
                input$atac.fragments$datapath,
                input$min.atac.cells,
                input$min.feat.cells)
    }
  )
  
  seurat_obj <- reactive({add_annotations(seurat_obj())})

  # compute nucleosome signal score per cell
  seurat_obj <- reactive({NucleosomeSignal(object = seurat_obj())})

  # compute TSS enrichment score per cell
  seurat_obj <- reactive({TSSenrichment(object = seurat_obj(), fast = FALSE)})

  # add blacklist ratio and fraction of reads in peaks
  seurat_obj <- reactive({
    pct <- seurat_obj()$peak_region_fragments / seurat_obj()$passed_filters * 100
    seurat_obj$pct_reads_in_peak <- pct
    })

  seurat_obj <-  reactive({
    data <- seurat_obj()
    ratio <- data$blacklist_region_fragments / data$peak_region_fragments
    data$blacklist_ratio <- ratio
    data
    })

  seurat_obj <- reactive({
    tss <- ifelse(seurat_obj()$TSS.enrichment > input$tss.enrichment, 'High', 'Low')
    seurat_obj$high.tss <- tss
    })

  output$tss <- renderPlot(TSSPlot(seurat_obj(), group.by = 'high.tss') + NoLegend())

  seurat_obj <- reactive({
    minimum <- paste("NS >", as.character(input$min.nuc.signal))
    maximum <- paste("NS < ", as.character(input$min.nuc.signal))
    data <- seurat_obj()
    nuc <- ifelse(data$nucleosome_signal >
             input$min.nuc.signal, minimum, maximum)
    data$nucleosome_group <- nuc
    data
    })

  output$nuc <- renderPlot(FragmentHistogram(object = seurat_obj(), group.by = 'nucleosome_group'))

  output$qcmetrics <- renderPlot(VlnPlot(
    object = seurat_obj(),
    features = c('pct_reads_in_peaks', 'peak_region_fragments',
                 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
    pt.size = 0.1,
    ncol = 5
  ))


  # observe({
  #   if(input$run.atac){
  #      seurat$filtered<- filter_dataset(
  #       data, min.peak = input$min.peaks.frag, max.peak = input$max.peaks.frag,
  #       pct.blacklist = input$blacklist, pct.reads = input$pct.reads,
  #       nuc.signal = input$min.nuc.signal, tss.enrich = input$tss.enrichment)
  #     
  #   }
  # })
  
}





























