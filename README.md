
# A R Shiny App to ease single-cell RNA transcriptomics analysis


## Aim
-------
The aim of this project is to build a minimalist, but flexible R Shiny application 
to aid wet lab scientist to analyse their single cell data. 

The analysis to be included are;

1. Single cell data clustering
2. Single cell data integration
3. Single cell ATAC-seq analysis
4. Integration of scRNA-seq and scATAC-seq data
5. Integration of cite-Seq data
6. Automatic cell labelling

## Status
1. The complete clustering workflow is complete
2. Comming soon: automatic labelling and doublet removal

## Packages
-----------
1. Seurat
2. Signac (scATAC-seq)
3. DoubletFinder (Removal of doublets)
4. SingleR (Automatic labelling of cells)

## Usage
--------
Currently the program can only be used by closing to the repository and running locally

The next step would be to provide it as a docker image to ease its distribution.
