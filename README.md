
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
5. Multimodal analysis using cite-seq
6. Automatic cell labelling


Supported data formats
----------------------
1. 10X Cell Ranger (.HDF5)
2. loom
3. rds
4. AnnData (.h5ad)
5. .mtx, name your files as; features_*.tsv, matrix_*.mtrx and genes_*.tsv


## Status
---------
- [X] scRNA-seq clustering workflow
- [X] Automatic cell labelling
- [ ] scATAC-seq workflow
- [ ] Doublet removal
- [ ] scRNA-seq data integration
- [ ] scRNA-seq + scATAC seq data integration
- [ ] Multimodal analysis using cite-seq

## Packages
-----------
1. Seurat
2. Signac (scATAC-seq)
3. DoubletFinder (Removal of doublets)
4. SingleR (Automatic labelling of cells)
5. celldex (Autmatic labelling of cells)

## Usage
--------
Currently the program can only be used by cloning the repository and running locally

The next step would be to provide it as a docker image to ease its distribution.
