
source("src/rna_seq.R")
source("src/atac_seq.R")


server <- function(input, output, session) {
  
  rnaseq_analysis(input, output, session)
  atacseq_analysis(input, output, session)
  
}






