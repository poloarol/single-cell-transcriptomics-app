
source("rna_seq.R")


server <- function(input, output, session) {
  
  rnaseq_analysis(input, output, session)
  
}






