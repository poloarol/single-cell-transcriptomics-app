

if (!require('pacman')) install.packages("pacman")

ui <- fluidPage(
  
  titlePanel("A Shiny app to facilate analysis of single-cell transcriptomics data"),
  
  navbarPage(
    "single-cell transcriptomics",
    tabPanel("Single-Cell Transcriptomics",
             mainPanel(
               
               h2("Quality Control and Cell Selection"),
               
               fluidRow(
                 column(4, 
                        h4("Project Name"),
                        textInput("proj_name", NULL, value = "pbmc3k", width = NULL, placeholder = NULL)
                 ),
                 column(4,
                        h4("Min. Number of Cells"),
                        textInput("min.cells", NULL, value = 3, width = NULL, placeholder = NULL)
                 ),
                 column(4,
                        h4("Min. Number of features"),
                        textInput("min.feats", NULL, value = 200, width = NULL, placeholder = NULL)
                 )
               ),
               
               fluidRow(
                 column(6, "QC Metrics",plotOutput("metrics")),
                 column(6, "Feature-Feature Relationships", plotOutput("features"))
               ),
               
               width = 12
               
               
               # ("")
               
             )),
    tabPanel("Single-Cell Intergration",
             mainPanel()),
    tabPanel("Spatial Transcriptomics Analysis",
             mainPanel()),
    navbarMenu("More",
               "Repository")
  )
  
)