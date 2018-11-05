########################
# Load the data
########################

data(named_tf)
data(named_metadata_tf)
data(encode690)


########################
# Define UI for the application
########################

shinyUI(fluidPage(
  
  titlePanel("visByGene"),
  sidebarLayout(position = "left",
                sidebarPanel(width=3,
                             selectInput("eQTL", "Select eQTL Collection", eQTLdata),
                             selectInput("fp", "Select FootPrint Collection", fpdata),
                             selectInput("hs", "Select Hypersensitivity Collection", hsdata),
                             selectInput("transcriptionFactor", "Select FIMO Transcription Factor Collection", tfdata),
                             textInput("geneName", "Enter Gene of interest", value="ORMDL3"),
                             numericInput("windowSize", "Window Size", 1000),
                             #actionButton("goButton", "Go!"),
                             #actionButton("addButton", "Add!"),
                             textInput("downloadName", "Download Name"),
                             downloadButton("downloadData", "Download Data")
                ),
                
                
                mainPanel(
                  tabsetPanel(id="inTabset",
                              #tabPanel("Scored Motifs in Transcribed Region", value="panel1", DT::dataTableOutput("mytable1")),
                              tabPanel("Gene Model", value="panel3", plotlyOutput("ggPlot")),
                              #tabPanel("Encode TF Model", value="panel4", plotOutput("tfgenePlot")),
                              tabPanel("eQTL Model", value="panel1", plotOutput("eqtlplot")),
                              tabPanel("Footprint Model", value="panel2", plotOutput("fpplot")),
                              tabPanel("DNAse Hypersensitivity Model", value="panel3", plotOutput("hsplot")),
                              tabPanel("FIMO TF Model", value="panel4", plotOutput("tfplot")),
                              tabPanel("visByGene", value="panel5", plotOutput("visByGene")),
                              #tabPanel("Plot", value="panel5", 
                                       #fluidRow(
                                        #plotOutput("fpplot"),
                                        #plotOutput("eqtlplot")
                              #)),
                              tabPanel("Metadata",value="panel6", DT::dataTableOutput("mytable2")) 
                 
                  
                  
                )
  )
)
))