########################
# Load the data
########################

data(named_tf)
data(named_metadata_tf)
data(encode690)
eQTLdata = system.file("shinyApps/eQTL.csv",package = "visByGene")
fpdata = system.file("shinyApps/FP.csv",package = "visByGene")
hsdata = system.file("shinyApps/HS.csv",package = "visByGene")
tfdata = system.file("shinyApps/TF.csv",package = "visByGene")

########################
# Define UI for the application
########################

shinyUI(fluidPage(
  
  titlePanel("visByGene"),
  sidebarLayout(position = "left",
                sidebarPanel(width=3,
                             selectInput("eQTL", "Select eQTL Collection", eQTLdata),
                             #actionButton("goButton", "Go!"),
                             #actionButton("addButton", "Add another eQTL"),
                             selectInput("fp", "Select FootPrint Collection", fpdata),
                             #selectInput("hs", "Select Hypersensitivity Collection", hsdata),
                             selectInput("transcriptionFactor", "Select FIMO Transcription Factor Collection", tfdata),
                             selectInput("encodeTF", "Select Encode Transcription Factor",encode690$target ),
                             textInput("geneName", "Enter Gene of interest", value="ORMDL3")
                              
                ),
                
                
                mainPanel(
                  tabsetPanel(id="inTabset",
                              #tabPanel("Gene Model", value="panel1", plotlyOutput("ggPlot")),
                              tabPanel("eQTL Model", value="panel2", plotOutput("eqtlplot")),
                              tabPanel("Footprint Model", value="panel3", plotOutput("fpplot")),
                              #tabPanel("DNAse Hypersensitivity Model", value="panel3", plotOutput("hsplot")),
                              tabPanel("FIMO TF Model", value="panel4", plotOutput("tfplot")),
                              tabPanel("Encode TF Model", value="panel5", plotOutput("tfgenePlot")),
                              #tabPanel("CEBPB binding nead ORMDL3", value="panel5", plotlyOutput("plotlySymbol")),
                              #tabPanel("ggplot", value="panel9", plotlyOutput("ggplot")),
                              tabPanel("visByGene", value="panel6", plotOutput("visByGene")),
                              tabPanel("TF Metadata",value="panel7", DT::dataTableOutput("mytable2"))                  )
                )
  )
)
)