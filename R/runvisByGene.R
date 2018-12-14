######## Shiny APP ########
### visByGene ###

#'shiny interface to Kimbie's CISBP FIMO scan TF bed files
#'@import shiny
#'@export
runvisByGene<-function(){
  myapp = system.file("shinyApps",package="visByGene")
  runApp(myapp)
}