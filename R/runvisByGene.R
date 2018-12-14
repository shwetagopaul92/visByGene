######## Shiny APP ########
### visByGene ###

#'shiny interface to eQTL,HS,FP,TF files
#'@import shiny
#'@export
runvisByGene<-function(){
  myapp = system.file("shinyApps",package="visByGene")
  runApp(myapp)
}