######## Shiny APP ########
### visByGene ###

#'shiny interface for visualizing datasets
#'@import shiny
#'@param host character string with inet ip address
#'@export
runvisByGene<-function(host){
  myapp = system.file("shinyApps",package="visByGene")
  runApp(myapp, host=host)
}
