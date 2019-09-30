# visByGene
Visualization can often be the key to discovering patterns and distilling the massive amount of biological information that rise from the rapid growth of high throughput sequencing experiments. The R package visByGene is a preliminary attack on synthesizing available experiment and annotation data to understand mechanisms of disease using R/Bioconductor (www.bioconductor.org) interfaces to high-resolution genomic data assembled in MongoDB, Bioconductor annotation packages, and tabix indexed files. We developed visByGene to utilize R/Bioconductor packages such as ggplot2, plotly, ensembldb, AnnotationHub, mongolite, GenomicScores, and others to provide users with specialized interactive shiny based visualizations. These visualizations allow users to easily query hundreds of Gb of data in thousands of data resources such as eQTL studies, DNAse I hypersensitivity regions, and transcription factor binding site data and visualize the results as annotated gene models. These queries are gene based with a scalable window parameter. Additional annotations can be overlayed onto the gene model by selecting new experimental results, and the combined results can be easily downloaded as a csv for later use. Further work is being done to utilize the hover text and plot linking features of plotly graphics to enrich the user experience.

## To get started with our tool from within Channing Division of Network Medicine. 

The tool would soon be available to use from any network. Currently it supports mongodb running within a cluster in Channing Network alone. 

- Install "visByGene"

From Github :

- Clone the github repository
> cd visByGene/

- Fire up R

```{
  > library(devtools)
  > build(".")
  > install(".")
  > library(visBygene)
  > runvisByGene()
```


<img src="vignettes/www/visByGene.png" width="850" height="600">



