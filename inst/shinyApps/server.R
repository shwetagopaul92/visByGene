########################
# Functions
########################

visByGene <- function(myeqtl, myfp, mytf, mysymbol){
  require(Gviz)
  require(TFutils)
  require(GenomicRanges)
  require(mongolite)
  mygene = genemodForGviz(sym=mysymbol,resource= EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  mychr = seqnames(mygene)[1]
  mygeneRange = GRanges(mychr, IRanges(start=min(start(mygene)),end=max(end(mygene))))
  mystart = min(start(mygene))
  myend = max(end(mygene))
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = "hg19", chromosome = mychr)
  biomtrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = mychr, symbol=mysymbol, name=mysymbol)
  #baseplot = plotTracks(list(gtrack,itrack,biomtrack), sizes=c(1,1,3))
  
  #FOR TF data 
  tfcoll = mongo(url="mongodb://172.27.105.48",
                 collection = mytf,
                 db = "txregnet"
  )
  mytfquery=paste0('{',
                   '"chr":"',mychr,'",',
                   '"start" : {"$gte":',mystart,'},',
                   '"end" : {"$lte":',myend,'}}')
  restf=tfcoll$find(mytfquery)
  mysubtf=GRanges(restf$chr, IRanges(restf$start, restf$end), mcols=restf)
  #tf = Rsamtools::TabixFile(paste0("/udd/reshg/tbifiles/tabix_all_tf_new/",mytf,".02_sort.bed.gz"))
  #mysub = importFIMO(tf, mygeneRange)
  genome(mysubtf)="hg19"
  atrack1 <- AnnotationTrack(mysubtf, name = mytf)
  dtrack1=DataTrack(data=mysubtf$mcols.pvalue, start=start(mysubtf), end=end(mysubtf),
                    chromosome = mychr, genome = "hg19", name = "TF", background.panel = "#FFFEDB",
                    background.title = "darkblue")
  #tfp = plotTracks(dtrack)
  # FOR FP data 
  fpcoll = mongo(url="mongodb://172.27.105.48",
                 collection = myfp,
                 db = "txregnet"
  )
  myfpquery=paste0('{',
                   '"chr":"',mychr,'",',
                   '"start" : {"$gte":',mystart,'},',
                   '"end" : {"$lte":',myend,'}}')
  resfp=fpcoll$find(myfpquery)
  mysubfp=GRanges(resfp$chr, IRanges(resfp$start, resfp$end), mcols=resfp)
  #tf = Rsamtools::TabixFile(paste0("/udd/reshg/tbifiles/tabix_all_tf_new/",mytf,".02_sort.bed.gz"))
  #mysub = importFIMO(tf, mygeneRange)
  genome(mysubfp)="hg19"
  atrack1 <- AnnotationTrack(mysubfp, name = myfp)
  dtrack2=DataTrack(data=mysubfp$mcols.stat, start=start(mysubfp), end=end(mysubfp),
                    chromosome = mychr, genome = "hg19", name = "FP")
  #fpp = plotTracks(dtrack)
  # FOR eQTL data 
  eqtlcoll = mongo(url="mongodb://172.27.105.48",
                   collection = myeqtl,
                   db = "txregnet"
  )
  mychr2=as.numeric(x=sub(mychr,pattern="chr",replacement=""))
  myeqtlquery=paste0('{',
                     '"chr":',mychr2,',',
                     '"snp_pos" : {"$gte":',mystart,'},',
                     '"snp_pos" : {"$lte":',myend,'}}')
  reseqtl=eqtlcoll$find(myeqtlquery)
  mysubeqtl=GRanges(reseqtl$chr, IRanges(reseqtl$snp_pos, width = 1), mcols=reseqtl)
  #tf = Rsamtools::TabixFile(paste0("/udd/reshg/tbifiles/tabix_all_tf_new/",mytf,".02_sort.bed.gz"))
  #mysub = importFIMO(tf, mygeneRange)
  genome(mysubeqtl)="hg19"
  atrack1 <- AnnotationTrack(mysubeqtl, name = myeqtl)
  dtrack3=DataTrack(data=mysubeqtl$mcols.qvalue, start=start(mysubeqtl), end=end(mysubeqtl),
                    chromosome = mychr, genome = "hg19", name = "eQTL",  background.panel = "#FFFEDB",
                    background.title = "darkblue")
  #eqtlp = plotTracks(dtrack)
  visplot = plotTracks(list(gtrack,biomtrack,dtrack1,dtrack2,dtrack3,itrack), sizes=c(3,7,5,5,5,10))
  
}

plotTF <- function(mytf, mysymbol){
  require(Gviz)
  require(mongolite)
  mygene = TFutils::genemodForGviz(sym=mysymbol,resource= EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  mychr = GenomicRanges::seqnames(mygene)[1]
  mygeneRange = GenomicRanges::GRanges(mychr, IRanges(start=min(start(mygene)),end=max(end(mygene))))
  mystart = min(GenomicRanges::start(mygene))
  myend = max(GenomicRanges::end(mygene))
  tfcoll = mongo(url="mongodb://172.27.105.48",
                 collection = mytf,
                 db = "txregnet"
  )
  myquery=paste0('{',
                 '"chr":"',mychr,'",',
                 '"start" : {"$gte":',mystart,'},',
                 '"end" : {"$lte":',myend,'}}')
  res=tfcoll$find(myquery)
  validate(
    need(res != "", "Please select another transcription factor")
  )
  mysub=GenomicRanges::GRanges(res$chr, IRanges(res$start, res$end), mcols=res)
  #tf = Rsamtools::TabixFile(paste0("/udd/reshg/tbifiles/tabix_all_tf_new/",mytf,".02_sort.bed.gz"))
  #mysub = importFIMO(tf, mygeneRange)
  genome(mysub)="hg19"
  atrack1 <- AnnotationTrack(mysub, name = mytf)
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = "hg19", chromosome = mychr)
  dtrack=DataTrack(data=mysub$mcols.pvalue, start=start(mysub), end=end(mysub),
                   chromosome = mychr, genome = "hg19", name = mytf)
  biomtrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = mychr, symbol=mysymbol, name=mysymbol)
  tfp = plotTracks(list(dtrack, gtrack, biomtrack, itrack),sizes=c(1,1,2,2))
}

plotFP <- function(myfp, mysymbol){
  require(Gviz)
  require(mongolite)
  mygene = TFutils::genemodForGviz(sym=mysymbol,resource= EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  mychr = GenomicRanges::seqnames(mygene)[1]
  mygeneRange = GenomicRanges::GRanges(mychr, IRanges(start=min(start(mygene)),end=max(end(mygene))))
  mystart = min(GenomicRanges::start(mygene))
  myend = max(GenomicRanges::end(mygene))
  fpcoll = mongo(url="mongodb://172.27.105.48",
                 collection = myfp,
                 db = "txregnet"
  )
  myquery=paste0('{',
                 '"chr":"',mychr,'",',
                 '"start" : {"$gte":',mystart,'},',
                 '"end" : {"$lte":',myend,'}}')
  res=fpcoll$find(myquery)
  validate(
    need(res != "", "Please select another Footprint")
  )
  mysub=GenomicRanges::GRanges(res$chr, IRanges(res$start, res$end), mcols=res)
  #tf = Rsamtools::TabixFile(paste0("/udd/reshg/tbifiles/tabix_all_tf_new/",mytf,".02_sort.bed.gz"))
  #mysub = importFIMO(tf, mygeneRange)
  genome(mysub)="hg19"
  atrack1 <- AnnotationTrack(mysub, name = myfp)
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = "hg19", chromosome = mychr)
  dtrack=DataTrack(data=mysub$mcols.stat, start=start(mysub), end=end(mysub),
                   chromosome = mychr, genome = "hg19", name = myfp)
  biomtrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = mychr, symbol=mysymbol, name=mysymbol)
  fpp = plotTracks(list(dtrack, gtrack, biomtrack, itrack),sizes=c(1,1,2,2))
}

ploteQTL <- function(myeqtl, mysymbol){
  require(Gviz)
  require(mongolite)
  mygene = TFutils::genemodForGviz(sym=mysymbol,resource= EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  mychr = GenomicRanges::seqnames(mygene)[1]
  mygeneRange = GenomicRanges::GRanges(mychr, IRanges(start=min(start(mygene)),end=max(end(mygene))))
  mystart = min(GenomicRanges::start(mygene))
  myend = max(GenomicRanges::end(mygene))
  eqtlcoll = mongo(url="mongodb://172.27.105.48",
                   collection = myeqtl,
                   db = "txregnet"
  )
  myquery=paste0('{',
                 '"chr":',17,',',
                 '"snp_pos" : {"$gte":',mystart,'},',
                 '"snp_pos" : {"$lte":',myend,'}}')
  res=eqtlcoll$find(myquery)
  validate(
    need(res != "", "Please select another eQTL")
  )
  mysub=GenomicRanges::GRanges(res$chr, IRanges(res$snp_pos, width = 1), mcols=res)
  #tf = Rsamtools::TabixFile(paste0("/udd/reshg/tbifiles/tabix_all_tf_new/",mytf,".02_sort.bed.gz"))
  #mysub = importFIMO(tf, mygeneRange)
  genome(mysub)="hg19"
  atrack1 <- AnnotationTrack(mysub, name = myeqtl)
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = "hg19", chromosome = mychr)
  dtrack=DataTrack(data=mysub$mcols.qvalue, start=start(mysub), end=end(mysub),
                   chromosome = mychr, genome = "hg19", name = "eqtl")
  biomtrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = mychr, symbol=mysymbol, name=mysymbol)
  eqtlp = plotTracks(list(dtrack, gtrack, biomtrack, itrack), sizes=c(1,1,2,2))
}

plotHS <- function(myhs, mysymbol){
  require(Gviz)
  require(GenomicRanges)
  mygene = TFutils::genemodForGviz(sym=mysymbol,resource= EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  mychr = GenomicRanges::seqnames(mygene)[1]
  mygeneRange = GenomicRanges::GRanges(mychr, IRanges(start=min(start(mygene)),end=max(end(mygene))))
  mystart = min(GenomicRanges::start(mygene))
  myend = max(GenomicRanges::end(mygene))
  hscoll = mongo(url="mongodb://172.27.105.48",
                 collection = myhs,
                 db = "txregnet"
  )
  myquery=paste0('{',
                 '"chrom":"',mychr,'",',
                 '"chromStart" : {"$gte":',mystart,'},',
                 '"chromEnd" : {"$lte":',myend,'}}')
  res=hscoll$find(myquery)
  validate(
    need(res != "", "Please select another Hypersensitivity collection")
  )
  mysub=GenomicRanges::GRanges(res$chrom, IRanges(res$chromStart, res$chromEnd), mcols=res)
  #tf = Rsamtools::TabixFile(paste0("/udd/reshg/tbifiles/tabix_all_tf_new/",mytf,".02_sort.bed.gz"))
  #mysub = importFIMO(tf, mygeneRange)
  genome(mysub)="hg19"
  atrack1 <- AnnotationTrack(mysub, name = myhs)
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = "hg19", chromosome = mychr)
  dtrack=DataTrack(data=mysub$mcols.qValue, start=start(mysub), end=end(mysub),
                   chromosome = mychr, genome = "hg19", name = myhs)
  biomtrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = mychr, symbol=mysymbol, name=mysymbol)
  hsp = plotTracks(list(dtrack, gtrack, biomtrack, itrack), sizes=c(1,1,2,2))
}

ggvisForSymbol = function (sym, resource = EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79, 
                           columnsKept = c("gene_id", "tx_id"), yval = 1, arrmm=1.5, viewtype="transcripts", ...) 
{
  exs = GenomicFeatures::exons(resource, filter = SymbolFilter(sym), columns = columnsKept, 
                               ...)
  if (viewtype == "exons") exs = unique(exs)
  rd = reduce(exs)
  fo = findOverlaps(rd, exs)
  gr = split(subjectHits(fo), queryHits(fo))
  pp = function(n) (seq_len(n)-1)/n
  st = start(exs)
  en = end(exs)
  if (viewtype == "exons") {
    ys = lapply(gr, function(x) pp(length(x)))
    yvs = unlist(ys) #1+(0:(nel-1))/nel
  }
  else if (viewtype == "transcripts") {
    tnms = exs$tx_id
    ft = factor(tnms)
    yvs = (as.numeric(ft)-1)/length(levels(ft))
  }
  else stop("viewtype not %in% c('exons', 'transcripts')")
  newdf = data.frame(st, en, yv = yvs, sym = sym)
  rng = range(exs)
  df = data.frame(range = c(start(rng), end(rng)), yval = rep(yval,2)) 
  strn = as.character(strand(exs)[1])
  ardir = ifelse(strn=="+", "last", "first")
  pl = ggplot(df, aes(x = range, y = yval)) + 
    geom_segment(aes(x = st, y = yv, xend = en, yend = yv, colour = sym),data = newdf, arrow=arrow(ends=ardir, length=unit(arrmm, "mm")))
  #ggplotly(pl)
  pl + xlab(as.character(seqnames(exs)[1]))
}



########################
# Define the server logic
########################

shinyServer(function(input, output, session) {
  
  #observeEvent(input$goButton, {
   # output$eqtlplot = renderPlot(ploteQTL(input$eQTL, input$geneName))
  #})
  #observeEvent(input$addButton, {
   # output$eqtlplot = renderPlot(ploteQTL2(input$eQTL, input$geneName))
  #})
  
  output$tfplot = renderPlot(plotTF(input$transcriptionFactor, input$geneName))
  
  output$fpplot = renderPlot(plotFP(input$fp, input$geneName))
  
  output$eqtlplot = renderPlot(ploteQTL(input$eQTL, input$geneName))
  
  output$tfgenePlot = renderPlot(enc690ByFactor(input$encodeTF, input$geneName))
  
  output$hsplot = renderPlot(plotHS(input$hs, input$geneName))
  
  output$visByGene = renderPlot(visByGene(input$eQTL, input$fp, input$transcriptionFactor, input$geneName))
  
  output$mytable2 = DT::renderDataTable(as.data.frame(metadata_tf), options=list(scrollX=TRUE,pageLength=25))
  
  
})