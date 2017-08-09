library(shiny)
library(shinycssloaders)

### Version - 2

# Only Take ChiP signals and RNA signals as input
# Prepare peaks internally from ChIP signals data

#options(error=function() dump.frames(to.file=TRUE, dumpto ="log.txt"))

ui<-shinyUI(fluidPage(
  titlePanel("inTAD - Enhancer Gene Association within TAD regions"),
  
  sidebarLayout(
    sidebarPanel(
      h5('Prepare 2 input files:'),
      p('1. H3K27ac ChIP-seq RPKM normalized signals within peaks'),
      p('Example: ', a(href = 'https://raw.githubusercontent.com/venuthatikonda/inTAD-shiny/master/chipSignal.txt', 'ChIP_Signal.txt')),
      # p('2. H3K27ac ChIP-seq peaks in bed format'),
      # p('Example: ', a(href = 'https://raw.githubusercontent.com/venuthatikonda/inTAD-shiny/master/chip-cordi.txt', 'ChIP_coordinates.txt')),
      p('3. RNA-seq log2(RPKM) normalized counts for Gencode annotaions'),
      p('Example: ', a(href = 'https://raw.githubusercontent.com/venuthatikonda/inTAD-shiny/master/RNA_seq_signal.txt', 'RNAseq_Signal.txt')),
      
      #p('Upload a file here'),
      fileInput('file1', 'H3K27ac CHIP-seq signals',
                accept = c('text/csv',
                           'text/comma-separated-values',
                           'text/tab-separated-values',
                           'text/plain',
                           '.csv',
                           '.tsv',
                           '.txt'
                )
      ),
      # fileInput('file2', 'ChiP-seq bed file', 
      #           accept = c('text/csv',
      #                      'text/comma-separated-values',
      #                      'text/tab-separated-values',
      #                      'text/plain',
      #                      '.csv',
      #                      '.tsv',
      #                      '.txt'
      #           )
      #),
      fileInput('file3', 'RNA-seq singal file',
                accept = c('text/csv',
                           'text/comma-separated-values',
                           'text/tab-separated-values',
                           'text/plain',
                           '.csv',
                           '.tsv',
                           '.txt'
                )
      ),
      
      actionButton('analyse', 'Perform analysis'),
      h4('Downstream analysis'),
      p('1. Download unfiltered data'),
      downloadButton('downloadData', 'Download'),
      p('2. Perform GO Enrichment analysis with enhancer asociated genes'),
      actionButton('goEnrich', 'Perform GO enrichment'),
      p('3. Download GO enrichment results'),
      downloadButton('downloadGO', 'Download GO')
    ),
    mainPanel(
      h2('Top ranked Enhancer-Gene correlations'),
      withSpinner(dataTableOutput('corResults')),
      
      # h2('Top GO terms enriched with Enhancer associated genes'),
      # withSpinner(plotOutput('goPlot')),
      
      h2('Enriched GO terms explained'),
      withSpinner(dataTableOutput('goResults')),
      
      h2("Top significantly enriched GO terms"),
      withSpinner(plotOutput('BarGO', width = "100%", height = "400px"))
      # h4('Few lines from ChIP-seq Bed file'),
      # dataTableOutput('chipBedLines'),
      # h4('Few lines from RNA-seq singal file'),
      # dataTableOutput('rnaSignalLines')
    )
  )
  
))

### server logic


server<-shinyServer(function(input, output){
  library(InTAD)
  library(GenomicRanges)
  library(rtracklayer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  # Gencode annotations
  gencode.txs=readRDS("/Users/venu/Desktop/work/pipelines/inTAD/inTAD-shiny/gencode.transcript.v19.rds")
  
  # ChIP-seq signal file
  ChipSignalPeaks<-eventReactive(input$analyse, {
    
    chipSingalFile<-input$file1
    chipSignalFileLoad=read.delim(chipSingalFile$datapath, header = TRUE)
    
    #prepare peaks
    
    chipPEAKS=makeGRangesFromDataFrame(chipSignalFileLoad, seqnames.field = "chr", start.field = "start", end.field = "end")
    chipPEAKSdf=subset(chipSignalFileLoad, select=c(chr,start,end))
    row.names(chipSignalFileLoad)<-as.character(chipPEAKS)
    chipSignalFileLoad=subset(chipSignalFileLoad, select=-c(chr,start,end))
    list(chipSignalFileLoad, chipPEAKSdf)
    # in above list [[1]] is prepared chip signal
    # and [[2]] is a dataframe of peaks - convert this into GRanges object and use with inTAD functions
    
  })
  
  # # ChIP-seq Bed fileas
  # dfChipBed<-eventReactive(input$analyse, {
  #   chipBedFile=input$file2
  #   chipBedFileLoad=read.table(chipBedFile$datapath, header = FALSE, sep="\t", stringsAsFactors = FALSE)
  #   chipBedFileLoad
  # })
  
  # RNA-seq signal file Gencode annotations
  dfRnaSignal=eventReactive(input$analyse, {
    rnaSignalFile=input$file3
    rnaSignalFileLoad=read.table(rnaSignalFile$datapath, sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = "Gene")
    rnaSignalFileLoad
  })
  
  ## inTAD analysis
  
  intadCor=eventReactive(input$analyse, {
    inTadSig=newSigInTAD(ChipSignalPeaks()[[1]], makeGRangesFromDataFrame(ChipSignalPeaks()[[2]], seqnames.field = "chr", start.field = "start", end.field = "end"), dfRnaSignal(), gencode.txs)
    inTadSig=filterGeneExpr(inTadSig, geneType = "protein_coding")
    inTadSig=combineInTAD(inTadSig, tadGR)
    
    corData<-findCorrelation(inTadSig)
    res=cbind(corData, BHadj=p.adjust(corData$pvalue, n=length(corData$pvalue)))
    res
    
  })
  
  # GO enrichment analysis with enhancer associated genes
  # NOTE ##############
  # as the example data is too low, I'm performing GO enrichment with raw RNA-seq gene data for demonstration purpose
  
  intadGO=eventReactive(input$goEnrich, {
    ensembl=unique(gsub("\\..*", "", row.names(dfRnaSignal())))
    uniprot=bitr(ensembl, fromType="ENSEMBL", toType="UNIPROT", OrgDb=org.Hs.eg.db)
    ego=enrichGO(gene=uniprot$UNIPROT, keytype='UNIPROT', OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05)
    ego
  })
  
  # intad result output
  output$corResults=renderDataTable( 
    
    head(intadCor())
    
  )
  
  # GO enrichment result output
  output$goResults=renderDataTable(
    head(subset(intadGO()@result, select=-c(geneID)))
  )
  
  # GO enrichment plot
  
  output$BarGO=renderPlot({
    barplot(intadGO(), showCategory=5)
  })
  
  # intad result download
  output$downloadData=downloadHandler(
    filename = function(){"inTAD_results.txt"},
    content = function(file){
      write.table(intadCor(), file, sep="\t", quote=FALSE, row.names = FALSE)
    }
  )
  # GO enrichment result download
  
  output$downloadGO=downloadHandler(
    filename = function(){"inTAD_GOenrichments.txt"},
    content=function(fileGO){
      write.table(intadGO(), fileGO, sep="\t", quote = FALSE, row.names = FALSE)
    }
  )
})

### RUN

shinyApp(ui, server)
