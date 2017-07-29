library(shiny)
library(shinycssloaders)

#options(error=function() dump.frames(to.file=TRUE, dumpto ="log.txt"))

ui<-shinyUI(fluidPage(
  titlePanel("inTAD - Enhancer Gene Association"),
  
  sidebarLayout(
    sidebarPanel(
      p('Prepare 3 input files:'),
      p('1. H3K27ac ChIP-seq signals 2. ChIP-seq Bed 3. RNA-seq signals'),
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
      fileInput('file2', 'ChiP-seq bed file', 
                accept = c('text/csv',
                           'text/comma-separated-values',
                           'text/tab-separated-values',
                           'text/plain',
                           '.csv',
                           '.tsv',
                           '.txt'
                )
      ),
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
      
      actionButton('analyse', 'Perform analysis')
      
    ),
    mainPanel(
      h4('Top ranked Enhancer-Gene correlations'),
      dataTableOutput('corResults')
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
  
  # Gencode annotations
  gencode.txs=readRDS("/Users/venu/Desktop/work/pipelines/inTAD/inTAD-shiny/gencode.transcript.v19.rds")
  
  # ChIP-seq signal file
  dfChipSignal<-eventReactive(input$analyse, {
    
    chipSingalFile<-input$file1
    chipSignalFileLoad=read.table(chipSingalFile$datapath, header = TRUE, sep="\t", stringsAsFactors = FALSE, row.names = "region")
    chipSignalFileLoad
  })
  
  # ChIP-seq Bed fileas
  dfChipBed<-eventReactive(input$analyse, {
    chipBedFile=input$file2
    chipBedFileLoad=read.table(chipBedFile$datapath, header = FALSE, sep="\t", stringsAsFactors = FALSE)
    chipBedFileLoad
  })
  
  # RNA-seq signal file Gencode annotations
  dfRnaSignal=eventReactive(input$analyse, {
    rnaSignalFile=input$file3
    rnaSignalFileLoad=read.table(rnaSignalFile$datapath, sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = "Gene")
    rnaSignalFileLoad
  })
  
  ## inTAD analysis
  
  intadCor=eventReactive(input$analyse, {
    inTadSig=newSigInTAD(dfChipSignal(), makeGRangesFromDataFrame(dfChipBed(), seqnames.field = "V1", start.field = "V2", end.field = "V3"), dfRnaSignal(), gencode.txs)
    inTadSig=filterGeneExpr(inTadSig, geneType = "protein_coding")
    inTadSig=combineInTAD(inTadSig, tadGR)
    
    corData<-findCorrelation(inTadSig)
    res=cbind(corData, BHadj=p.adjust(corData$pvalue, n=length(corData$pvalue)))
    res
    
  })
  
  output$corResults=renderDataTable( 
    
    intadCor()
    
    
  )
  
})

### RUN

shinyApp(ui, server)
