
library(shiny)
library(shinycssloaders)


hs=data.frame(row.names = c('Homo sapiens'))
hs$code=c('hsa')
hs$db=c("org.Hs.eg.db")
hs=as.matrix(hs)

# UI
ui <- shinyUI(fluidPage(
  
  # Title
  titlePanel("inTAD"),
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      p('Prepare 3 input files'),
      p('1. H3K27ac ChIP-seq peaks with RPKM normalized signal'),
      p('2. H3K27ac ChIP-seq peaks'),
      p('3. Gencode RNA-seq log2(RPKM)'),
      
      selectInput("org", "Organism:", 
                  choices=row.names(hs)),
      
      # ChIP-seq signals upload
      fileInput('file1', 'H3K27ac ChIP-seq signals file',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                )
      ),
      # chip-seq bed file
      fileInput('file2', 'H3K27ac ChIP-seq bed file',
                accept=c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                )),
      
      # RNA-seq file
      fileInput('file3', 'RNA-seq file',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                )
      ),
      
      p('Select correlation cut-off'),
      numericInput('corValue', 'Correlation value', 0.6, min=0, max=1, step=0.1),

      p('Select p-value cut-off'),
      numericInput('pVal', 'P-Value cut-off', 0.05, min = 0, max=1, step=0.01),
      
      actionButton("analysis", "Perform Analysis")
      
    ),
    
    
    # Show some sample results
    mainPanel(
      # withSpinner(plotOutput("distPlot",height = '700px')),
      # p(),
      # verbatimTextOutput("pv"),
     
     h4("Top 10 highly correlated and significant Enhancer - gene associations"),
     dataTableOutput('tablecorr')
     
     #column(12,withSpinner(dataTableOutput('tablecorr'),proxy.height = '80px')),
     
      #column(12,dataTableOutput('table'))
     )
)))

server <- shinyServer(function(input, output){
  library(InTAD)
  library(rtracklayer)
  
  options(shiny.maxRequestSize = 9*1024^2)
  showmax=2000
  
  # receive files
  
  ress=reactive({
    
    # ChIP peaks RPKM file
    chipSignal=input$file1
    chipSignal.load=read.table(chipSignal$datapath, header=TRUE, row.names="region")
    
    # ChIP bed file
    chipBed=input$file2
    chipBed.load=import.bed(chipBed$datapath)
    
    # RNA log2(RPKM) file
    rnaSignal=input$file3
    rnaSignal.load=read.table(rnaSignal$datapath, header = TRUE, row.names = "Gene")

    #make sure rownames are same as GR object (0 or 1 based index)    
    rownames(chipSignal.load) <- as.character(chipBed.load)
    
    # load Gencode annotations from github repo - only transcripts
    gencode.txs=readRDS("gencode.transcript.v19.rds")
    inTADSig=newSigInTAD(chipSignal.load, chipBed.load, rnaSignal.load, gencode.txs)
    inTADSig=filterGeneExpr(inTADSig, geneType = "protein_coding")
    inTADSig <- combineInTAD(inTADSig, tadGR)
    
    corData<-findCorrelation(inTADSig)
    
    output$tablecorr=renderDataTable({head(rnaSignal.load)})
    
  })
  
  # output$tablecorr=renderDataTable({
  #   corData
  #   })
})

# Run shiny app
shinyApp(ui = ui, server = server)











