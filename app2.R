

library(shiny)
library(shinycssloaders)

# user interface

ui=shinyUI(fluidPage(
  sidebarLayout(
    sidebarPanel(
      p('Prepare 3 input files:'),
      p('1. H3K27ac ChIP-seq signal RPKM'),
      p('2. H3K27ac CHIP-seq peaks in bed'),
      p('3. RNA-seq Gencode log2(RPKM)'),
     
      # 1. ChIp-seq singal RPKM 
      fileInput('file1', 'ChIP-seq RPKM signal file',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                )),
      
      # 2. ChiP-seq bed
      fileInput('file2', 'ChIP-seq peaks in bed format',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                )),
      fileInput('file3', 'RNA-seq Gencode log2(RPKM)',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                ))),
      mainPanel(
        h3('ChIP-seq signal data'),
        dataTableOutput('chipSignal'),
        h3('RNA-seq signal data'),
        dataTableOutput('rnaSignal'))
      
    
    
  )
  
))

#### server script

server=function(input,output){
  ress=reactive({
    
    #load chip signal data
    inFile.chip=input$file1
    chipSignal.load=read.table(inFile.chip$datapath, sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = "region")
    
    # load chip peak data bed file
    ## write this code later
    
    # load rna signal data
    infile.rna=input$file3
    rnaSignal.load=read.table(infile.rna, sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = "Gene")
    
    # output$chipSignal=renderDataTable({head(chipSignal.load)})
    # output$rnaSignal=renderDataTable({head(rnaSignal.load)})
    })
  
  output$chipSignal=renderDataTable(ress()$chipSignal.load, options = list(pageLength = 5))
  output$rnaSignal=renderDataTable(ress()$rnaSignal.load, options=list(pageLength=5))

    
}



# run app
shinyApp(ui, server)










