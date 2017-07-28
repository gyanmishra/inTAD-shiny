library(shiny)

ui=fluidPage("Render mtcars", 
             dataTableOutput('bTable'))


serverScript=function(input, output){
  
  file1<-read.table("chip-signal.txt", header = TRUE, stringsAsFactors = FALSE)
  
  output$bTable<-renderDataTable({head(file1)})
}

app<-shinyApp(ui, serverScript)

runApp(app)