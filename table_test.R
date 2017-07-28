library(shiny)

ui <- fluidPage("RenderTable Example",
                tableOutput("bTable") #table of information on beds (beds)
)


serverShouldWork <- function(input, output) {
  output$bTable <- renderTable({iris}, caption=paste("test"),
                               caption.placement = getOption("xtable.caption.placement", "top"),
                               caption.width = getOption("xtable.caption.width", NULL)
  )
}

serverDoesWork <- function(input, output) {
  output$bTable <- renderTable({iris}, caption="test")
}



app1 <- shinyApp(ui, serverDoesWork)
#app2 <- shinyApp(ui, serverShouldWork)

runApp(app1)
#runApp(app2)