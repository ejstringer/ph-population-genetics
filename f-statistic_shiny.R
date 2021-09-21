library(shiny)
library(tidyverse)
library(emR)


ui <- fluidPage(
  
  sidebarLayout(
    input
    ,
    
    mainPanel(
      
    )
  )
)

server <- function(input, output, session) {
  
 
  
}
shinyApp(ui = ui, server = server)
