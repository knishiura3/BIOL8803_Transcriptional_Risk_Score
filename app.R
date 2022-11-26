#options(shiny.port = 8080)
library(shiny)

#App functionality to include: user input of GWAS data

ui <- fluidPage(
  
  # App title ----
  titlePanel("Coloc App"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      #sliderInput(inputId = "bins",
       #           label = "Number of bins:",
        #          min = 1,
         #         max = 50,
          #        value = 30),
      
      #Input: allow user text entry.
      textInput("OpenGWAS_study", "Enter OpenGWAS ID & press submit:", width = NULL, placeholder = NULL),
      actionButton("submit_gwasid", label = "Submit")
      
      #Input: Allow user to select one of the eQTL mapping databases
      #selectInput(
        #"choice",
        #"Pick a eQTL mapping database",
        #c("TIGER", "Blood eQTL"),
        #selected = NULL,
        #multiple = FALSE,
        #selectize = TRUE,
        #width = NULL,
        #size = NULL
      #)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      #plotOutput(outputId = "distPlot")
      
    )
  )
)

server <- function(input, output) {
  
  observeEvent(input$submit_gwasid, {
    gwas_dataset <<- input$OpenGWAS_study
    #source("blah.R")
    source("eQTL_query_parquet.r")
  })
  #output$distPlot <- renderPlot({
    
    #x    <- faithful$waiting
    #bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    #hist(x, breaks = bins, col = "#75AADB", border = "white",
     #    xlab = "Waiting time to next eruption (in mins)",
      #   main = "Histogram of waiting times")
    
  #})
  
}
  
  




#can be run with runApp("colocApp") when colocApp is the directory.
shinyApp(ui = ui, server = server)