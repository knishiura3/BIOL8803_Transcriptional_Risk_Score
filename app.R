library(shiny)
library(shinythemes)
library(reticulate)

#App functionality to include: user input of GWAS data

ui <- fluidPage(theme = shinytheme("superhero"),
  
  # App title ----
  titlePanel("Colocalization App"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      #Describe the app here.
      p("This app takes in RNAseq data and GWAS to generate a plot of colocalized variants"),
      
      # Input: Slider for the number of bins ----
      #sliderInput(inputId = "bins",
                  #label = "Number of bins:",
                 # min = 1,
                 # max = 50,
                 # value = 30),
      
      #Input: allow user text entry.
      textInput("gwasID", "Enter a GWAS ID:", value = "", width = NULL, placeholder = NULL),
      
      #Input: Allow user to select one of the eQTL mapping databases
      #Find a light way to include the blood database in the app.
      selectInput(
        "choice",
        "Pick a eQTL mapping database",
        c("Blood eQTL", "TIGER"),
        selected = NULL,
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL
      ),
      
      tags$a(href="https://genome.ucsc.edu/", "Find a GWAS ID Here!")
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      #Output: Score and Histogram.
      #p("Here's where the results will go!")
      
      #Output different aspects with tabs
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput(outputId = "distPlot")),
                  tabPanel("Score Results", p("TRS is 73"))
      ),
                  
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot")
      
    )
  )
)

server <- function(input, output) {
  
  
  #Using renderPlot makes the app reactive to user input. When inputs change
  #The outputs are re-executed automatically.
  output$distPlot <- renderPlot({
    
    x    <- faithful$waiting
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(x, breaks = bins, col = "#75AADB", border = "white",
         xlab = "Waiting time to next eruption (in mins)",
         main = "Histogram of waiting times")

    #Run python code within Shiny app.
    
    #Build th
    py_run_script()
    
    
  })
  
}
  



#Can also be run with runApp("TRSApp") when TRSApp is a directory in the 
#working directory with this app.R file within it.
shinyApp(ui = ui, server = server)

#To exit out of the app, press escape or click the stop sign in the upper right.