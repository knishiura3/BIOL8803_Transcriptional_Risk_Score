library(shiny)
library(shinythemes)
#library(reticulate)


#App functionality to include: user input of GWAS ID and eQTL browser; Output plot of 

ui <- fluidPage(theme = shinytheme("superhero"),
  
  # App title ----
  titlePanel("Colocalization App"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      #Describe the app here.
      p("This web app takes in eQTL data and GWAS hits to generate a plot of LD variants."),
      #("This app takes in eQTL data and GWAS hits to generate a plot of colocalized variants."),
      
      
      #Input: Allow user to select one of the eQTL mapping databases
      #Find a light way to include the blood database in the app.
      selectInput(
        "choice",
        "Pick a eQTL mapping database",
        c("Blood eQTL", "Second PLACE", "THIRD PLACE"),
        selected = NULL,
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL
      ),
      
      #Input: allow user text entry and provide link to GWAS IDs.
      textInput("gwasID", "Enter a GWAS ID:", value = "ieu-b-30", width = NULL, placeholder = "ex: ieu-b-30"),
      tags$a(href="https://gwas.mrcieu.ac.uk/datasets/?gwas_id__icontains=ukb-b&sort=-gwas_id&page=73", "Find a GWAS ID at this link!")
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      #Output different aspects with tabs
      tabsetPanel(type = "tabs",
                  tabPanel("Results", 
                           textOutput("resultHeader"),
                           imageOutput("peaks"),
                           #uiOutput("peaks"),
                           
                           #Offer method to download a readout of the results.
                           # Button
                           downloadButton("downloadData", "Download Top eQTLs (once plots appear)")
                  ),
                      
                  # tabPanel("About", p("This app reads in "),
                  #          p("The flow chart below describes the backend algorithm.")
                  #          #INCLUDE FLOW CHART HERE
                  #          #img(src = "")
                  # ),
                  tabPanel("Datasets",
                           p("The primary eQTL database is TIGAR's blood eQTL database, and that download can be found along with more supplementary info on our Github page."),
                           p("GWAS IDs are queried from the OpenGWASProject servers, which is a curated set of data with similar organization and quality.")
                           ),
                  tabPanel("Contact Us", 
                           p("This R Shiny Web App was developed by a team of Georgia Institute of Technology students: Andy Chea, Colin Naughton, Kenji Nishiura, and Jasmyn Pellebon. You can reach us at the project's Github page"),
                           tags$a(href="https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score", "Project GitHub")
                           )
      )
      
      
      #Display plot of LD peaks in the results tab.
      #imageOutput("peaks")
      
    )
  )
)

server <- function(input, output) {
  
  
  #Using renderPlot() and similar functions makes the app reactive to user input. When inputs change
  #The outputs are re-executed automatically.
  
  #Build the eQTL database and then query it based on user input.Formerly .py: py_run_script()
  
  #Run R script by Colin and Kenji, using the text input to query the GWAS database.
  #choice variable contains the eQTL choice.
  #gwasID contains the ID to query.
  output$resultHeader <- renderText({paste("Colocalized results of ", input$gwasID, sep = "")})
  #output$resultHeader <- renderText({"Hello!"})
  

  
  #Output the LD plot produced as an image.
  output$peaks <- renderImage({
    
    #Save gwasID as a reactive variable so only it is recalculated within the algorithm.
    gwasInput <- reactive({
      input$gwasID
    })
    
    
    #KENJI/COLIN Code below
    #Reference LD Panels
    #   need to set this to wherever you want the LD panel stored
    dir_ld <- "/ld"
    
    # Copied from the documentation at https://mrcieu.github.io/gwasglue/index.html
    #
    # Updated 1000 genomes LD reference panels (multiple populations):
    # http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz (1.5GB)
    #   <contents>/
    #   âââ AFR.bed
    #   âââ AFR.bim
    #   âââ AFR.fam
    #   âââ AMR.bed
    #   âââ AMR.bim
    #   âââ AMR.fam
    #   âââ EAS.bed
    #   âââ EAS.bim
    #   âââ EAS.fam
    #   âââ EUR.bed
    #   âââ EUR.bim
    #   âââ EUR.fam
    #   âââ SAS.bed
    #   âââ SAS.bim
    #   âââ SAS.fam
    
    
    # declare global variables
    gwas_dataset <<- gwasInput()
    dummy_dataset <<- "ieu-a-7"
    dir_eqtl <<- "data/eqtls"
    dir_eqtlmaf <<- "data/eqtl_MAF"
    eqtl_outdir <<- "top_eqtls/"
    
    source('eQTL_query_parquet.r')
    
    list(src = outfile, alt = "this is alt text")
  

    
    
      
    
  })
  
  #Create Download Handler to provide results
  output$downloadData <- downloadHandler(
    filename = "resultsFile.txt",
    content = function(file) {
      write.table(read.table("top_eqtls/eQTLs_colocalized_w_GWAS.txt", header = TRUE, sep = "\t"), file) #Might need to make this code more reactive later to update with eQTL file as it changes.
      #write.table(top_eqtl_table, file) #Might need to make this code more reactive later
    }
  )
  
  
  #Attempt to show every image.
  #Complicated regex pattern that also breaks R compilation: chr\d+_gwas\d+_pos\d+_H4_\d+\.?\d+\.png
  # df_img <- data.frame(img_path = list.files(pattern = "chr.+png", full.names = TRUE))
  # n <- nrow(df_img)
  # 
  # observe({
  #   for (i in 1:n)
  #   {
  #     print(i)
  #     local({
  #       my_i <- i
  #       imagename = paste0("img", my_i)
  #       print(imagename)
  #       output[[imagename]] <-
  #         renderImage({
  #           list(src = file.path(df_img$img_path[my_i]),
  #                alt = "Image failed to render")
  #         }, deleteFile = FALSE)
  #     })
  #   }
  # })
  # 
  # 
  # output$peaks <- renderUI({
  # 
  #   #Might need to move algorithm code into this render function.
  #   
  #   
  #   image_output_list <-
  #     lapply(1:n,
  #            function(i)
  #            {
  #              imagename = paste0("img", i)
  #              imageOutput(imagename)
  #            })
  # 
  #   do.call(tagList, image_output_list)
  # })
  
  
}
  



#Can also be run with runApp("colocApp") when colocApp is a directory in the 
#working directory with this app.R file within it.
shinyApp(ui = ui, server = server)

#To exit out of the app, press escape or click the stop sign in the upper right.