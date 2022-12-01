library(shiny)
library(shinythemes)
library(shinycssloaders)
library(DT)

# in a sidebar, take in input from user: 
# dropdown menu for eQTL database
# text input for GWAS ID with actionbutton trigger

ui <- fluidPage(
        # theme = shinytheme("superhero"),
                titlePanel("Colocalization App"),
                sidebarLayout(
                  sidebarPanel(width=3,
                    selectInput("eqtl_input", "Pick a eQTL mapping database",
                                c("eQTLGen"),
                                selected = NULL,
                                multiple = FALSE,
                                selectize = TRUE,
                                width = NULL,
                                size = NULL),
                    textInput("gwas_input", "Enter a GWAS ID:", width = NULL, placeholder = "ex: ieu-b-30"),
                    # tags$a(href="https://gwas.mrcieu.ac.uk/datasets/?gwas_id__icontains=ukb-b&sort=-gwas_id&page=73", "Find a GWAS ID at this link!"),
                    # dropdown of numeric integers between 1-22 for chromosomes
                    selectInput("chr_input", "Pick a chromosome",
                                c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22),
                                selected = TRUE,
                                multiple = TRUE,
                                selectize = TRUE,
                                width = NULL,
                                size = NULL),
                         
                    actionButton("submit", "Submit")
                    ),
                    mainPanel(width=9,
                        # render text of gwasID and choice
                        textOutput("eqtl_choice"),
                        textOutput("gwas_choice"),
                        textOutput("chromosome_choice"),
                        textOutput("gwas_dataset"),
                        # render datatable with the output of ieugwasr::gwasinfo()
                        withSpinner(DT::dataTableOutput("gwas_table")

                    )
                )
))

server <- function(input, output, session) {
    # use ieugwasr::gwasinfo() to get the gwas info
    gwas_info <- ieugwasr::gwasinfo()
    # return the gwas info as a datatable and only keep certain columns
    gwas_info <- gwas_info[, c("id", "category","subcategory", "trait", "population", "ontology", "consortium", "author", "year")]
    # print all the unique values of 'population'
    # print(unique(gwas_info$population))
    output$gwas_table <- DT::renderDataTable(
        gwas_info, filter="top", rownames=FALSE, height=800, selection = 'single',options = list(dom='Blfrtip',pageLength = 50, scrollX = TRUE, scrollY = 300, scrollCollapse = TRUE, autoWidth = TRUE, columnDefs = list(list(className = 'dt-center', targets = '_all'))))
    # when a row is selected, update the textinput labelled gwas_input with the id column value of the selected row
    observeEvent(input$gwas_table_rows_selected, {
        gwas_dataset <<- gwas_info[input$gwas_table_rows_selected, "id"]$id
        # ancestry <<- gwas_info[input$gwas_table_rows_selected, "population"]$population
        # print(gwas_dataset)
        # print(ancestry)
        updateTextInput(session, "gwas_input", value = gwas_dataset)
    })
    # when the submit button is clicked, 
    # set gwas_dataset to contents of textinput labelled gwas_input, 
    # and save vector of chromosomes to variable
    observeEvent(input$submit, {
        gwas_dataset <<- input$gwas_input
        chromosomes <<- input$chr_input
        print(gwas_dataset)
        print(chromosomes)
        dir_ld <<- "ld"
        dir_eqtl <<- "data/eqtls"
        dir_eqtlmaf <<- "data/eqtl_MAF"
        eqtl_outdir <<- "top_eqtls/"
        source('eQTL_query_parquet.r')
    })
}

shinyApp(ui, server)







