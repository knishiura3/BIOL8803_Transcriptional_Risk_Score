library(shiny)
library(shinythemes)
library(shinycssloaders)
library(DT)
library(slickR)

# in a sidebar, take in input from user: 
# dropdown menu for eQTL database
# text input for GWAS ID with actionbutton trigger

ui <- fluidPage(
    # theme = shinytheme("superhero"),
    titlePanel("Colocalization App"),
    # custom styling to position progressbar
    tags$head(tags$style(".shiny-notification {position: fixed; top: 50%; left: 50%; font-size: 20px}")),
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
        actionButton("submit", "Submit"),
        # create actionbutton to refresh the image list
        # actionButton("refresh", "Refresh")
        ),
        mainPanel(
            width=9,
            # hide gwas table if submiit button is pressed
            conditionalPanel(
                condition = "input.submit == 0",
                withSpinner(DT::dataTableOutput("gwas_table"))
            ),
            # tags$div(
                slickROutput("slickr", width="500px"),
                # style = "margin-left:auto;margin-right:auto;"
                # ) 

            # conditionalPanel(
            #     condition = "input.submit == 1",
            #     tags$div(
            #     slickROutput("slickr", width="500px"),
            #     style = "margin-left:100px;"
            #     )   
            # ),
        )
    )
)

server <- function(input, output, session) {
    # use ieugwasr::gwasinfo() to query API for the possible gwas datasets
    gwasinfo_df <- ieugwasr::gwasinfo()
    # str(gwasinfo_df)
    # return the gwas info as a datatable and only keep certain columns
    gwasinfo_render_df <- gwasinfo_df[, c("id", "category","subcategory", "trait", "population", "ontology", "consortium", "author", "year")]
    # print all the unique values of 'population'
    # print(unique(gwas_info$population))
    output$gwas_table <- DT::renderDataTable(
        gwasinfo_render_df, filter="top", rownames=FALSE, height=800, selection = 'single',options = list(dom='Blfrtip',pageLength = 50, scrollX = TRUE, scrollY = 300, scrollCollapse = TRUE, autoWidth = TRUE, columnDefs = list(list(className = 'dt-center', targets = '_all'))))
    # when a row is selected, update the textinput labelled gwas_input with the id column value of the selected row
    observeEvent(input$gwas_table_rows_selected, {
        gwas_dataset <<- gwasinfo_render_df[input$gwas_table_rows_selected, "id"]$id
        # ancestry <<- gwas_info[input$gwas_table_rows_selected, "population"]$population
        # print(gwas_dataset)
        # print(ancestry)
        updateTextInput(session, "gwas_input", value = gwas_dataset)
    })
    # when the refresh button is clicked,
    # list the images in the plots directory again and update the image list
    # observeEvent(input$refresh, {
    #     # list the images in the plots directory
    #     imgs <- list.files("plots/100000", pattern=".png", full.names = TRUE)
    #     str(imgs)
    #     output[["slickr"]] <<- renderSlickR({
    #     slickR(imgs) + settings(slidesToShow = 1, slidesToScroll = 1, lazyLoad = 'anticipated', dots = TRUE, arrows = TRUE, infinite = FALSE)
    #     })
    # })
    plot_dir = "plots"
    # print the working directory
    imgs <- list.files(plot_dir, pattern=".png", full.names = TRUE)
    str(imgs)
    output[["slickr"]] <- renderSlickR({
        slickR(imgs) + settings(slidesToShow = 1, slidesToScroll = 1, lazyLoad = 'anticipated', dots = TRUE, arrows = TRUE, infinite = FALSE)
    })
    
    # when the submit button is clicked, execute pipeline and refresh the image viewer between each generated plot
    observeEvent(input$submit, {
        gwas_dataset <<- input$gwas_input
        chromosomes <<- input$chr_input
        dir_ld <<- "ld"
        dir_eqtl <<- "data/eqtls"
        dir_eqtlmaf <<- "data/eqtl_MAF"
        eqtl_outdir <<- "top_eqtls"
        
        window_size <<- as.integer(100000)
        source('ieugwasr_to_coloc_modified.r')
        source("gassocplot_modified.r")
        source('eQTL_query_parquet.r')
        con <- initialize_db()
        # query API for top GWAS hits for all chromosomes
        top_all <- ieugwasr::tophits(gwas_dataset)
        for (chromosome in chromosomes) {
            # filter the top GWAS hits for the current chromosome
            top <- top_all %>%
                filter(chr == chromosome) %>%
                arrange(p)
                # wrap the loop execution in withProgress

            withProgress(
                message='Please wait',
                detail=glue('Processing chr {chromosome}:'),
                value=0, {
                for (tophit_idx in seq_len(nrow(top))) {
                    out_PPH4_rawResult <- gather_and_format_gwas_eqtl_in_region(con, chromosome, top[tophit_idx,])
                    # assign first element of list to out
                    out <- out_PPH4_rawResult[[1]]
                    # str(out)
                    # assign second element of list to PP_H4
                    PP_H4 <- out_PPH4_rawResult[[2]]
                    # continue to the next window if posterior probability of H4 is less than 0.5
                    if (PP_H4 < 0.50) {
                        print(glue("PP_H4 < 0.50, skipping and continuing to next GWAS top hit"))
                        setProgress(tophit_idx / nrow(top), detail = glue('Processing chr {chromosome}: region {tophit_idx} out of {nrow(top)}'))
                        next
                    }

                    # print the message below if debug_mode is TRUE
                    if (debug_mode) {
                        print("running coloc_to_gassocplot:")
                    }
                    
                    # API rejects requests if >500 rsids, so need to run plink locally
                    if (length(out[[2]]$pos) >= 500) {
                        print("too many rsids, running plink locally")
                        # input to coloc_to_gassocplot is list of rsids (should be identical in gwas/eqtl data at this stage): out[[1]]$snp
                        # choices for ancestry are AMR, AFR, EAS, EUR, SAS
                        # note: a bit slow first time because the plink_bin function will download/install plink if it's not already installed.
                        temp <- coloc_to_gassocplot(out, bfile = paste0(dir_ld, "/EUR"), plink_bin = genetics.binaRies::get_plink_binary())
                    } else {
                        print("running coloc_to_gassocplot via API")
                        # query the API if <500 rsids
                        temp <- coloc_to_gassocplot(out)
                    }
                    # get the rsids matching the gassocplot plot labels
                    top_marker_gwas <- extract_top_markers_from_gassocplot_output(temp)[1]
                    top_marker_eqtl <- extract_top_markers_from_gassocplot_output(temp)[2]
                    
                    result_raw <- out_PPH4_rawResult[[3]]

                    top_table <- get_top_marker_raw_data(top_marker_gwas, top_marker_eqtl, result_raw, PP_H4)
                    top_table_filename <- "eQTLs_colocalized_w_GWAS.txt"
                    # write the header only once
                    if (!file.exists(glue(eqtl_outdir, top_table_filename))) {
                        write.table(top_table, glue(eqtl_outdir,top_table_filename), sep = "\t", row.names = FALSE, quote = FALSE,col.names = TRUE, append = FALSE)
                    } else {
                        write.table(top_table, glue(eqtl_outdir,top_table_filename), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
                    }    

                    plot_associated_signals_and_save(tophit_idx, chromosome, PP_H4, temp, plot_dir)

                    
                    # increment progress bar
                    print(glue('working on tophit_idx {tophit_idx} out of {nrow(top)}'))
                    setProgress(tophit_idx / nrow(top), detail = glue('Processing chr {chromosome}: region {tophit_idx} out of {nrow(top)}'))
                    
                    } # close loop over each gwas top hit
                # redraw the main panel explicitly
                imgs <<- list.files(plot_dir, pattern=".png", full.names = TRUE)
                str(imgs)
                output[["slickr"]] <<- renderSlickR({
                    slickR(imgs) + settings(slidesToShow = 1, slidesToScroll = 1, lazyLoad = 'anticipated', dots = TRUE, arrows = TRUE, infinite = FALSE)
                })
            }) # close withProgress
        } # close loop over chromosomes
        
        clean_up(con)        
    })
}

shinyApp(ui, server)







