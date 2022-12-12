library(shiny)
library(shinythemes)
library(shinycssloaders)
library(DT)
library(slickR)

# in a sidebar, take in input from user: 
# dropdown menu for eQTL database
# text input for GWAS ID with submit button

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
        # dropdown of numeric integers between 1-22 for chromosomes
        selectInput("chr_input", "Pick a chromosome",
            c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22),
            selected = TRUE,
            multiple = TRUE,
            selectize = TRUE,
            width = NULL,
            size = NULL),
        actionButton("submit", "Submit"),
        ),
        mainPanel(
            width=9,
            # hide gwas table if submiit button is pressed
            conditionalPanel(
                condition = "input.submit == 0",
                withSpinner(DT::dataTableOutput("gwas_table"))
            ),
            slickROutput("slickr", width="500px"),
            # add vertical whitespace
            tags$br(),
            tags$br(),
            tags$br(),
            tags$br(),
            tags$br(),
            # placeholder for error text when 0 colocalized eqtls are found
            textOutput(outputId = "error"),
            # format error text
            tags$head(tags$style("#error{color: red;
                                font-size: 20px;
                                font-style: italic;}"
                )
            ),
            # create a download button for the plots and top_eqtls in user_output
            conditionalPanel(
                condition = "input.submit >= 1",
                # download buttons for plot and top eqtls
                downloadButton("dl_plots", "Download Plots"),
                downloadButton("dl_colocalized", "Download Colocalized GWAS/eQTLs")
                )
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

    # # print the working directory
    # imgs <- list.files(glue::glue(plot_dir,'/',user_outdir), pattern=".png", full.names = TRUE)
    # # str(imgs)
    # output[["slickr"]] <- renderSlickR({
    #     slickR(imgs) + settings(slidesToShow = 1, slidesToScroll = 1, lazyLoad = 'anticipated', dots = TRUE, arrows = TRUE, infinite = FALSE)
    # })
    

    # when the submit button is clicked, execute pipeline and refresh the image viewer between each generated plot
    observeEvent(input$submit, {
        # clear error message 
        output$error <- renderText({""})
        gwas_dataset <<- input$gwas_input
        chromosomes <<- input$chr_input
        dir_ld <<- "ld"
        dir_eqtl <<- "data/eqtls"
        dir_eqtlmaf <<- "data/eqtl_MAF"
        eqtl_outdir <<- "top_eqtls"
        plot_dir <<- "plots"
        # if it doesn't exist, create a directory named user_outdir in plots and top_eqtls
        # if it exists already, delete files from previous run
        # if user_outdir is defined and the plot_dir/user_outdir exists, delete the files in that directory
        if (exists("user_outdir") && dir.exists(glue::glue(plot_dir,'/', user_outdir))) {
            unlink(glue::glue(plot_dir,'/', user_outdir, '/*'))
            unlink(glue::glue(eqtl_outdir,'/', user_outdir, '/*'))
        }

        # generate unique id using timestamp
        user_outdir <<- as.character(Sys.time())
        user_outdir <<- gsub(" ", "_", user_outdir)
        user_outdir <<- gsub(":", "_", user_outdir)
        user_outdir <<- gsub("-", "_", user_outdir)

        # if user subdirectories don't exist yet, create them
        if (!dir.exists(glue::glue(plot_dir,'/', user_outdir))) {
            dir.create(glue::glue(plot_dir,'/', user_outdir))
        }
        if (!dir.exists(glue::glue(eqtl_outdir,'/', user_outdir))) {
            dir.create(glue::glue(eqtl_outdir,'/', user_outdir))
        }

        
        window_size <<- as.integer(100000)
        source('ieugwasr_to_coloc_modified.r')
        source("gassocplot_modified.r")
        source('eQTL_query_parquet.r')
        con <- initialize_db()
        # query API for top GWAS hits for all chromosomes
        top_all <- ieugwasr::tophits(gwas_dataset)

        # initialize counter at 0 each time submit button is pressed
        colocalized_counter <- 0

        for (chromosome in chromosomes) {    
            # filter the top GWAS hits for the current chromosome
            top <- top_all %>%
                filter(chr == chromosome) %>%
                arrange(p)
            # initialize progressbar for each chromosome separately
            withProgress(
                message='Please wait',
                detail=glue::glue('Processing chr {chromosome}:'),
                value=0, {
                for (tophit_idx in seq_len(nrow(top))) {
                    # calls customized ieugwasr_to_coloc function to take in eQTL data as direct input
                    out_PPH4_rawResult <- gather_and_format_gwas_eqtl_in_region(con, chromosome, top[tophit_idx,])
                    # unpack returned objects
                    out <- out_PPH4_rawResult[[1]]
                    PP_H4 <- out_PPH4_rawResult[[2]]
                    # continue to the next window if posterior probability of H4 is less than 0.5
                    if (PP_H4 < 0.50) {
                        print(glue::glue("PP_H4 < 0.50, skipping and continuing to next GWAS top hit"))
                        setProgress(tophit_idx / nrow(top), detail = glue::glue('Processing chr {chromosome}: region {tophit_idx} out of {nrow(top)}'))
                        next
                    }
                    # API rejects requests if >500 rsids, so need to run plink locally
                    if (length(out[[2]]$pos) >= 500) {
                        print("too many rsids, running plink locally")
                        # input to coloc_to_gassocplot is list of rsids (should be identical in gwas/eqtl data at this stage): out[[1]]$snp
                        # choices for ancestry are AMR, AFR, EAS, EUR, SAS
                        # note: a bit slow first time because the plink_bin function will download/install plink if it's not already installed.
                        temp <- coloc_to_gassocplot(out, bfile = paste0(dir_ld, "/EUR"), plink_bin = '/projects/team1/bin/mambaforge/envs/plink/bin/plink')
                        # local debugging (comment out line above, uncomment line below)
                        # temp <- coloc_to_gassocplot(out, bfile = paste0(dir_ld, "/EUR"), plink_bin = genetics.binaRies::get_plink_binary())
                        # increment counter for number of colocalized eQTLs
                        colocalized_counter <- colocalized_counter + 1
                    } else {
                        print("running coloc_to_gassocplot via API")
                        # query the API if <500 rsids
                        temp <- coloc_to_gassocplot(out)
                        # increment counter for number of colocalized eQTLs
                        colocalized_counter <- colocalized_counter + 1
                    }
                    # get the rsids matching the gassocplot plot labels
                    top_marker_gwas <- extract_top_markers_from_gassocplot_output(temp)[1]
                    top_marker_eqtl <- extract_top_markers_from_gassocplot_output(temp)[2]
                    
                    result_raw <- out_PPH4_rawResult[[3]]

                    top_table <<- get_top_marker_raw_data(top_marker_gwas, top_marker_eqtl, result_raw, PP_H4)
                    top_table_filename <<- "eQTLs_colocalized_w_GWAS.txt"
                    top_table_path <<- glue::glue(eqtl_outdir, '/', user_outdir, '/', top_table_filename)
                    # print(top_table_path)
                    # write the header only once
                    if (!file.exists(top_table_path)) {
                        write.table(top_table, top_table_path, sep = "\t", row.names = FALSE, quote = FALSE,col.names = TRUE, append = FALSE)
                    } else {
                        write.table(top_table, top_table_path, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
                    }    
                    # call customized gassocplot function w/ -log(pval) ceiling removed
                    plot_associated_signals_and_save(tophit_idx, chromosome, PP_H4, temp, glue::glue(plot_dir,'/', user_outdir))

                    
                    # increment progress bar
                    print(glue::glue('working on tophit_idx {tophit_idx} out of {nrow(top)}'))
                    setProgress(tophit_idx / nrow(top), detail = glue::glue('Processing chr {chromosome}: region {tophit_idx} out of {nrow(top)}'))
                    
                    } # close loop over each gwas top hit
                # redraw the main panel explicitly
                imgs <<- list.files(glue::glue(plot_dir,'/', user_outdir), pattern=".png", full.names = TRUE)
                # str(imgs)
                output[["slickr"]] <<- renderSlickR({
                    slickR(imgs) + settings(slidesToShow = 1, slidesToScroll = 1, lazyLoad = 'anticipated', dots = TRUE, arrows = TRUE, infinite = FALSE)
                })
            }) # close withProgress
        } # close loop over chromosomes
        # close the database connection
        clean_up(con)        
        # if colocalized_counter is 0, print message
        if (colocalized_counter == 0) {
            print("No colocalized eQTLs found. Try again with different chromosomes or GWAS datasets.")
            # set error text when 0 colocalized eqtls are found
            output$error <- renderText({
                "No colocalized eQTLs found. Try again with different chromosomes or GWAS datasets."
            })
        } else {
            print(glue::glue("Number of colocalized eQTLs: {colocalized_counter}"))
        }
        # if top_table is not defined, print message
        # if (!exists("top_table")) {
        #     print("top_table not defined")
        #     # set error text when 0 colocalized eqtls are found
        #     output$error <- renderText({
        #         "No colocalized eQTLs found. Try again with different chromosomes or GWAS datasets."
        #     })
        # } else {
        #     # print number of observations in top_table
        #     print(top_table.nrow)
        #     # str(top_table)
        #     # print(dim(top_table))
        # }
        # get basenames of the paths in imgs
        plot_names <- basename(imgs)
        # zip plots and store in user directory
        zip::zip(zipfile = 'plots.zip', files = plot_names, root = glue::glue(plot_dir,'/', user_outdir), recurse = FALSE)
        zip_path <<- glue::glue(plot_dir,'/', user_outdir, '/plots.zip')
        
    }) # close submit button trigger
    # when dl_colocalized or dl_plots is clicked, send existing file to user
    output$dl_colocalized <- downloadHandler(
        filename = function(){
            paste("colocalized","txt",sep=".")
        },
        content = function(con){
            file.copy(top_table_path, con)
    })
    output$dl_plots <- downloadHandler(
        filename = function(){
            paste("plots","zip",sep=".")
        },
        content = function(con){
            file.copy(zip_path, con)
    })


}

shinyApp(ui, server)
