library(shiny)
library(shinythemes)
#library(reticulate)

#Kenji/Colin imports below
# packages for parquet querying
library(arrow)
library(duckdb)
library(fs)
library(tidyverse) #Error here
library(DBI)
library(glue)
library(tictoc)


# packages for coloc
suppressPackageStartupMessages(suppressWarnings({
  library(gwasglue)
  library(gassocplot)
  library(dplyr)
  library(coloc)
  library(ggplot2)
  library(httpgd)
}))



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
                           #box(uiOutput("peaks"))
                           
                           #Offer method to download a readout of the results.
                           # Button
                           downloadButton("downloadData", "Download")
                  ),
                      
                  tabPanel("About", p("This app reads in "),
                           p("The flow chart below describes the backend algorithm.")
                           #INCLUDE FLOW CHART HERE
                           #img(src = "")
                  ),
                  tabPanel("Datasets", p(align = "center", "The primary eQTL database is __, and that can be found along with more supplementary info on our Github page.")),
                  tabPanel("Contact Us", p("This R Shiny Web App was developed by a team of Georgia Institute of Technology students: Andy Chea, Colin Naughton, Kenji Nishiura, and Jasmyn Pellebon. You can reach us at...")
                           )
      ),
      
      
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
    
    
    # study ID (needs to have value assigned by user)
    gwas_dataset <- gwasInput()
    dummy_dataset <- "ieu-a-7"
    # query API with study ID
    gwasinfo(id = as.character(gwas_dataset))
    
    dir_eqtl <- "data/eqtls"
    dir_eqtlmaf <- "data/eqtl_MAF"
    eqtl_outdir <- "top_eqtls/"
    
    # open parquet files
    ds_eQTL <- arrow::open_dataset(dir_eqtl, partitioning = "SNPChr")
    ds_eqtlMAF <- arrow::open_dataset(dir_eqtlmaf, partitioning = "hg19_chr")
    
    # clean up if previous connection if still open (only happens if script is interrupted)
    if (exists("con")) {
      duckdb_unregister(con, "eqtlTable")
      duckdb_unregister(con, "mafTable")
      dbDisconnect(con)
    }
    
    # connect duckdb to open parquet files
    con <- dbConnect(duckdb::duckdb())
    
    # register the datasets as DuckDB table, and give it a name
    duckdb::duckdb_register_arrow(con, "eqtlTable", ds_eQTL)
    duckdb::duckdb_register_arrow(con, "mafTable", ds_eqtlMAF)
    
    
    # define window size for coloc
    window_size <- as.integer(100000)
    
    # write header line to output file
    # write to a log file, timestamp, chr, gwas_pos, number of eQTLs in window, PP_H4
    write(paste0("Time", "\t", "chr", "\t", "pos_gwas", "\t", "num_eqtls_in_window", "\t", "PP_H4"), file = glue("coloc_log_window_{window_size}.txt"), append = FALSE)
    
    debug_mode=FALSE
    
    # to keep memory usage in check, work on one chromosome at a time
    for (chromosome in 1:22) {
      if (debug_mode) {
        if (chromosome < 6) {
          next
        }
      }
      
      # get the top GWAS hits for the current chromosome
      top <- ieugwasr::tophits(gwasInput()) %>%
        filter(chr == chromosome) %>%
        arrange(p)
      print(paste0(dim(top)[1], " GWAS top hits identified on chromosome ", chromosome))
      
      
      # iterate over gwas top hits for current chromosome
      for (tophit in seq_len(nrow(top))) {
        
        # if (debug_mode) {
        #     if (tophit < 15) {
        #         next
        #     }
        # }
        
        # start timer
        tic(glue("time elapsed for the {tophit}-th following GWAS top hit:"))
        
        print(top[tophit, ])
        # print the rsid field
        #print(top[tophit, "rsid"])
        
        pos_gwas <- top[tophit, ]$position
        lower <- pos_gwas - window_size
        
        # set floor of 0 for lower bound of window
        if (lower < 0) {
          lower <- 0
        }
        upper <- pos_gwas + window_size
        chrpos <- paste0(chromosome, ":", lower, "-", upper)
        
        # print the message below if debug_mode is TRUE
        if (debug_mode) {
          print("starting ieugwasr_to_coloc:")
        }
        
        out <- ieugwasr_to_coloc(
          id1 = as.character(gwas_dataset),
          id2 = as.character(dummy_dataset),
          chrompos = as.character(chrpos),
          type1 = "quant",
          # type2 = "cc" # dummy dataset
        )
        # drop dummy dataset2 from out
        out <- out[1]
        
        # define the SQL queries
        subquery_eqtl <- glue::glue(
          # window function query for top eQTLs, defined as:
          #   for each locus, determined by values of chromosome & position,
          #   sort descending by absolute value Zcore and take only top row)
          "
            SELECT * FROM(
                SELECT *, ROW_NUMBER() OVER (PARTITION BY SNPPos ORDER BY abs(Zscore) DESC) AS 'row'
                FROM eqtlTable
                WHERE SNPChr = {chromosome} AND SNPPos BETWEEN {lower} AND {upper}
                )
            WHERE row = 1
            "
        )
        
        subquery_maf <- glue::glue(
          "
            SELECT * FROM mafTable
            WHERE hg19_chr = {chromosome}
            "
        )
        
        # print the message below if debug_mode is TRUE
        if (debug_mode) {
          print("executing SQL queries:")
        }
        
        # execute the SQL queries
        result_eqtl <- dbGetQuery(con, subquery_eqtl)
        result_maf <- dbGetQuery(con, subquery_maf)
        
        # left join the eQTL and MAF tables
        result <- left_join(result_eqtl, result_maf,
                            by = c("SNPChr" = "hg19_chr", "SNPPos" = "hg19_pos"),
                            suffix = c(".eqtl", ".maf")
        ) %>%
          # keep columns Pvalue, SNP, SNPPos, Zscore, NrSamples, AlleleB_all
          dplyr::select(Pvalue, NrSamples, AlleleB_all, SNP.eqtl, Zscore, SNPPos, ) %>%
          # add back column SNPChr
          tibble::add_column(chr = as.character(chromosome), .before = "SNPPos") %>%
          # rename AlleleB_all to MAF
          dplyr::rename(MAF = AlleleB_all, snp = SNP.eqtl, pvalues = Pvalue, N = NrSamples, z = Zscore, pos = SNPPos) %>%
          # estimate beta/varbeta (citation?):
          #   beta = Zscore / sqrt(2 * Pvalue * (1 - Pvalue) * (NrSamples + Zscore**2))
          #   betaSE:  betaSE = 1 / sqrt(2 * Pvalue * (1 - Pvalue) * (NrSamples + Zscore**2))
          #   varbeta:  betaSE**2
          # create column of nulls for beta/varbeta
          dplyr::mutate(pos = as.integer(pos))
        
        # strip name from out for compatibility with join function
        names(out) <- NULL
        out <- as.data.frame(out)
        
        # ensure that only loci present in both GWAS and eQTL datasets are considered
        joined <- inner_join(out, result, by = "pos", suffix = c(".gwas", ".eqtl"))
        
        # separate joined dataframe back out to gwas and eqtl datasets
        result <- dplyr::select(joined, ends_with("eqtl"), pos) %>%
          dplyr::rename(pvalues = pvalues.eqtl, N = N.eqtl, MAF = MAF.eqtl, snp = snp.eqtl, z = z.eqtl, chr = chr.eqtl)
        out <- dplyr::select(joined, ends_with("gwas"), beta, varbeta, pos, id) %>%
          dplyr::rename(pvalues = pvalues.gwas, N = N.gwas, MAF = MAF.gwas, snp = snp.gwas, z = z.gwas, chr = chr.gwas)
        
        # convert dataframe to nested list of lists
        result <- as.list(result)
        out <- as.list(out)
        
        # note:  "type" and "id" fields are NOT lists, just strings
        result$type <- "quant"
        out$type <- "quant"
        result$id <- "eQTLGen_cis-eQTL"
        out$id <- out$id[1]
        
        # extract list of rsids to feed to coloc_to_gassocplot()
        rsid_list <- out[["snp"]]
        result <- result[c("pvalues", "N", "MAF", "type", "snp", "z", "chr", "pos", "id")]
        
        
        # if there are any null values in result print the NULL values
        # if (any(is.null(result$beta))) {
        #     print("NULL values in result:")
        #     print(result[is.null(result$beta)])
        # }
        # if (any(is.null(out$beta))) {
        #     print("NULL values in out:")
        #     print(out[is.null(out$beta)])
        # }
        
        result <- list(result)
        out <- list(out)
        
        # add name to outer list to match ieugwasr_to_coloc format
        names(result) <- "dataset2"
        
        # assemble gwas and eqtl datasets into one 'out' object in expected format
        out <- list(out[[1]], result[[1]])
        names(out) <- c("dataset1", "dataset2")
        
        # print the message below if debug_mode is TRUE
        if (debug_mode) {
          print("running Coloc:")
        }
        
        # run coloc
        res <- coloc::coloc.abf(out[[1]], result[[1]])
        
        PP_H4 <- res$summary[["PP.H4.abf"]]
        
        # write to a log file, timestamp, chr, gwas_pos, number of eQTLs in window, PP_H4
        write(paste0(Sys.time(), "\t", chromosome, "\t", pos_gwas, "\t", length(out[[1]]$pos), "\t", PP_H4), file = glue("coloc_log_window_{window_size}.txt"), append = TRUE)
        
        if (PP_H4 < 0.50) {
          print(glue("PP_H4 < 0.50 for eQTLs near GWAS locus chr{chromosome}:pos{pos_gwas}, skipping and continuing to next GWAS top hit"))
          next
        }
        H4 <- round(as.numeric(PP_H4), digits = 2)
        
        # get the minimum Pvalue from the eQTL table
        min_pval <- min(result_eqtl$Pvalue)
        # get the rsids for all the rows tied for the minimum Pvalue
        top_eqtl_table <- result_eqtl %>%
          dplyr::filter(Pvalue == min_pval) 
        
        str(top_eqtl_table)
        
        # if it doesn't exist, create directory
        if (!dir.exists(eqtl_outdir)) {
          dir.create(eqtl_outdir, recursive = TRUE)
        }
        # if the file doesn't exist already, write rows with headers from eQTL table to file when p-value is tied for minimum. else, append to file
        if (!file.exists(glue(eqtl_outdir,"eQTLs_colocalized_w_GWAS.txt"))) {
          write.table(top_eqtl_table, glue(eqtl_outdir,"eQTLs_colocalized_w_GWAS.txt"), sep = "\t", row.names = FALSE, quote = FALSE,col.names = TRUE, append = FALSE)
        } else {
          write.table(top_eqtl_table, glue(eqtl_outdir,"eQTLs_colocalized_w_GWAS.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        }    
        
        # print the message below if debug_mode is TRUE
        if (debug_mode) {
          print("running coloc_to_gassocplot:")
        }
        
        # API rejects requests if >500 rsids, so need to run plink locally
        if (length(out[[1]]$pos) >= 500) {
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
        # construct plot
        theplot <- gassocplot::stack_assoc_plot(temp$markers, temp$z, temp$corr, traits = temp$traits)
        
        base_dir <- glue("plots/{as.integer(window_size)}")
        # output path for saving figure
        outfile <- glue(base_dir, "/chr{chromosome}_gwas{sprintf('%03d', tophit)}_pos{pos_gwas}_H4_{H4}.png")
        # if it doesn't exist, create directory
        if (!dir.exists(base_dir)) {
          dir.create(base_dir, recursive = TRUE)
        }
        # save plot w/ base R functions
        png(filename = outfile)
        plot(theplot)
        dev.off()
        
        # function to save plot to file (not working:  output is blank square)
        # stack_assoc_plot_save(theplot, outfile, 2, width = 3, height = 3, dpi = 500)
        
        # stop timer
        # toc(log = TRUE)
      } # close loop over each gwas top hit
    } # close loop over each chromosome
    
    # clean up
    duckdb_unregister(con, "eqtlTable")
    duckdb_unregister(con, "mafTable")
    dbDisconnect(con, shutdown = TRUE)
    
    #COLIN/KENJI CODE ABOVE
    
    
    
    
    list(src = outfile, alt = "this is alt text")
      
    
  })
  
  
  #Attempt to show every image.
  # df_img <- data.frame(id = c(1:5), img_path = c("h1000.png", "h2000.png", "h3000.png", "h4000.png", "h000.png"))
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
  #           list(src = file.path('www', df_img$img_path[my_i]), 
  #                width = "100%", height = "55%",
  #                alt = "Image failed to render")
  #         }, deleteFile = FALSE)
  #     })
  #   }
  # })
  # 
  # 
  # output$peaks <- renderUI({
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
  
  
  
  #Create Download Handler to provide results
  output$downloadData <- downloadHandler(
    filename = resultsFile.txt,
    content = function(file) {
      write.table(read.table("top_eqtls/eQTLs_colocalized_w_GWAS.txt", header = TRUE, sep = "\t"), file) #Might need to make this code more reactive later to update with eQTL file as it changes.
    }
  )
  # 
  
}
  



#Can also be run with runApp("colocApp") when colocApp is a directory in the 
#working directory with this app.R file within it.
shinyApp(ui = ui, server = server)

#To exit out of the app, press escape or click the stop sign in the upper right.