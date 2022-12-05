# packages for parquet querying
library(arrow)
library(duckdb)
library(fs)
library(tidyverse)
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

source('ieugwasr_to_coloc_modified.r')
source("gassocplot_modified.r")

initialize_db <- function() {
    # open parquet files
    ds_eQTL <- arrow::open_dataset(dir_eqtl, partitioning = "SNPChr")
    ds_eqtlMAF <- arrow::open_dataset(dir_eqtlmaf, partitioning = "hg19_chr")

    # clean up if previous connection if still open (only happens if script is interrupted)
    if (exists("con")) {
        duckdb::duckdb_unregister(con, "eqtlTable")
        duckdb::duckdb_unregister(con, "mafTable")
        duckdb::dbDisconnect(con)
    }

    # connect duckdb to open parquet files
    con <- duckdb::dbConnect(duckdb::duckdb())

    # register the datasets as DuckDB table, and give it a name
    duckdb::duckdb_register_arrow(con, "eqtlTable", ds_eQTL)
    duckdb::duckdb_register_arrow(con, "mafTable", ds_eqtlMAF)
    return(con)
}

gather_and_format_gwas_eqtl_in_region <- function(gwas_hit) {
    # str(gwas_hit)
    # get the position of the current GWAS hit
    gwas_pos <- gwas_hit$position
    lower <- gwas_pos - window_size
    # set floor of 0 for lower bound of window
    if (lower < 0) {
        lower <- 0
    }
    upper <- gwas_pos + window_size
    chrpos <- paste0(chromosome, ":", lower, "-", upper)

    # window function query for top eQTLs, defined as:
    #   for each locus, determined by values of chromosome & position,
    #   sort descending by absolute value Zcore and take only top row)
    subquery_eqtl <- glue::glue(
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
    # execute the SQL queries
    result_eqtl <- dbGetQuery(con, subquery_eqtl)
    result_maf <- dbGetQuery(con, subquery_maf)
    print('SQL queries executed')
    # left join the eQTL and MAF tables to get full set of relevant eqtl data
    result_raw <- left_join(result_eqtl, result_maf,
        by = c("SNPChr" = "hg19_chr", "SNPPos" = "hg19_pos"),
        suffix = c(".eqtl", ".maf")
    ) 
    print('join complete')
    # filter and rename columns
    result <- result_raw %>%
        # keep columns Pvalue, SNP, SNPPos, Zscore, NrSamples, AlleleB_all
        dplyr::select(Pvalue, NrSamples, AlleleB_all, SNP.eqtl, Zscore, SNPPos, AssessedAllele,OtherAllele) %>%
        # add back column SNPChr
        tibble::add_column(chr = as.character(chromosome), .before = "SNPPos") %>%
        # rename AlleleB_all to MAF
        dplyr::rename(eaf = AlleleB_all, rsid = SNP.eqtl, p = Pvalue, n = NrSamples, z = Zscore, position = SNPPos, ea = AssessedAllele, nea = OtherAllele)
    
    # str(result)
    # harmonize gwas and eqtl datasets and structure for downstream analysis  
    out <- ieugwasr_to_coloc_custom(
        id1= gwas_dataset, 
        result, # eqtl data
        type2 = "quant", 
        chrompos = as.character(chrpos)
    )
    
    # str(out)
    # run coloc
    res <- coloc::coloc.abf(out[[1]], out[[2]])
    # str(res)
    # extract H4 from coloc output
    PP_H4 <- round(as.numeric(res$summary[["PP.H4.abf"]]), digits = 2)
    
    return(list(out, PP_H4, result_raw))
}

extract_top_markers_from_gassocplot_output <- function(temp) {
    # str(temp)
    # code modified from gassocplot for identifying top markers
    # https://github.com/jrs95/gassocplot/blob/master/R/figures.R#L387-L390
    mlog10p_gwas <- -(log(2) + pnorm(-abs(temp$z[[1]]), log.p=TRUE))/log(10)
    mlog10p_eqtl <- -(log(2) + pnorm(-abs(temp$z[[2]]), log.p=TRUE))/log(10)
    # add as column of temp object because of shared indices
    temp$stats_gwas <- mlog10p_gwas
    temp$stats_eqtl <- mlog10p_eqtl
    # get the values and indices of the the maximum value of temp$stats (or elements of values/indices tied for max)
    max_pval_gwas <- max(temp$stats_gwas)
    max_pval_eqtl <- max(temp$stats_eqtl)
    # str(max_pval_eqtl)
    max_pval_indices_gwas <- which(temp$stats_gwas == max_pval_gwas)
    max_pval_indices_eqtl <- which(temp$stats_eqtl == max_pval_eqtl)
    # str(max_pval_indices_eqtl)
    # get the rsid (temp$markers$marker) at the index position of top marker
    top_marker_gwas <- temp$markers$marker[max_pval_indices_gwas]
    top_marker_eqtl <- temp$markers$marker[max_pval_indices_eqtl]

    return(list(top_marker_gwas, top_marker_eqtl))
}

# take top markers as input and return a raw data row for each colocalized gwas and eqtl pair
get_top_marker_raw_data <- function(top_marker_gwas, top_marker_eqtl, result_raw, PP_H4) {
    # create a table of the row in result_raw that corresponds to the top marker and append H4 as a column
    # and only keep desired columns
    top_eqtl_table <- result_raw[result_raw$SNP.eqtl %in% top_marker_eqtl, c("GeneChr","Pvalue", "SNP.eqtl", "SNPPos", "AssessedAllele", "OtherAllele", "Zscore", "Gene", "GenePos", "NrCohorts", "NrSamples", "FDR", "BonferroniP", "AlleleA", "AlleleB", "allA_total", "allAB_total", "allB_total", "AlleleB_all")]
    # rename GeneChr to Chr
    colnames(top_eqtl_table)[colnames(top_eqtl_table) == "GeneChr"] <- "Chr"
    # loop over and add .eqtl suffix to columns Pvalue, SNPPos, AssessedAllele, OtherAllele, Zscore, AlleleA, AlleleB, allA_total, allAB_total, allB_total, AlleleB_all
    for (i in 1:ncol(top_eqtl_table)) {
        # skip column name that already has .eqtl suffix or if Chr
        if (grepl(".eqtl", colnames(top_eqtl_table)[i]) | colnames(top_eqtl_table)[i] == "Chr") {
            next
        }
        colnames(top_eqtl_table)[i] <- paste0(colnames(top_eqtl_table)[i], ".eqtl")
    }

    # queries GWAS table w/ rsid for top marker
    top_gwas <- ieugwasr::associations(top_marker_gwas,gwas_dataset)
    # only keep columns position, beta, se, p, n, id, rsid, ea, nea, eaf, trait
    top_gwas <- top_gwas[, c("position", "beta", "se", "p", "n", "id", "rsid", "ea", "nea", "eaf", "trait")]
    colnames(top_gwas) <- paste0(colnames(top_gwas), ".gwas")
    # combine top_gwas and top_eqtl_table side-by-side in new table
    top_table <- cbind(top_gwas, top_eqtl_table)
    # move Chr column to first position
    top_table <- top_table[, c("Chr", colnames(top_table)[!colnames(top_table) %in% "Chr"])]
    top_table$H4 <- PP_H4
    return(top_table)
}

plot_associated_signals_and_save <- function(temp) {
    # construct plot
    theplot <- stack_assoc_plot_custom(temp$markers, temp$z, temp$corr, traits = temp$traits)
    
    plot_dir <- glue("plots/{as.integer(window_size)}")
    # output path for saving a figure for each colocalized pair of GWAS/eQTL 
    outfile <- glue(plot_dir, "/chr{chromosome}_gwastophit{sprintf('%03d', tophit_idx)}_H4_{PP_H4}.png")
    # if it doesn't exist, create directory
    if (!dir.exists(plot_dir)) {
        dir.create(plot_dir, recursive = TRUE)
    }
    # save plot w/ base R functions
    png(filename = outfile)
    plot(theplot)
    dev.off()
}

# function to clean up
clean_up <- function() {
    # clean up
    print("cleaning up")
    duckdb_unregister(con, "eqtlTable")
    duckdb_unregister(con, "mafTable")
    dbDisconnect(con, shutdown = TRUE)
    print("done")
}




# Reference LD Panels
#   need to set this to wherever you want the LD panel stored
# dir_ld <- "/home/kenji/BIOL8803/BACKUP_BIOL8803_Transcriptional_Risk_Score/ld_reference"

# Copied from the documentation at https://mrcieu.github.io/gwasglue/index.html
#
# Updated 1000 genomes LD reference panels (multiple populations):
# http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz (1.5GB)
#   <contents>/
#   ├── AFR.bed
#   ├── AFR.bim
#   ├── AFR.fam
#   ├── AMR.bed
#   ├── AMR.bim
#   ├── AMR.fam
#   ├── EAS.bed
#   ├── EAS.bim
#   ├── EAS.fam
#   ├── EUR.bed
#   ├── EUR.bim
#   ├── EUR.fam
#   ├── SAS.bed
#   ├── SAS.bim
#   └── SAS.fam


debug_mode=FALSE


# if debug mode is TRUE, set these values
if (debug_mode) {
    # study ID (needs to have value assigned by user)
    gwas_dataset <- "ieu-b-30"
    # dummy_dataset <- "ieu-a-7"
    dir_ld <- "/home/kenji/BIOL8803/BACKUP_BIOL8803_Transcriptional_Risk_Score/ld_reference"
    # query API with study ID
    gwasinfo(id = as.character(gwas_dataset))
    chromosomes <- c(21)

    dir_eqtl <- "/home/kenji/BIOL8803/BACKUP2_IOL8803_Transcriptional_Risk_Score/eqtls_merged"
    dir_eqtlmaf <- "/home/kenji/BIOL8803/BACKUP_BIOL8803_Transcriptional_Risk_Score/eqtl_MAF"

    eqtl_outdir <- "top_eqtls/"
}

# # open parquet files
# ds_eQTL <- arrow::open_dataset(dir_eqtl, partitioning = "SNPChr")
# ds_eqtlMAF <- arrow::open_dataset(dir_eqtlmaf, partitioning = "hg19_chr")

# # clean up if previous connection if still open (only happens if script is interrupted)
# if (exists("con")) {
#     duckdb::duckdb_unregister(con, "eqtlTable")
#     duckdb::duckdb_unregister(con, "mafTable")
#     duckdb::dbDisconnect(con)
# }

# # connect duckdb to open parquet files
# con <- duckdb::dbConnect(duckdb::duckdb())

# # register the datasets as DuckDB table, and give it a name
# duckdb::duckdb_register_arrow(con, "eqtlTable", ds_eQTL)
# duckdb::duckdb_register_arrow(con, "mafTable", ds_eqtlMAF)
con <- initialize_db()

# define window size for coloc
window_size <- as.integer(100000)

# write header line to output file
# split chrpos on colon
# pos_gwas <- stringr::str_split(chrpos, ":")[2]
# write to a log file, timestamp, chr, gwas_pos, number of eQTLs in window, PP_H4
# write(paste0("Time", "\t", "chr", "\t", "pos_gwas", "\t", "num_eqtls_in_window", "\t", "PP_H4"), file = glue("coloc_log_window_{window_size}.txt"), append = FALSE)

top_all <- ieugwasr::tophits(gwas_dataset)
# to keep memory usage in check, work on one chromosome at a time
# for (chromosome in 1:22) {
for (chromosome in chromosomes) {

    # get the top GWAS hits for the current chromosome
    top <- top_all %>%
        filter(chr == chromosome) %>%
        arrange(p)
    print(paste0(dim(top)[1], " GWAS top hits identified on chromosome ", chromosome))

    # iterate over gwas top hits for current chromosome
    for (tophit_idx in seq_len(nrow(top))) {
        out_PPH4_rawResult <- gather_and_format_gwas_eqtl_in_region(top[tophit_idx,])
        # assign first element of list to out
        out <- out_PPH4_rawResult[[1]]
        # str(out)
        # assign second element of list to PP_H4
        PP_H4 <- out_PPH4_rawResult[[2]]
        # continue to the next window if posterior probability of H4 is less than 0.5
        if (PP_H4 < 0.50) {
            print(glue("PP_H4 < 0.50, skipping and continuing to next GWAS top hit"))
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
        # str(temp$z[[1]])
        # str(temp$z[[2]])
        # # code modified from gassocplot for identifying top markers
        # # https://github.com/jrs95/gassocplot/blob/master/R/figures.R#L387-L390
        # mlog10p_gwas <- -(log(2) + pnorm(-abs(temp$z[[1]]), log.p=TRUE))/log(10)
        # mlog10p_eqtl <- -(log(2) + pnorm(-abs(temp$z[[2]]), log.p=TRUE))/log(10)
        # # add as column of temp object because of shared indices
        # temp$stats_gwas <- mlog10p_gwas
        # temp$stats_eqtl <- mlog10p_eqtl
        # # get the values and indices of the the maximum value of temp$stats (or elements of values/indices tied for max)
        # max_pval_gwas <- max(temp$stats_gwas)
        # max_pval_eqtl <- max(temp$stats_eqtl)
        # # str(max_pval_eqtl)
        # max_pval_indices_gwas <- which(temp$stats_gwas == max_pval_gwas)
        # max_pval_indices_eqtl <- which(temp$stats_eqtl == max_pval_eqtl)
        # # str(max_pval_indices_eqtl)
        # # get the rsid (temp$markers$marker) at the index position of top marker
        # top_marker_gwas <- temp$markers$marker[max_pval_indices_gwas]
        # top_marker_eqtl <- temp$markers$marker[max_pval_indices_eqtl]
        # str(top_marker_gwas)
        # str(top_marker_eqtl)

        result_raw <- out_PPH4_rawResult[[3]]
        # # create a table of the row in result_raw that corresponds to the top marker and append H4 as a column
        # # and only keep desired columns
        # top_eqtl_table <- result_raw[result_raw$SNP.eqtl %in% top_marker_eqtl, c("GeneChr","Pvalue", "SNP.eqtl", "SNPPos", "AssessedAllele", "OtherAllele", "Zscore", "Gene", "GenePos", "NrCohorts", "NrSamples", "FDR", "BonferroniP", "AlleleA", "AlleleB", "allA_total", "allAB_total", "allB_total", "AlleleB_all")]
        # # rename GeneChr to Chr
        # colnames(top_eqtl_table)[colnames(top_eqtl_table) == "GeneChr"] <- "Chr"
        # # loop over and add .eqtl suffix to columns Pvalue, SNPPos, AssessedAllele, OtherAllele, Zscore, AlleleA, AlleleB, allA_total, allAB_total, allB_total, AlleleB_all
        # for (i in 1:ncol(top_eqtl_table)) {
        #     # skip column name that already has .eqtl suffix or if Chr
        #     if (grepl(".eqtl", colnames(top_eqtl_table)[i]) | colnames(top_eqtl_table)[i] == "Chr") {
        #         next
        #     }
        #     colnames(top_eqtl_table)[i] <- paste0(colnames(top_eqtl_table)[i], ".eqtl")
        # }

        # # queries GWAS table w/ rsid for top marker
        # top_gwas <- ieugwasr::associations(top_marker_gwas,gwas_dataset)
        # # only keep columns position, beta, se, p, n, id, rsid, ea, nea, eaf, trait
        # top_gwas <- top_gwas[, c("position", "beta", "se", "p", "n", "id", "rsid", "ea", "nea", "eaf", "trait")]
        # colnames(top_gwas) <- paste0(colnames(top_gwas), ".gwas")
        # # combine top_gwas and top_eqtl_table side-by-side in new table
        # top_table <- cbind(top_gwas, top_eqtl_table)
        # # move Chr column to first position
        # top_table <- top_table[, c("Chr", colnames(top_table)[!colnames(top_table) %in% "Chr"])]
        # top_table$H4 <- PP_H4
        top_table <- get_top_marker_raw_data(top_marker_gwas, top_marker_eqtl, result_raw, PP_H4)
        top_table_filename <- "eQTLs_colocalized_w_GWAS.txt"
        # write the header only once
        if (!file.exists(glue(eqtl_outdir, top_table_filename))) {
            write.table(top_table, glue(eqtl_outdir,top_table_filename), sep = "\t", row.names = FALSE, quote = FALSE,col.names = TRUE, append = FALSE)
        } else {
            write.table(top_table, glue(eqtl_outdir,top_table_filename), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        }    

        plot_associated_signals_and_save(temp)
        # # construct plot
        # theplot <- stack_assoc_plot_custom(temp$markers, temp$z, temp$corr, traits = temp$traits)
        
        # plot_dir <- glue("plots/{as.integer(window_size)}")
        # # output path for saving a figure for each colocalized pair of GWAS/eQTL 
        # outfile <- glue(plot_dir, "/chr{chromosome}_gwastophit{sprintf('%03d', tophit_idx)}_H4_{PP_H4}.png")
        # # if it doesn't exist, create directory
        # if (!dir.exists(plot_dir)) {
        #     dir.create(plot_dir, recursive = TRUE)
        # }
        # # save plot w/ base R functions
        # png(filename = outfile)
        # plot(theplot)
        # dev.off()

        # refresh the slick image viewer after each plot is generated
        # imgs <<- list.files(plot_dir, pattern=".png", full.names = TRUE)
        # str(imgs)
        # output[["slickr"]] <<- renderSlickR({
        # slickR(imgs) + settings(slidesToShow = 1, slidesToScroll = 1, lazyLoad = 'anticipated', dots = TRUE, arrows = TRUE, infinite = FALSE)
        # })
        
    } # close loop over each gwas top hit
    # print(glue("finished chr{chromosome}"))
} # close loop over each chromosome

clean_up()
