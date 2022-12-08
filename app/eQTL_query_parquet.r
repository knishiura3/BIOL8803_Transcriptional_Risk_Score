# packages for parquet querying
library(arrow)
library(duckdb)
library(fs)
library(tidyverse)
library(DBI)
library(glue)

# packages for coloc
suppressPackageStartupMessages(suppressWarnings({
    library(gwasglue)
    library(gassocplot)
    library(dplyr)
    library(coloc)
    library(ggplot2)
    library(httpgd)

}))


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

gather_and_format_gwas_eqtl_in_region <- function(con, chromosome, gwas_hit) {
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

plot_associated_signals_and_save <- function(tophit_idx, chromosome, PP_H4, temp, plot_dir) {
    # this is needed to suppress automatic pdf output from gassocplot
    pdf(NULL)
    # construct plot
    theplot <- stack_assoc_plot_custom(temp$markers, temp$z, temp$corr, traits = temp$traits)
    # output path for saving a figure for each colocalized pair of GWAS/eQTL 
    outfile <- glue::glue(plot_dir, "/chr{chromosome}_gwastophit{sprintf('%03d', tophit_idx)}_H4_{PP_H4}.png")
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
clean_up <- function(con) {
    # clean up
    print("cleaning up")
    duckdb_unregister(con, "eqtlTable")
    duckdb_unregister(con, "mafTable")
    dbDisconnect(con, shutdown = TRUE)
    print("done")
}
