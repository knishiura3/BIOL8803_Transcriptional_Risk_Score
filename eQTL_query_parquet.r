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

# Reference LD Panels
#   need to set this to wherever you want the LD panel stored
ld_path <- "/home/kenji/BIOL8803/BIOL8803_Transcriptional_Risk_Score/ld_reference"

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

# this was not used,
#   1kg European reference panel for LD (legacy):
#   http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz

# coloc datasets
gwas_dataset <- "ieu-b-30"
dummy_dataset <- "ieu-a-7"
gwasinfo(id = as.character(gwas_dataset))

dir_eqtl <- "eqtls"
dir_eqtlmaf <- "eqtl_MAF"

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

# write header line to output file
# write to a log file, timestamp, chr, gwas_pos, number of eQTLs in window, PP_H4
write(paste0("Time", "\t", "chr", "\t", "pos_gwas", "\t", "num_eqtls_in_window", "\t", "PP_H4"), file = "coloc_log.txt", append = FALSE)

# define window size for coloc
window_size <- 100000

# to keep memory usage in check, work on one chromosome at a time
for (chromosome in 1:22) {
    # get the top GWAS hits for the current chromosome
    top <- ieugwasr::tophits("ieu-b-30") %>%
        filter(chr == chromosome) %>%
        arrange(p)
    print(paste0(dim(top)[1], " GWAS top hits identified on chromosome ", chromosome))


    # iterate over gwas top hits for current chromosome
    for (tophit in seq_len(nrow(top))) {
        # start timer
        tic("time elapsed for the following GWAS top hit:")

        print(top[tophit, ])

        pos_gwas <- top[tophit, ]$position
        lower <- pos_gwas - window_size
        upper <- pos_gwas + window_size
        chrpos <- paste0(chromosome, ":", lower, "-", upper)

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
            dplyr::mutate(
                beta = z / sqrt(2 * pvalues * (1 - pvalues) * (N + z**2)),
                varbeta = (1 / sqrt(2 * pvalues * (1 - pvalues) * (N + z**2)))**2, .before = snp
            ) %>%
            dplyr::mutate(pos = as.integer(pos))

        # strip name from out for compatibility with join function
        names(out) <- NULL
        out <- as.data.frame(out)

        # ensure that only loci present in both GWAS and eQTL datasets are considered
        joined <- inner_join(out, result, by = "pos", suffix = c(".gwas", ".eqtl"))

        # separate joined dataframe back out to gwas and eqtl datasets
        result <- dplyr::select(joined, ends_with("eqtl"), pos) %>%
            dplyr::rename(pvalues = pvalues.eqtl, N = N.eqtl, MAF = MAF.eqtl, beta = beta.eqtl, varbeta = varbeta.eqtl, snp = snp.eqtl, z = z.eqtl, chr = chr.eqtl)
        out <- dplyr::select(joined, ends_with("gwas"), pos, id) %>%
            dplyr::rename(pvalues = pvalues.gwas, N = N.gwas, MAF = MAF.gwas, beta = beta.gwas, varbeta = varbeta.gwas, snp = snp.gwas, z = z.gwas, chr = chr.gwas)

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
        result <- result[c("pvalues", "N", "MAF", "beta", "varbeta", "type", "snp", "z", "chr", "pos", "id")]
        result <- list(result)
        out <- list(out)

        # add name to outer list to match ieugwasr_to_coloc format
        names(result) <- "dataset2"

        # assemble gwas and eqtl datasets into one 'out' object in expected format
        out <- list(out[[1]], result[[1]])
        names(out) <- c("dataset1", "dataset2")

        # run coloc
        res <- coloc::coloc.abf(out[[1]], result[[1]])

        PP_H4 <- res$summary[["PP.H4.abf"]]

        # write to a log file, timestamp, chr, gwas_pos, number of eQTLs in window, PP_H4
        write(paste0(Sys.time(), "\t", chromosome, "\t", pos_gwas, "\t", length(out[[1]]$pos), "\t", PP_H4), file = "coloc_log.txt", append = TRUE)

        if (PP_H4 < 0.50) {
            print(glue("PP_H4 < 0.50 for eQTLs near GWAS locus chr{chromosome}:pos{pos_gwas}, skipping and continuing to next GWAS top hit"))
            next
        }
        H4 <- round(as.numeric(PP_H4), digits = 2)

        # API rejects requests if >500 rsids.
        if (length(out[[1]]$pos) >= 500) {
            # input to coloc_to_gassocplot is list of rsids (should be identical in gwas/eqtl data at this stage): out[[1]]$snp
            # choices for ancestry are AMR, AFR, EAS, EUR, SAS
            # note: a bit slow first time because the plink_bin function will download/install plink if it's not already installed.
            temp <- coloc_to_gassocplot(out, bfile = paste0(ld_path, "/EUR"), plink_bin = genetics.binaRies::get_plink_binary())
        } else {
            # query the API if <500 rsids
            temp <- coloc_to_gassocplot(out)
        }
        # show plot
        theplot <- gassocplot::stack_assoc_plot(temp$markers, temp$z, temp$corr, traits = temp$traits)

        base_dir <- glue("plots/{as.integer(window_size)}")
        # output path for saving figure
        outfile <- glue(base_dir, "/chr{chromosome}_gwas{sprintf('%03d', tophit)}_pos{pos_gwas}_H4_{H4}.png")
        # if it doesn't exist, create a directory named coloc_output
        if (!dir.exists(base_dir)) {
            dir.create(base_dir, recursive = TRUE)
        }
        # # save plot w/ base R functions
        png(filename = outfile)
        plot(theplot)
        dev.off()

        # function to save plot to file (not working:  output is blank square)
        # stack_assoc_plot_save(theplot, outfile, 2, width = 3, height = 3, dpi = 500)

        # stop timer
        toc(log = TRUE)
    } # close loop over each gwas top hit
} # close loop over each chromosome

# clean up
duckdb_unregister(con, "eqtlTable")
duckdb_unregister(con, "mafTable")
dbDisconnect(con, shutdown = TRUE)
