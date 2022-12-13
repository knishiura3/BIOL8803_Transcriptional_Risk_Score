library(arrow)

dir_output <- "data/eqtl_MAF"

# eqtl data loaded from https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz and unzipped
# load eqtl data
ds_maf <- open_dataset("data/eqtl_MAF/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt", format = "tsv")

# creates a separate parquet file for each chromosome
arrow::write_dataset(ds_maf, dir_output, partitioning = "hg19_chr", format = "parquet")

file.remove("data/eqtl_MAF/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt")