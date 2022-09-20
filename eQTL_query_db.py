print("Usage:")
print(
    "python3 eQTL_query_db.py <GWAS hit file> <eQTL DB> <Window size (bp)> <output directory>"
)


import sqlite3
import datetime
import os
from sys import argv
from math import sqrt

from alive_progress import alive_bar


class eqtl_DB:
    def __init__(self, db_name):
        # Create a database
        self.db = db_name
        self.connection = None
        self.cursor = None

    def connect(self):
        # Connect to the database and create a cursor
        self.connection = sqlite3.connect(self.db)
        self.cursor = self.connection.cursor()

    def close_connection(self):
        # Clean up and close the connection
        self.cursor.close()
        self.cursor = None
        self.connection.close()
        self.connection = None

    # Query the database for all rows within 50kb of given gwas position
    def query_interval(self, chromosome, lower_end, upper_end):
        # note:  BETWEEN operator is inclusive of both endpoints, sort by rsid
        query = "SELECT * FROM eqtlTable WHERE SNPChr = ? AND SNPPos BETWEEN ? AND ? ORDER BY SNP, Pvalue"
        self.cursor.execute(query, (chromosome, lower_end, upper_end))
        return self.cursor.fetchall()

    # query the database for specific SNP by chromosome and position
    def query_position(self, chromosome, position):
        # note:  BETWEEN operator is inclusive of both endpoints, sort by rsid
        query = "SELECT * FROM eqtlTable WHERE SNPChr = ? AND SNPPos = ? ORDER BY SNP, Pvalue"
        self.cursor.execute(query, (chromosome, position))
        # change this to fetchone() if data is validated to only have no duplicates
        return self.cursor.fetchall()


def main():
    # take input/output variables from command line arguments, if not all provided, use hardcoded defaults
    input_gwas_dir = str(argv[1]) if len(argv) > 1 else "gwas_hits"
    # input_file_gwas = str(argv[1]) if len(argv) > 1 else "gwas_top_hits.tsv"
    db_name = str(argv[2]) if len(argv) > 2 else "eQTLs_full.db"
    # window = int(argv[3]) if len(argv) > 3 else int(50000)
    # output_dir = str(argv[4]) if len(argv) > 4 else "output_full"
    output_dir = str(argv[3]) if len(argv) > 3 else "output_full"

    manager = eqtl_DB(db_name)
    manager.connect()

    # loop over the files in the directory with SNPs of interest near GWAS top hits
    for dirpath, dirnames, filenames in os.walk(input_gwas_dir):
        # getting total number of gwas top hits to use for progress bar
        total_regions = len(filenames)
        # initialize progress bar
        with alive_bar(total_regions) as bar:
            for filename in filenames:
                if filename.endswith("_gwas.tsv"):
                    input_file_gwas = os.path.join(dirpath, filename)
                    (
                        pval_rank_gwashit,
                        chr_gwas_tophit,
                        lower_pos_gwas_tophit,
                        upper_pos_gwas_tophit,
                        _,
                    ) = filename.split("_")
                    # declare unpacked variables as integers
                    # add leading 0 to pval_rank_gwashit so that it will always be 3 digits for output file name
                    pval_rank_gwashit_str = str(pval_rank_gwashit).zfill(3)
                    pval_rank_gwashit_int = int(pval_rank_gwashit)
                    chr_gwas_tophit = int(chr_gwas_tophit)
                    lower_pos_gwas_tophit = int(lower_pos_gwas_tophit)
                    upper_pos_gwas_tophit = int(upper_pos_gwas_tophit)

                    # print("Processing file: " + input_file_gwas)

                    # with open(input_file_gwas, "r") as fh:
                    #     # count total lines in file for progress bar
                    #     total_lines = sum(1 for line in fh)

                    with open(input_file_gwas, "r") as fh:
                        total_eQTL_count = 0
                        header_check = None
                        # skip header line
                        next(fh)
                        # loop over each SNP in the gwas file
                        for line in fh:
                            segs = line.strip().split("\t")
                            # pvalues_gwas = float(segs[0])
                            # N_gwas = int(segs[1])
                            MAF_gwas = float(segs[2])
                            # beta_gwas = float(segs[3])
                            # varbeta_gwas = float(segs[4])
                            # type_gwas = str(segs[5])
                            # snp_gwas = str(segs[6])
                            # z_gwas = float(segs[7])
                            chr_gwas = int(segs[8])
                            pos_gwas = int(segs[9])
                            # id_gwas = str(segs[10])

                            # capture time info for logging
                            current_date = (
                                f"{datetime.datetime.now().strftime('%Y-%m-%d')}"
                            )
                            current_date_and_time = f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"

                            # query the eQTL DB for eQTLs in the exact position of GWAS row
                            results = manager.query_position(chr_gwas, pos_gwas)

                            # query the eQTL DB for all eQTLs falling within the window around a SNP locus and save to a list
                            # results = manager.query_interval(
                            #     chr_gwas_tophit,
                            #     lower_pos_gwas_tophit,
                            #     upper_pos_gwas_tophit,
                            # )

                            # create an output directory if it doesn't exist already
                            if not os.path.exists(f"{output_dir}"):
                                os.makedirs(f"{output_dir}")
                            # write eQTLs to tab-delimited file in output directory
                            with open(
                                f"{output_dir}/{pval_rank_gwashit_str}_{chr_gwas_tophit}_{lower_pos_gwas_tophit}_{upper_pos_gwas_tophit}_eQTLs.tsv",
                                "a",
                            ) as out:
                                # checks if header has been written to file yet
                                if header_check is None:
                                    header = f"pvalues\tN\tMAF\tbeta\tvarbeta\ttype\tsnp\tz\tchr\tpos\tid\n"
                                    out.write(header)
                                    header_check = True
                                # total_eQTL_count = 0
                                SNP_set = set()
                                output_list = []
                                # looping over eQTL DB query results
                                for result in results:
                                    # assign result to variables
                                    Pvalue = float(result[0])
                                    # discard eQTL rows where p-value = 1
                                    if Pvalue == 1:
                                        continue
                                    SNP = str(result[1])
                                    SNPChr = int(result[2])
                                    SNPPos = int(result[3])
                                    # AssessedAllele = str(result[4])
                                    # OtherAllele = str(result[5])
                                    Zscore = float(result[6])
                                    # Gene = str(result[7])
                                    # GeneSymbol = str(result[8])
                                    # GeneChr = int(result[9])
                                    # GenePos = int(result[10])
                                    # NrCohorts = int(result[11])
                                    NrSamples = int(result[12])
                                    # FDR = float(result[13])
                                    # BonferroniP = float(result[14])

                                    # calculated variables
                                    # source:  Zhu, Z. et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets.
                                    # betaSE <- 1/sqrt(2*p*(1-p)*(n+z^2))
                                    # beta <- z/sqrt(2*p*(1-p)*(n+z^2))
                                    # try calculating betaSE and if there's an error, print the values of Pvalue, NrSamples, and Zscore
                                    try:
                                        betaSE = 1 / sqrt(
                                            2
                                            * Pvalue
                                            * (1 - Pvalue)
                                            * (NrSamples + Zscore**2)
                                        )
                                        beta = Zscore / sqrt(
                                            2
                                            * Pvalue
                                            * (1 - Pvalue)
                                            * (NrSamples + Zscore**2)
                                        )
                                    except:
                                        print(
                                            f"Error calculating beta/betaSE:  Pvalue={Pvalue} NrSamples={NrSamples} Zscore={Zscore}"
                                        )
                                    # set MAF to eaf from gwas hits, if not defined, set to 0.5
                                    MAF_eqtl = MAF_gwas if MAF_gwas != "NA" else 0.5

                                    # variance of beta is simply standard error squared
                                    varbeta = betaSE**2

                                    type = "quant"

                                    # since SQL query results are already sorted by SNP and p-value,
                                    # we only need to keep the first row for each unique SNP, add to list, sort list by SNP,
                                    # and then loop through list to write to file

                                    if SNP not in SNP_set:
                                        SNP_set.add(SNP)
                                        output_list.append(
                                            [
                                                Pvalue,
                                                NrSamples,
                                                MAF_eqtl,
                                                beta,
                                                varbeta,
                                                type,
                                                SNP,
                                                Zscore,
                                                SNPChr,
                                                SNPPos,
                                                db_name.split(".")[0],
                                                "\n",
                                            ]
                                        )

                                        total_eQTL_count += 1
                                    else:
                                        continue

                                    # export raw eQTL rows
                                    # out.write(
                                    #     f"{result[0]}\t{result[1]}\t{result[2]}\t{result[3]}\t{result[4]}\t{result[5]}\t{result[6]}\t{result[7]}\t{result[8]}\t{result[9]}\t{result[10]}\t{result[11]}\t{result[12]}\t{result[13]}\t{result[14]}\n"
                                    # )

                                # list comprehension to write tab-separated, terminating each line w/ newline
                                out.write(
                                    "\n".join(
                                        "\t".join(map(str, row)) for row in output_list
                                    )
                                )
                        # increment progressbar counter
                        bar()
                        # write to log
                        with open(
                            f"{current_date}_gwas_hits_in_eqtl_regions.log",
                            "a",
                        ) as logger:
                            logger.write(
                                # f"{current_date_and_time} {total_eQTL_count} eQTLs in region +/- {window}bp of gwas hit (rsid chr pos):  {rsid} {chr} {position}\n"
                                f"{current_date_and_time} {total_eQTL_count} eQTLs in region : #{pval_rank_gwashit_int} GWAS top hit at chr {chr_gwas_tophit} in region {lower_pos_gwas_tophit}:{upper_pos_gwas_tophit}\n"
                            )

    manager.close_connection()


if __name__ == "__main__":
    main()
