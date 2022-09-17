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
    # Return a list of tuples
    def query_interval(self, chromosome, position, window):
        # note:  BETWEEN operator is inclusive of both endpoints, sort by rsid
        query = "SELECT * FROM eqtlTable WHERE SNPChr = ? AND SNPPos BETWEEN ? AND ? ORDER BY SNP, Pvalue"
        self.cursor.execute(query, (chromosome, position - window, position + window))
        return self.cursor.fetchall()


def main():
    # take input/output variables from command line arguments, if not all provided, use hardcoded defaults
    input_file_gwas = str(argv[1]) if len(argv) > 1 else "gwas_top_hits.tsv"
    db_name = str(argv[2]) if len(argv) > 2 else "eQTLs_fdr.db"
    window = int(argv[3]) if len(argv) > 3 else int(50000)
    output_dir = str(argv[4]) if len(argv) > 4 else "output_fdr"

    manager = eqtl_DB(db_name)
    manager.connect()

    # count total lines in file for progress bar
    with open(input_file_gwas, "r") as fh:
        for count, line in enumerate(fh):
            pass
        line_total = count + 1

    # initialize progress bar
    with alive_bar(line_total) as bar:
        # loop over each gwas hit and query the eQTL DB for all eQTLs within a certain window
        with open(input_file_gwas, "r") as fh:

            for line in fh:

                segs = line.strip().split("\t")

                chr = int(segs[0])
                position = int(segs[1])
                beta = float(segs[2])
                se = float(segs[3])
                varbeta_gwas = se**2
                p = float(segs[4])
                n = int(segs[5])
                id = str(segs[6])
                rsid = str(segs[7])
                ea = str(segs[8])
                nea = str(segs[9])
                eaf = float(segs[10])
                trait = str(segs[11])

                # calculate zscore from beta, se, and p-value
                zscore_gwas = beta / se

                # type? not sure about this one
                type_gwas = "quant"

                # write counts to a log file in current directory
                current_date = f"{datetime.datetime.now().strftime('%Y-%m-%d')}"
                current_date_and_time = (
                    f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
                )

                # create an output directory if it doesn't exist already
                if not os.path.exists(f"{output_dir}"):
                    os.makedirs(f"{output_dir}")

                # write GWAS hit to tab-delimited file in output directory
                with open(
                    f"{output_dir}/region_{rsid}_{chr}_{position}_GWAS.tsv", "w"
                ) as out:
                    header = (
                        f"pvalues\tN\tMAF\tbeta\tvarbeta\ttype\tsnp\tz\tchr\tpos\tid\n"
                    )
                    out.write(header)
                    out.write(
                        f"{p}\t{n}\t{eaf}\t{beta}\t{varbeta_gwas}\t{type_gwas}\t{rsid}\t{zscore_gwas}\t{chr}\t{position}\tieu-b-30\n"
                    )

                # query the eQTL DB for all eQTLs falling within the window around a GWAS hit locus and save to a list
                results = manager.query_interval(chr, position, window)

                # write eQTLs to tab-delimited file in output directory
                with open(
                    f"{output_dir}/region_{rsid}_{chr}_{position}_eQTLs.tsv", "w"
                ) as out:
                    header = (
                        f"pvalues\tN\tMAF\tbeta\tvarbeta\ttype\tsnp\tz\tchr\tpos\tid\n"
                    )
                    out.write(header)
                    total_eQTL_count = 0
                    SNP_set = set()
                    for result in results:
                        # assign result to variables
                        Pvalue = float(result[0])
                        # discard eQTL rows where p-value = 1
                        if Pvalue == 1:
                            continue
                        SNP = str(result[1])
                        SNPChr = int(result[2])
                        SNPPos = int(result[3])
                        AssessedAllele = str(result[4])
                        OtherAllele = str(result[5])
                        Zscore = float(result[6])
                        Gene = str(result[7])
                        GeneSymbol = str(result[8])
                        GeneChr = int(result[9])
                        GenePos = int(result[10])
                        NrCohorts = int(result[11])
                        NrSamples = int(result[12])
                        FDR = float(result[13])
                        BonferroniP = float(result[14])

                        # calculated variables
                        # source:  Zhu, Z. et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets.
                        # betaSE <- 1/sqrt(2*p*(1-p)*(n+z^2))
                        # beta <- z/sqrt(2*p*(1-p)*(n+z^2))
                        # try calculating betaSE and if there's an error, print the values of Pvalue, NrSamples, and Zscore
                        try:
                            betaSE = 1 / sqrt(
                                2 * Pvalue * (1 - Pvalue) * (NrSamples + Zscore**2)
                            )
                            beta = Zscore / sqrt(
                                2 * Pvalue * (1 - Pvalue) * (NrSamples + Zscore**2)
                            )
                        except:
                            print(
                                f"Error calculating beta/betaSE:  Pvalue={Pvalue} NrSamples={NrSamples} Zscore={Zscore}"
                            )
                        # set MAF to eaf from gwas hits, if not defined, set to 0.5
                        MAF = eaf if eaf != "NA" else 0.5

                        # variance of beta is simply standard error squared
                        varbeta = betaSE**2

                        # type? not sure about this one
                        type = "quant"

                        # since SQL query results are already sorted by SNP and p-value, we only need to keep the first row for each unique SNP
                        if SNP not in SNP_set:
                            SNP_set.add(SNP)
                            out.write(
                                f"{Pvalue}\t{NrSamples}\t{MAF}\t{beta}\t{varbeta}\t{type}\t{SNP}\t{Zscore}\t{SNPChr}\t{SNPPos}\t{db_name.split('.')[0]}\n"
                            )
                            total_eQTL_count += 1
                        else:
                            continue

                        # export raw eQTL rows
                        # out.write(
                        #     f"{result[0]}\t{result[1]}\t{result[2]}\t{result[3]}\t{result[4]}\t{result[5]}\t{result[6]}\t{result[7]}\t{result[8]}\t{result[9]}\t{result[10]}\t{result[11]}\t{result[12]}\t{result[13]}\t{result[14]}\n"
                        # )

                with open(
                    f"{current_date}_gwas_hits_in_eqtl_regions.log",
                    "a",
                ) as logger:
                    logger.write(
                        f"{current_date_and_time} {total_eQTL_count} eQTLs in region +/- {window}bp of gwas hit (rsid chr pos):  {rsid} {chr} {position}\n"
                    )

                # update progress bar
                bar()

    manager.close_connection()


if __name__ == "__main__":
    main()
