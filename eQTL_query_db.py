import sqlite3
import datetime
import os
import click
import pandas as pd
import numpy as np
from math import sqrt
import argparse

from alive_progress import alive_bar


class eqtl_DB:
    def __init__(self, db_name):
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


# @click.command()
# @click.argument("input_gwas_dir", type=click.Path(exists=True), default="gwas_hits")
# @click.argument("db_name", type=click.Path(exists=True), default="eQTLs.db")
# @click.argument("output_dir", type=click.Path(), default="output")
def main(input_file, db_name, output_dir):

    manager = eqtl_DB(db_name)
    manager.connect()

    filename = input_file
    if filename.endswith("_gwas.tsv"):
        input_file_gwas = filename
        filename_suffix = filename.split("/")[1]
        (
            pval_rank_gwashit,
            chr_gwas_tophit,
            lower_pos_gwas_tophit,
            upper_pos_gwas_tophit,
            _,
        ) = filename_suffix.split("_")
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
            output_list = []
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
                current_date = f"{datetime.datetime.now().strftime('%Y-%m-%d')}"
                current_date_and_time = (
                    f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
                )

                # query the eQTL DB for eQTLs in the exact position of GWAS row
                try:
                    results = [manager.query_position(chr_gwas, pos_gwas)[0]]
                except IndexError:
                    print(f"no eQTLs found for {chr_gwas}:{pos_gwas}")
                    continue

                # read manager.query_position(chr_gwas, pos_gwas into dataframe, only keep columns used for calculations
                df = pd.DataFrame(
                    results,
                    columns=[
                        "Pvalue",
                        "SNP",
                        "SNPChr",
                        "SNPPos",
                        "AssessedAllele",
                        "OtherAllele",
                        "Zscore",
                        "Gene",
                        "GeneSymbol",
                        "GeneChr",
                        "GenePos",
                        "NrCohorts",
                        "NrSamples",
                        "FDR",
                        "BonferroniP",
                    ],
                )
                # Only keep columns Pvalue, SNP, SNPChr, SNPPos, Zscore, NrSamples
                df = df[
                    [
                        "Pvalue",
                        "SNP",
                        "SNPChr",
                        "SNPPos",
                        "Zscore",
                        "NrSamples",
                    ]
                ]
                # drop rows where Pvalue = 1
                # df = df[df.Pvalue != 1]
                # vectorized calculation of betaSE:  betaSE = 1 / sqrt(2 * Pvalue * (1 - Pvalue) * (NrSamples + Zscore**2))
                df["betaSE"] = 1 / np.sqrt(
                    2
                    * df["Pvalue"]
                    * (1 - df["Pvalue"])
                    * (df["NrSamples"] + df["Zscore"] ** 2)
                )
                # vectorized calculation of beta:  beta = Zscore / sqrt(2 * Pvalue * (1 - Pvalue) * (NrSamples + Zscore**2))
                df["beta"] = df["Zscore"] / np.sqrt(
                    2
                    * df["Pvalue"]
                    * (1 - df["Pvalue"])
                    * (df["NrSamples"] + df["Zscore"] ** 2)
                )
                # vectorized assignment of MAF_eqtl = MAF_gwas if MAF_gwas != "NA" else 0.5
                df["MAF_eqtl"] = MAF_gwas if MAF_gwas != "NA" else 0.5

                # vectorized calculation of varbeta:  betaSE**2
                df["varbeta"] = df["betaSE"] ** 2
                df["type"] = "quant"
                df["db_name"] = db_name.split(".")[0]
                # reorder columns in desired output format
                df = df[
                    [
                        "Pvalue",
                        "NrSamples",
                        "MAF_eqtl",
                        "beta",
                        "varbeta",
                        "type",
                        "SNP",
                        "Zscore",
                        "SNPChr",
                        "SNPPos",
                        "db_name",
                    ]
                ]
                # append df["Pvalue"] to output list without a header
                output_list.append(df.values.tolist()[0])

            # increment progressbar counter
            # bar()
            # write header and output_list to file
            with open(
                f"{output_dir}/{pval_rank_gwashit_str}_{chr_gwas_tophit}_{lower_pos_gwas_tophit}_{upper_pos_gwas_tophit}_eQTLs.tsv",
                "w",
            ) as out:
                header = [
                    "Pvalue",
                    "NrSamples",
                    "MAF_eqtl",
                    "beta",
                    "varbeta",
                    "type",
                    "SNP",
                    "Zscore",
                    "SNPChr",
                    "SNPPos",
                    "db_name",
                ]
                out.write("\t".join(header) + "\n")
                out.write("\n".join("\t".join(map(str, row)) for row in output_list))
            # write to log
            total_eQTL_count = len(output_list)
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
    # @click.argument("input_gwas_dir", type=click.Path(exists=True), default="gwas_hits")
    # @click.argument("db_name", type=click.Path(exists=True), default="eQTLs.db")
    # @click.argument("output_dir", type=click.Path(), default="output")
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True)
    parser.add_argument("--db_name", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    args = parser.parse_args()
    main(args.input_file, args.db_name, args.output_dir)

# python3 eQTL_build_db-2.py --input_file {} --db_name eQTLs.db --output_dir output
