print("Usage:")
print("python3 eQTL_in_GWAS_hit_regions.py <GWAS hit file> <eQTL file>")


import sqlite3
import datetime
import os

from alive_progress import alive_bar
import time

# import pandas as pd
# import numpy as np


class eqtl_DB:
    def __init__(self, db_name):
        # Create a database
        self.db = db_name
        self.connection = None
        self.cursor = None

        self.taxa = None

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
        # note:  BETWEEN operator is inclusive of both endpoints
        query = "SELECT * FROM eqtlTable WHERE SNPChr = ? AND SNPPos BETWEEN ? AND ?"
        self.cursor.execute(query, (chromosome, position - window, position + window))
        return self.cursor.fetchall()


# hardcoded test files
db_name = "eQTLs.db"
input_file_gwas = "gwas_hits.txt"

# size of interval
window = 50000

manager = eqtl_DB(db_name)
manager.connect()

# count total lines in file for progress bar
with open(input_file_gwas, "r") as fh:
    for count, line in enumerate(fh):
        pass
    line_total = count + 1

# loop over each gwas hit and query the eQTL DB for all eQTLs within a certain window
with open(input_file_gwas, "r") as fh:
    # initialize progress bar
    with alive_bar(line_total) as bar:

        for line in fh:

            segs = line.strip().split(" ")

            rsid = str(segs[0])
            chromosome = int(segs[1])
            position = int(segs[2])

            results = manager.query_interval(chromosome, position, window)

            # write counts to a log file
            with open(f"gwas_hits_in_eqtl_regions.log", "a") as logger:
                logger.write(
                    f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {len(results)} eQTLs in region +/- {window}bp of gwas hit (rsid chr pos):  {line}"
                )

            # create an output directory if it doesn't exist already
            if not os.path.exists("output"):
                os.makedirs("output")

            # write results to file tab-delimited
            with open(
                f"output/eQTLs_in_region_{rsid}_{chromosome}_{position}.tsv", "w"
            ) as out:
                for result in results:
                    out.write(
                        f"{result[0]}\t{result[1]}\t{result[2]}\t{result[3]}\t{result[4]}\t{result[5]}\t{result[6]}\t{result[7]}\t{result[8]}\t{result[9]}\t{result[10]}\t{result[11]}\t{result[12]}\t{result[13]}\t{result[14]}\n"
                    )
            # update progress bar
            bar()

manager.close_connection()
