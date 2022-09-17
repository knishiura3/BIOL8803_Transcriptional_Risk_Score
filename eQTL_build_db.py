print("Usage:")
print("python3 eQTL_build_db.py <eQTL input .gz file> <output DB name>")


import sqlite3
import gzip
from sys import argv


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

    # convert eQTL file to SQL database
    def eQTL_to_sql(self, input_file_eQTL, separator="\t"):
        # Create default table
        table_sql = "CREATE TABLE IF NOT EXISTS eqtlTable (Pvalue REAL, SNP TEXT, SNPChr INTEGER, SNPPos INTEGER, AssessedAllele TEXT, OtherAllele TEXT, Zscore REAL, Gene TEXT, GeneSymbol TEXT, GeneChr INTEGER, GenePos INTEGER, NrCohorts INTEGER, NrSamples INTEGER, FDR REAL, BonferroniP REAL)"
        self.cursor.execute(table_sql)
        self.connection.commit()

        fh = gzip.open(input_file_eQTL, "rt")

        # skip header
        fh.readline()

        count = 0
        current_lines = []

        for line in fh:

            segs = line.strip().split(separator)

            Pvalue = float(segs[0])
            SNP = str(segs[1])
            SNPChr = int(segs[2])
            SNPPos = int(segs[3])
            AssessedAllele = str(segs[4])
            OtherAllele = str(segs[5])
            Zscore = float(segs[6])
            Gene = str(segs[7])
            GeneSymbol = str(segs[8])
            GeneChr = int(segs[9])
            GenePos = int(segs[10])
            NrCohorts = int(segs[11])
            NrSamples = int(segs[12])
            FDR = float(segs[13])
            BonferroniP = float(segs[14])

            row = (
                Pvalue,
                SNP,
                SNPChr,
                SNPPos,
                AssessedAllele,
                OtherAllele,
                Zscore,
                Gene,
                GeneSymbol,
                GeneChr,
                GenePos,
                NrCohorts,
                NrSamples,
                FDR,
                BonferroniP,
            )

            current_lines.append(row)
            count += 1

            if len(current_lines) >= 10**6:
                self.cursor.executemany(
                    "INSERT OR REPLACE INTO eqtlTable VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    current_lines,
                )
                self.connection.commit()
                count = 0
                current_lines = []

        fh.close()

        if len(current_lines) >= 0:
            self.cursor.executemany(
                "INSERT OR REPLACE INTO eqtlTable VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                current_lines,
            )
            self.connection.commit()
            count = 0
            current_lines = []

        print(f"Finished loading {input_file_eQTL} into {self.db}. Indexing...")

        # create index to speed up filtering
        self.cursor.execute(
            "CREATE INDEX IF NOT EXISTS pos_index ON eqtlTable (SNPChr, SNPPos)"
        )
        self.connection.commit()


def main():

    # take input/output variables from command line arguments, otherwise use hardcoded defaults
    input_file_eqtl = (
        argv[1]
        if len(argv) > 1
        else "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
    )
    db_name = argv[2] if len(argv) > 2 else "eQTLs.db"

    manager = eqtl_DB(db_name)
    manager.connect()
    # Note:  this takes ~30 minutes for full summary statistics file, DB takes up 13.7 GB
    # but only takes a few minutes for FDR filtered file, DB takes up 1.1 GB
    manager.eQTL_to_sql(input_file_eqtl)

    manager.close_connection()


if __name__ == "__main__":
    main()
