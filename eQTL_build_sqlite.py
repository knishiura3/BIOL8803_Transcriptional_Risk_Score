import sqlite3
import gzip
import click
import time
from datetime import timedelta


class eqtl_DB:
    def __init__(self, db_name):
        self.db = db_name
        self.connection = None
        self.cursor = None

    def connect(self):
        # Connect to the database and create a cursor
        self.connection = sqlite3.connect(self.db)
        self.cursor = self.connection.cursor()

    def adjust_settings(self):
        # self.cursor.execute("PRAGMA page_size = 32768;")
        # self.connection.commit()
        self.cursor.execute("PRAGMA journal_mode = WAL;")
        self.connection.commit()
        self.cursor.execute("PRAGMA synchronous = NORMAL;")
        self.connection.commit()
        # self.cursor.execute("PRAGMA locking_mode = exclusive;")
        # self.connection.commit()

    def close_connection(self):
        # Clean up and close the connection
        self.cursor.close()
        self.cursor = None
        self.connection.close()
        self.connection = None

    # convert eQTL file to SQL database
    def eQTL_to_sql(self, input_file_eQTL, separator="\t"):

        table_sql = "CREATE TABLE IF NOT EXISTS eqtlTable (Pvalue REAL, SNP TEXT, SNPChr INTEGER, SNPPos INTEGER, Zscore REAL, NrSamples INTEGER)"
        # table_sql = "CREATE TABLE IF NOT EXISTS eqtlTable (Pvalue REAL, SNP TEXT, SNPChr INTEGER, SNPPos INTEGER, AssessedAllele TEXT, OtherAllele TEXT, Zscore REAL, Gene TEXT, GeneSymbol TEXT, GeneChr INTEGER, GenePos INTEGER, NrCohorts INTEGER, NrSamples INTEGER, FDR REAL, BonferroniP REAL)"
        self.cursor.execute(table_sql)
        self.connection.commit()

        # if the input file ends with .gz, use gzip.open. otherwise, just use open
        if input_file_eQTL.endswith(".gz"):
            fh = gzip.open(input_file_eQTL, "rt")
        else:
            fh = open(input_file_eQTL, "rt")

        # skip header
        fh.readline()

        count = 0
        current_lines = []

        for line in fh:

            segs = line.strip().split(separator)

            # only build DB with subset of columns since we don't need all of them
            Pvalue = float(segs[0])
            SNP = str(segs[1])
            SNPChr = int(segs[2])
            SNPPos = int(segs[3])
            # AssessedAllele = str(segs[4])
            # OtherAllele = str(segs[5])
            Zscore = float(segs[6])
            # Gene = str(segs[7])
            # GeneSymbol = str(segs[8])
            # GeneChr = int(segs[9])
            # GenePos = int(segs[10])
            # NrCohorts = int(segs[11])
            NrSamples = int(segs[12])
            # FDR = float(segs[13])
            # BonferroniP = float(segs[14])

            row = (
                float(segs[0]),  # Pvalue
                str(segs[1]),  # SNP
                int(segs[2]),  # SNPChr
                int(segs[3]),  # SNPPos
                # AssessedAllele,
                # OtherAllele,
                float(segs[6]),  # Zscore
                # Gene,
                # GeneSymbol,
                # GeneChr,
                # GenePos,
                # NrCohorts,
                int(segs[12]),  # NrSamples
                # FDR,
                # BonferroniP,
            )

            current_lines.append(row)
            count += 1
            # print count every 100000 lines
            if count % 100000 == 0:
                print(f"Inserted {count} rows into eqtlTable", end="\r", flush=True)

            if len(current_lines) >= 10**7:
                self.cursor.executemany(
                    "INSERT OR REPLACE INTO eqtlTable VALUES (?, ?, ?, ?, ?, ?)",
                    # "INSERT OR REPLACE INTO eqtlTable VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    current_lines,
                )
                self.connection.commit()
                count = 0
                # clear printed progress
                print(end="\x1b[2K", flush=True)
                current_lines = []

        fh.close()

        if len(current_lines) >= 0:
            self.cursor.executemany(
                "INSERT OR REPLACE INTO eqtlTable VALUES (?, ?, ?, ?, ?, ?)",
                # "INSERT OR REPLACE INTO eqtlTable VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
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

    def eQTL_MAF_to_sql(self, input_file_MAF, separator="\t"):

        table_sql = "CREATE TABLE IF NOT EXISTS mafTable (SNP TEXT, chr INTEGER, pos INTEGER, maf REAL, PRIMARY KEY (chr, pos))"
        self.cursor.execute(table_sql)
        self.connection.commit()

        fh = open(input_file_MAF, "rt")

        # skip header
        fh.readline()

        count = 0
        current_lines = []

        for line in fh:

            segs = line.strip().split(separator)

            # if any element of segs is 'NA', skip this row
            if any([x == "NA" for x in segs]):
                continue

            # SNP = str(segs[0])
            # chr = int(segs[1])
            # pos = int(segs[2])
            # maf = float(segs[8])

            row = (
                str(segs[0]),
                int(segs[1]),
                int(segs[2]),
                float(segs[8]),
            )

            current_lines.append(row)
            count += 1
            # print count every 100000 lines
            if count % 100000 == 0:
                print(f"Inserted {count} rows into mafTable", end="\r", flush=True)

            if len(current_lines) >= 10**7:
                self.cursor.executemany(
                    "INSERT OR REPLACE INTO mafTable VALUES (?, ?, ?, ?)",
                    current_lines,
                )
                self.connection.commit()
                count = 0
                # clear printed progress
                print(end="\x1b[2K", flush=True)
                current_lines = []

        fh.close()

        if len(current_lines) >= 0:
            self.cursor.executemany(
                "INSERT OR REPLACE INTO mafTable VALUES (?, ?, ?, ?)",
                current_lines,
            )
            self.connection.commit()
            count = 0
            current_lines = []

        print(f"Finished loading {input_file_MAF} into {self.db}")
        # create index to speed up filtering
        # self.cursor.execute(
        #     "CREATE INDEX IF NOT EXISTS rsid_index ON mafTable (chr, pos)"
        # )
        # self.connection.commit()


# use click to t
@click.command()
@click.option(
    "--input_eqtl",
    "-i",
    type=click.Path(exists=True),
    help="Input eQTL file",
    default="2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
    # default="subset.txt.gz",
)
@click.option(
    "--output_db_name",
    "-o",
    type=click.Path(),
    help="Output database name",
    default="eQTLs.db",
)
@click.option(
    "--input_maf",
    "-m",
    type=click.Path(exists=True),
    help="Input MAF file",
    default="2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt",
)
def main(input_eqtl, input_maf, output_db_name):
    manager = eqtl_DB(output_db_name)
    manager.connect()
    manager.adjust_settings()

    eqtl_start = time.monotonic()
    manager.eQTL_to_sql(input_eqtl)
    eqtl_end = time.monotonic()
    print(
        f"Finished loading eQTLs into sqlite db in {timedelta(seconds=eqtl_end - eqtl_start)} h:mm:ss"
    )

    delete_start = time.monotonic()
    manager.delete_duplicate_rsid_in_eqtlTable_keep_lowest_Pvalue()
    delete_end = time.monotonic()
    print(
        f"Finished deleting duplicate rsids in {timedelta(seconds=delete_end - delete_start)} h:mm:ss"
    )

    maf_start = time.monotonic()
    manager.eQTL_MAF_to_sql(input_maf)
    maf_end = time.monotonic()
    print(
        f"Finished loading MAFs into sqlite db in {timedelta(seconds=maf_end - maf_start)} h:mm:ss"
    )

    manager.close_connection()


if __name__ == "__main__":
    main()
