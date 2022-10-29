import time
from datetime import timedelta

import pandas as pd
import pyarrow.parquet as pq
import pyarrow as pa


def append_df_to_parquet_table(dataframe, filepath, algo, writer):

    table = pa.Table.from_pandas(dataframe)
    print(table.schema)
    if writer is None:
        writer = pq.ParquetWriter(filepath, table.schema, compression=algo)
    writer.write_table(table=table)
    return writer


def main():
    compression_algo = ("gzip","snappy")
    parquet_writer = None
    input_file = (
        "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"
    )

    for algo in compression_algo:

        start_time = time.monotonic()

        df_chunks = pd.read_csv(input_file, sep="\t", chunksize=10**6)
        for chunk in df_chunks:
            parquet_writer = append_df_to_parquet_table(
                chunk, f"eQTLs_{algo}.parquet", algo, parquet_writer
            )

        end_time = time.monotonic()
        print(f"Time taken with {algo}: {timedelta(seconds=end_time - start_time)}")
    if parquet_writer:
        parquet_writer.close()


if __name__ == "__main__":
    main()
