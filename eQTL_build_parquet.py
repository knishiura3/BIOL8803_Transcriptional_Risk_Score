import dask.dataframe as dd
import time
from datetime import timedelta
import sys

out_dir = "data/eqtls"


subset_cols=[
        "Pvalue",
        "SNP",
        "SNPChr",
        "SNPPos",
        "AssessedAllele",
        "OtherAllele",
        "Zscore",
        "Gene",
        "GeneChr",
        "GenePos",
        "NrCohorts",
        "NrSamples",
        "FDR",
        "BonferroniP",
    ]

eqtl_sumstats = sys.argv[1]
start_time = time.monotonic()
df = dd.read_csv(
    eqtl_sumstats,
    sep="\t",
    usecols=lambda col: col in set(subset_cols),
    dtype={
        "Pvalue": "float64",
        "SNP": "str",
        "SNPChr": "int32",
        "SNPPos": "int64",
        "AssessedAllele": "str",
        "OtherAllele": "str",
        "Zscore": "float64",
        "Gene": "str",
        "GeneChr": "str",
        "GenePos": "int64",
        "NrCohorts": "int64",
        "NrSamples": "int64",
        "FDR": "float64",
        "BonferroniP": "float64",
    },
)

df.repartition(partition_size="1024MB").to_parquet(
    out_dir,
    partition_on=["SNPChr"],
    engine="pyarrow",
    write_index=False,
    compression="snappy",
    append=True
)
end_time = time.monotonic()
print(f"Time taken: {timedelta(seconds=end_time - start_time)}")
