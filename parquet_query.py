from duckdb import query
import time
from datetime import timedelta
from numpy import mean, std

replicates = 3

runtime = tuple()
for replicate in range(replicates):
    print(f"working on replicate {replicate+1}", end="\r", flush=True)
    start_time = time.monotonic()
    query(
        """
        SELECT *
        FROM 'eQTLs_snappy.parquet'
        LIMIT 10000000
        """
    ).fetchall()
    end_time = time.monotonic()
    runtime += (timedelta.total_seconds(timedelta(seconds=end_time - start_time)),)

# print the average runtime rounded to 3 decimal points
print(
    f"\nn={replicates}\naverage (seconds):  {round(mean(runtime), 3)}\nstdev (seconds):  {round(std(runtime), 3)}"
)
