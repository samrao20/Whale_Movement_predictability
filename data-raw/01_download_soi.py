"""Download the monthly Southern Oscillation Index (NOAA CPC) -> data/soi_monthly.csv

The CPC file has two tables: raw Tahiti-Darwin pressure anomaly, then the
STANDARDIZED data. The standardized one is the conventional SOI and is used here.

Run from the repo root:  python data-raw/01_download_soi.py
"""
import pandas as pd
import requests

URL = "https://www.cpc.ncep.noaa.gov/data/indices/soi"

lines = requests.get(URL, timeout=60).text.splitlines()
start = next(i for i, l in enumerate(lines) if "STANDARDIZED" in l)

recs = []
for line in lines[start + 1:]:
    if len(line) >= 4 and line[:4].isdigit():
        year = int(line[:4])
        # fixed width: 4-char year then 12 x 6-char values (they can touch)
        for m in range(12):
            v = float(line[4 + 6 * m: 10 + 6 * m])
            if v > -900:  # -999.9 = missing
                recs.append((year, m + 1, v))

df = pd.DataFrame(recs, columns=["year", "month", "soi"])
df.to_csv("data/soi_monthly.csv", index=False)
print(df.tail(), len(df))
