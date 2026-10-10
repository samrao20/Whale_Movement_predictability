"""IMOS OceanCurrent GSLA delayed-mode geostrophic currents -> data/currents_cells.csv

Samples every 5th day (~6 per month) via OpenDAP, averages u/v onto the 2-degree
cells of the Water sheet, then takes monthly means. Resumable: per-year partial
results are cached in data-raw/cache/ (git-ignored).

Columns: u_ms (east), v_ms (north; negative = southward, e.g. East Australian
Current), speed_ms, eke (mean of 0.5*(u'^2+v'^2) about the monthly mean, a
mesoscale-eddy-activity proxy), n_days.

Run from repo root:  python data-raw/03_download_currents.py [first_year last_year]
Needs: pandas, xarray, netCDF4, requests
"""
import os, re, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
import requests
import xarray as xr

CAT = "https://thredds.aodn.org.au/thredds/catalog/IMOS/OceanCurrent/GSLA/DM/{y}/catalog.xml"
DAP = "https://thredds.aodn.org.au/thredds/dodsC/"
XLSX = "data-raw/Humpback_East_Australia_data.xlsx"
CACHE = "data-raw/cache"
STEP = 5

y0 = int(sys.argv[1]) if len(sys.argv) > 2 else 2002
y1 = int(sys.argv[2]) if len(sys.argv) > 2 else 2025
os.makedirs(CACHE, exist_ok=True)

cells = (pd.read_excel(XLSX, sheet_name="Water", usecols=["cell_centre_lat", "cell_centre_lon"])
         .drop_duplicates())
lats = np.sort(cells.cell_centre_lat.unique()); lons = np.sort(cells.cell_centre_lon.unique())


def day_to_cell(ds):
    """average the 0.2-degree u/v onto 2-degree cells -> dict {(lat,lon): (u,v)}"""
    # one remote read of the whole region, then slice locally
    s = ds[["UCUR", "VCUR"]].isel(TIME=0).sel(
        LATITUDE=slice(lats.min() - 1, lats.max() + 1),
        LONGITUDE=slice(lons.min() - 1, lons.max() + 1)).load()
    rows = []
    for la in lats:
        for lo in lons:
            b = s.sel(LATITUDE=slice(la - 1, la + 1), LONGITUDE=slice(lo - 1, lo + 1))
            rows.append((la, lo, float(b.UCUR.mean(skipna=True)), float(b.VCUR.mean(skipna=True))))
    return rows


def fetch(path):
    for attempt in range(3):
        try:
            with xr.open_dataset(DAP + path) as ds:
                date = pd.Timestamp(ds.TIME.values[0])
                return [(date.year, date.month, *r) for r in day_to_cell(ds)]
        except Exception as e:  # network hiccup: retry
            err = e
    print("FAILED", path, err, flush=True)
    return []


for year in range(y0, y1 + 1):
    out = f"{CACHE}/currents_{year}.csv"
    if os.path.exists(out):
        continue
    xml = requests.get(CAT.format(y=year), timeout=60).text
    paths = sorted(set(re.findall(r'urlPath="([^"]+\.nc)"', xml)))
    paths = [p for i, p in enumerate(paths) if i % STEP == 0]
    with ThreadPoolExecutor(8) as ex:
        res = [r for chunk in ex.map(fetch, paths) for r in chunk]
    pd.DataFrame(res, columns=["year", "month", "cell_centre_lat", "cell_centre_lon", "u", "v"]).to_csv(out, index=False)
    print(year, len(paths), "files", flush=True)

d = pd.concat([pd.read_csv(f"{CACHE}/currents_{y}.csv") for y in range(y0, y1 + 1)])
k = ["year", "month", "cell_centre_lat", "cell_centre_lon"]
g = d.groupby(k)
m = g[["u", "v"]].mean()
m["n_days"] = g.size()
d = d.join(m, on=k, rsuffix="_m")
d["ke"] = 0.5 * ((d.u - d.u_m) ** 2 + (d.v - d.v_m) ** 2)
m["eke"] = d.groupby(k).ke.mean()
m["speed_ms"] = np.hypot(m.u, m.v)
m = m.rename(columns={"u": "u_ms", "v": "v_ms"}).dropna(subset=["u_ms"]).reset_index()
m.to_csv("data/currents_cells.csv", index=False)
print(m.describe().T[["count", "min", "mean", "max"]])
