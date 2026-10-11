"""ESA CCI monthly ocean-colour chlorophyll-a (multi-sensor, 0.04 deg, NOAA CoastWatch ERDDAP
dataset pmlEsaCCI60OceanColorMonthly) matched to the satellite-tag points
-> data/tag_points_chl.csv  (whale_id, timestamp, chl_mg_m3, chl_log10)

Each point takes the month's value at the nearest ~0.2 deg cell (fields are strided by 5).
Two boxes are fetched per month because the tracks cross the dateline.

Run from repo root after 06_tag_imos_sst.py:  python data-raw/08_tag_chlorophyll.py
"""
import os, tempfile
import numpy as np, pandas as pd, requests, xarray as xr

B = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/pmlEsaCCI60OceanColorMonthly"
p = pd.read_csv("data/tag_points_env.csv", parse_dates=["timestamp"])[["whale_id", "timestamp", "year", "month", "latitude", "longitude"]]
p["whale_id"] = p.whale_id.astype(str)

def box(t, lo0, lo1):
    q = f"{B}.nc?chlor_a[({t}):1:({t})][(-70):5:(-15)][({lo0}):5:({lo1})]"
    for a in range(3):
        try:
            r = requests.get(q, timeout=300); r.raise_for_status(); break
        except Exception:
            if a == 2: raise
    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as f:
        f.write(r.content)
    try:
        with xr.open_dataset(f.name) as d:
            return d.chlor_a.isel(time=0).load()
    finally:
        os.unlink(f.name)

out = p.assign(chl_mg_m3=np.nan)
for (y, m), idx in p.groupby(["year", "month"]).groups.items():
    t = f"{y}-{m:02d}-01T00:00:00Z"
    # one field spanning 100E..180 and -180..-170 (the tracks cross the dateline)
    f = xr.concat([box(t, -180, -170), box(t, 100, 179.99)], dim="longitude")
    sub = p.loc[idx]
    la = xr.DataArray(sub.latitude.values, dims="i"); lo = xr.DataArray(sub.longitude.values, dims="i")
    near = f.sel(latitude=la, longitude=lo, method="nearest")
    # drop matches whose nearest cell is more than 0.5 deg away (points outside the fetched boxes)
    ok = (np.abs(near.latitude.values - sub.latitude.values) <= 0.5) & (np.abs(near.longitude.values - sub.longitude.values) <= 0.5)
    vals = np.where(ok, near.values, np.nan)
    out.loc[idx, "chl_mg_m3"] = vals
    print(y, m, len(idx), f"{np.isfinite(vals).mean():.0%} with chlorophyll", flush=True)

out["chl_log10"] = np.log10(out.chl_mg_m3.where(out.chl_mg_m3 > 0))
out[["whale_id", "timestamp", "chl_mg_m3", "chl_log10"]].to_csv("data/tag_points_chl.csv", index=False)
print(out.chl_mg_m3.describe().round(3).to_string())
