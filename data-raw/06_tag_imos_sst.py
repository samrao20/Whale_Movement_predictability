"""IMOS SRS GHRSST L3S-1m monthly SST (and SST gradient = front strength) matched to
the satellite-tag points -> data/tag_points_env.csv

For every state-space-model tag location, takes the month's SST at the nearest 0.2-deg
grid cell. sst_grad_C_per_100km is the local SST gradient magnitude, a front indicator.
Tags span ~100E to -176 (dateline), so two bounding boxes are fetched per month.

Run from repo root (needs data-raw/tags/tracks.csv):  python data-raw/06_tag_imos_sst.py
Needs: pandas, numpy, xarray, netCDF4, requests
"""
import os, re, tempfile, time, urllib.parse
import numpy as np, pandas as pd, requests, xarray as xr

BASE = "https://thredds.aodn.org.au"
CATV = BASE + "/thredds/catalog/IMOS/SRS/SST/ghrsst/L3S-1m/day/{y}/catalog.xml"

t = pd.read_csv("data-raw/tags/tracks.csv", low_memory=False)
t = t[t.modelled == True].copy()
t["whale_id"] = t["individual-local-identifier"].astype(str)
t["timestamp"] = pd.to_datetime(t.timestamp)
t["latitude"] = t["location-lat"]; t["longitude"] = t["location-long"]
t["year"] = t.timestamp.dt.year; t["month"] = t.timestamp.dt.month
t = t[["whale_id", "timestamp", "year", "month", "latitude", "longitude"]].reset_index(drop=True)

def fetch(fp, west, east, tries=3):
    q = urllib.parse.urlencode(dict(var="sea_surface_temperature", north=-15, south=-70, west=west,
        east=east, disableProjSubset="on", horizStride=10, addLatLon="true", accept="netcdf4"))
    for a in range(tries):
        try:
            r = requests.get(f"{BASE}/thredds/ncss/grid/{fp}?{q}", timeout=240); r.raise_for_status()
            with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as f:
                f.write(r.content)
            try:
                with xr.open_dataset(f.name) as ds:
                    return ds.load()
            finally:
                os.unlink(f.name)
        except Exception:
            if a == tries - 1: raise
            time.sleep(2 ** a)

def month_field(fp):
    parts = [fetch(fp, 100, 180), fetch(fp, -180, -170)]
    da = [p.sea_surface_temperature.isel(time=0) for p in parts]
    da = [d - 273.15 if float(d.max()) > 100 else d for d in da]
    # the server returns the west=-180 box as 180..190; shift to -180..-170 and drop the overlap
    da[1] = da[1].sel(lon=da[1].lon > 180).assign_coords(lon=da[1].lon[da[1].lon > 180] - 360)
    return xr.concat(da, dim="lon").sortby("lon")

def grad_100km(da):
    la, lo = da.lat.values, da.lon.values
    dy = np.gradient(la) * 111.0                        # km per row
    dx = np.gradient(lo)[None, :] * 111.0 * np.cos(np.deg2rad(la))[:, None]
    gy = np.gradient(da.values, axis=0) / dy[:, None]
    gx = np.gradient(da.values, axis=1) / dx
    return np.hypot(gx, gy) * 100.0

out = t.assign(sst_C=np.nan, sst_grad_C_per_100km=np.nan)
for (y, m), idx in t.groupby(["year", "month"]).groups.items():
    xml = requests.get(CATV.format(y=y), timeout=60).text
    fp = next((p for p in re.findall(r'urlPath="([^"]+\.nc)"', xml) if p.split("/")[-1].startswith(f"{y}{m:02d}")), None)
    if fp is None:
        print("no file", y, m); continue
    da = month_field(fp)
    g = xr.DataArray(grad_100km(da), coords=da.coords, dims=da.dims)
    sub = t.loc[idx]
    la = xr.DataArray(sub.latitude.values, dims="p"); lo = xr.DataArray(sub.longitude.values, dims="p")
    out.loc[idx, "sst_C"] = da.sel(lat=la, lon=lo, method="nearest", tolerance=0.5).values
    out.loc[idx, "sst_grad_C_per_100km"] = g.sel(lat=la, lon=lo, method="nearest", tolerance=0.5).values
    print(y, m, len(idx), "points", f"{out.loc[idx, 'sst_C'].notna().mean():.0%} with SST", flush=True)

out.to_csv("data/tag_points_env.csv", index=False)
print(out[["sst_C", "sst_grad_C_per_100km"]].describe().round(2).to_string())
