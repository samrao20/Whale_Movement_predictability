"""Bathymetry (NOAA ERDDAP ETOPO 1-arcmin, subsampled to 0.1 deg) summarised onto
the 2-degree grid cells used in the Water sheet -> data/bathymetry_cells.csv

Columns: depth_mean_m, depth_sd_m (roughness/slope proxy), frac_shelf (<200 m),
frac_land, dist_shelf_edge_km (cell centre to nearest 200 m contour point).

Run from the repo root:  python data-raw/02_download_bathymetry.py
"""
import io
import numpy as np
import pandas as pd
import requests

XLSX = "data-raw/Humpback_East_Australia_data.xlsx"
ERDDAP = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.csv"

cells = (pd.read_excel(XLSX, sheet_name="Water", usecols=["cell_centre_lat", "cell_centre_lon"])
         .drop_duplicates().reset_index(drop=True))
lat0, lat1 = cells.cell_centre_lat.min() - 1, cells.cell_centre_lat.max() + 1
lon0, lon1 = cells.cell_centre_lon.min() - 1, cells.cell_centre_lon.max() + 1

q = f"?altitude%5B({lat0}):6:({lat1})%5D%5B({lon0}):6:({lon1})%5D"
r = requests.get(ERDDAP + q, timeout=300)
r.raise_for_status()
b = pd.read_csv(io.StringIO(r.text), skiprows=[1])
b.columns = ["lat", "lon", "alt"]
print("bathymetry points:", len(b))

# shelf-edge points: ocean points within 200 m of the 200 m isobath
edge = b[(b.alt < 0) & (b.alt.between(-250, -150))][["lat", "lon"]].to_numpy()

def hav(lat1, lon1, lat2, lon2):
    p = np.pi / 180
    a = (np.sin((lat2 - lat1) * p / 2) ** 2
         + np.cos(lat1 * p) * np.cos(lat2 * p) * np.sin((lon2 - lon1) * p / 2) ** 2)
    return 2 * 6371 * np.arcsin(np.sqrt(a))

out = []
for _, c in cells.iterrows():
    s = b[(b.lat.between(c.cell_centre_lat - 1, c.cell_centre_lat + 1)) &
          (b.lon.between(c.cell_centre_lon - 1, c.cell_centre_lon + 1))]
    ocean = s[s.alt < 0]
    d = hav(c.cell_centre_lat, c.cell_centre_lon, edge[:, 0], edge[:, 1]).min()
    out.append(dict(
        cell_centre_lat=c.cell_centre_lat, cell_centre_lon=c.cell_centre_lon,
        depth_mean_m=ocean.alt.mean() if len(ocean) else np.nan,
        depth_sd_m=ocean.alt.std() if len(ocean) > 1 else np.nan,
        frac_shelf=(ocean.alt > -200).sum() / max(len(s), 1),
        frac_land=(s.alt >= 0).sum() / max(len(s), 1),
        dist_shelf_edge_km=d))

res = pd.DataFrame(out)
res.to_csv("data/bathymetry_cells.csv", index=False)
print(res.describe().T[["count", "min", "mean", "max"]])
