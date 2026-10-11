"""IMOS OceanCurrent GSLA (delayed-mode) matched to the satellite-tag points
-> data/tag_points_currents.csv  (columns: whale_id, timestamp, ucur_ms, vcur_ms, gsla_m)

GSLA is gridded sea-level anomaly and geostrophic currents (0.2 deg, daily, 1993-2025,
lat -60 to 10 only, so tag points south of -60 get no value). To keep it quick, one
file is sampled every 7 days and each tag point takes the nearest sampled day
(at most 3 days away) and nearest grid cell.

Run from repo root after 06_tag_imos_sst.py:  python data-raw/07_tag_imos_currents.py
"""
import re
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd, requests, xarray as xr

CAT = "https://thredds.aodn.org.au/thredds/catalog/IMOS/OceanCurrent/GSLA/DM/{y}/catalog.xml"
DAP = "https://thredds.aodn.org.au/thredds/dodsC/"

p = pd.read_csv("data/tag_points_env.csv", parse_dates=["timestamp"])
p["lon360"] = p.longitude % 360
p["day"] = p.timestamp.dt.normalize()

# sampled days: every 7th day over the tag period, keep those with points within 3 days
start, end = p.day.min(), p.day.max()
cand = pd.date_range(start, end, freq="7D")
days = np.array(cand.values, dtype="datetime64[D]")
pd_days = p.day.values.astype("datetime64[D]")
nearest = days[np.abs(pd_days[:, None] - days[None, :]).argmin(axis=1)]
p["sample_day"] = nearest
p = p[np.abs((pd_days - nearest).astype(int)) <= 3].copy()

files = {}
for y in sorted({int(str(d)[:4]) for d in p.sample_day.unique()}):
    xml = requests.get(CAT.format(y=y), timeout=60).text
    for u in re.findall(r'urlPath="([^"]+\.nc)"', xml):
        m = re.search(r"_(\d{8})T", u)
        files[np.datetime64(f"{m.group(1)[:4]}-{m.group(1)[4:6]}-{m.group(1)[6:]}")] = u

def work(day):
    sub = p[p.sample_day == day]
    sub = sub[sub.latitude >= -60]
    if sub.empty or day not in files:
        return []
    for _ in range(3):
        try:
            with xr.open_dataset(DAP + files[day]) as ds:
                s = ds[["UCUR", "VCUR", "GSLA"]].isel(TIME=0).sel(
                    LATITUDE=slice(-60, -10), LONGITUDE=slice(95, 185)).load()
            la = xr.DataArray(sub.latitude.values, dims="i"); lo = xr.DataArray(sub.lon360.values, dims="i")
            g = s.sel(LATITUDE=la, LONGITUDE=lo, method="nearest", tolerance=0.3)
            return list(zip(sub.index, g.UCUR.values, g.VCUR.values, g.GSLA.values))
        except Exception as e:
            err = e
    print("FAILED", day, err, flush=True)
    return []

udays = sorted(p.sample_day.unique())
print(len(udays), "sampled days", flush=True)
rows = []
with ThreadPoolExecutor(6) as ex:
    for i, r in enumerate(ex.map(work, udays), 1):
        rows += r
        if i % 10 == 0:
            print(i, "/", len(udays), flush=True)

res = pd.DataFrame(rows, columns=["idx", "ucur_ms", "vcur_ms", "gsla_m"]).set_index("idx")
out = p.join(res)[["whale_id", "timestamp", "ucur_ms", "vcur_ms", "gsla_m"]]
out.to_csv("data/tag_points_currents.csv", index=False)
print(len(out), "points;", out.ucur_ms.notna().sum(), "with currents")
