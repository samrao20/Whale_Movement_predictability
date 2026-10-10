import urllib.request, urllib.parse, re, time, io, os, sys
import xarray as xr, numpy as np, pandas as pd

BASE = "https://thredds.aodn.org.au"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data-raw", "sst_southern_2003_2022.csv")

def get(u):
    req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
    return urllib.request.urlopen(req, timeout=180).read()

# 1. catalogue
yearcat = {}
for y in range(2003, 2023):
    xml = get(f"{BASE}/thredds/catalog/IMOS/SRS/SST/ghrsst/L3S-1m/day/{y}/catalog.xml").decode("utf-8", "ignore")
    yearcat[y] = sorted(re.findall(r'urlPath="(IMOS/SRS/SST/ghrsst/L3S-1m/day/\d{4}/[^"]+\.nc)"', xml))

tasks = [(y, f) for y, fs in yearcat.items() for f in fs]
def ym_of(fp):
    fn = fp.split("/")[-1]
    return int(fn[:4]), int(fn[4:6])

rows, done = [], set()
if os.path.exists(OUT):
    prev = pd.read_csv(OUT)
    done = set(zip(prev.year, prev.month))
    rows = prev.to_dict("records")
print(len(tasks), "tasks,", len(done), "already done", flush=True)

def fetch_month(fp, retries=3):
    url = f"{BASE}/thredds/ncss/grid/{fp}?" + urllib.parse.urlencode(
        dict(var="sea_surface_temperature", north=-38, west=143, east=170, south=-60,
             disableProjSubset="on", horizStride=10, addLatLon="true", accept="netcdf4"))
    for a in range(retries):
        try:
            return xr.open_dataset(io.BytesIO(get(url)))
        except Exception:
            if a == retries - 1: raise
            time.sleep(2 ** a)

t0 = time.time()
for k, (y, fp) in enumerate(tasks):
    yy, mm = ym_of(fp)
    if (yy, mm) in done:
        continue
    try:
        ds = fetch_month(fp)
        da = ds.sea_surface_temperature.isel(time=0)
        if float(da.max()) > 100:
            da = da - 273.15
        df = da.to_dataframe(name="sst").reset_index().dropna(subset=["sst"])
        df["cell_centre_lat"] = np.floor((df.lat + 1) / 2) * 2 - 1
        df["cell_centre_lon"] = np.floor((df.lon + 1) / 2) * 2 - 1
        g = df.groupby(["cell_centre_lat", "cell_centre_lon"], as_index=False).sst.mean()
        g["year"] = yy; g["month"] = mm
        rows.extend(g.to_dict("records"))
    except Exception as e:
        print("FAIL", yy, mm, e, flush=True)
    if (k + 1) % 20 == 0:
        pd.DataFrame(rows).to_csv(OUT, index=False)
        print(f"{k+1}/{len(tasks)} ({time.time()-t0:.0f}s)", flush=True)

pd.DataFrame(rows).to_csv(OUT, index=False)
print("DONE", len(rows), "rows", f"{time.time()-t0:.0f}s", flush=True)
