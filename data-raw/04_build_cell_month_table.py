"""Join sightings + water + SOI + bathymetry + currents into one cell-month table
-> data/cell_month_2002_2015.csv

One row per (2-degree cell, year, month) of the Water sheet, 2002-2015 (years where
chlorophyll is available). Whale columns are counts of sightings falling in that
cell-month; zeros mean "no sighting recorded", NOT confirmed absence.

Run from the repo root (after scripts 01-03):  python data-raw/04_build_cell_month_table.py
"""
import numpy as np
import pandas as pd

XLSX = "data-raw/Humpback_East_Australia_data.xlsx"
Y0, Y1 = 2002, 2015
K = ["cell_centre_lat", "cell_centre_lon", "year", "month"]

# ---- water grid ---------------------------------------------------------------
water = pd.read_excel(XLSX, sheet_name="Water").rename(
    columns={"salinity_psu_SODA_to2015": "salinity_psu"})
water = water[water.year.between(Y0, Y1)].copy()

# SST anomaly = SST minus that cell's mean for the same calendar month (Y0-Y1)
clim = water.groupby(["cell_centre_lat", "cell_centre_lon", "month"]).sea_surface_temp_C.transform("mean")
water["sst_anom_C"] = water.sea_surface_temp_C - clim
water["chl_log10"] = np.log10(water.chlorophyll_mg_per_m3.where(water.chlorophyll_mg_per_m3 > 0))

# previous-month values within each cell (prey/water respond with a lag)
water["t"] = water.year * 12 + water.month
water = water.sort_values(K)
prev = water[["cell_centre_lat", "cell_centre_lon", "t", "sea_surface_temp_C", "chl_log10"]].copy()
prev["t"] += 1
prev = prev.rename(columns={"sea_surface_temp_C": "sst_lag1_C", "chl_log10": "chl_log10_lag1"})
water = water.merge(prev, on=["cell_centre_lat", "cell_centre_lon", "t"], how="left").drop(columns="t")

# ---- sightings -> cells -------------------------------------------------------
h = pd.read_excel(XLSX, sheet_name="Humpback")
n_all = len(h)
h = h.dropna(subset=["year", "month"])
h = h[h.year.between(Y0, Y1)].copy()
n_period = len(h)
h["cell_centre_lat"] = 2 * np.floor(h.latitude / 2) + 1   # cell centres are odd degrees
h["cell_centre_lon"] = 2 * np.floor(h.longitude / 2) + 1
h["year"] = h.year.astype(int); h["month"] = h.month.astype(int)

sight = h.groupby(K).agg(
    n_sightings=("record_id", "size"),
    n_individuals=("individual_count", "sum"),          # only records with a count
    n_with_count=("individual_count", "count"),
    n_datasets=("dataset", "nunique")).reset_index()

df = water.merge(sight, on=K, how="left")
lost = n_period - int(df.n_sightings.sum(skipna=True))   # sightings in cells not in Water grid
for c in ["n_sightings", "n_individuals", "n_with_count", "n_datasets"]:
    df[c] = df[c].fillna(0).astype(int)
df["present"] = (df.n_sightings > 0).astype(int)

# observer-effort proxy: all sightings (any cell) that year, and in that month of the year
df = df.merge(h.groupby("year").size().rename("effort_year"), on="year", how="left")
df = df.merge(h.groupby(["year", "month"]).size().rename("effort_year_month"), on=["year", "month"], how="left")
df[["effort_year", "effort_year_month"]] = df[["effort_year", "effort_year_month"]].fillna(0).astype(int)

# ---- other covariates ---------------------------------------------------------
df = df.merge(pd.read_csv("data/soi_monthly.csv"), on=["year", "month"], how="left")
df = df.merge(pd.read_csv("data/bathymetry_cells.csv"), on=["cell_centre_lat", "cell_centre_lon"], how="left")
cur = pd.read_csv("data/currents_cells.csv").rename(columns={"n_days": "cur_n_days"})
df = df.merge(cur, on=K, how="left")

df = df.sort_values(K).reset_index(drop=True)
df.to_csv("data/cell_month_2002_2015.csv", index=False)

print(f"sightings total {n_all}; with year+month in {Y0}-{Y1}: {n_period}; "
      f"falling outside the Water grid cells: {lost}")
print(f"rows {len(df)}, cell-months with a sighting {df.present.sum()} ({df.present.mean():.1%})")
print(df.isna().mean().round(3)[lambda s: s > 0].to_string())
