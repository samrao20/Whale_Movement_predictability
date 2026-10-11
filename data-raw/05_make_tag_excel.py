"""Satellite tag tracks (Andrews-Goff et al. 2023, Movebank, CC0) -> data/humpback_satellite_tags.xlsx

Input (downloaded from the Movebank Data Repository, doi 10.5441/001/1.294):
  data-raw/tags/tracks.csv, data-raw/tags/reference.csv
Run from repo root: python data-raw/05_make_tag_excel.py
"""
import pandas as pd

t = pd.read_csv("data-raw/tags/tracks.csv", low_memory=False)
r = pd.read_csv("data-raw/tags/reference.csv")

t["whale_id"] = t["individual-local-identifier"].astype(str)
t["timestamp"] = pd.to_datetime(t["timestamp"])
t["is_modelled"] = t["modelled"].eq(True)
t["outlier"] = t["algorithm-marked-outlier"].eq(True)
ref = r.assign(whale_id=r["animal-id"].astype(str))[["whale_id", "animal-sex", "animal-life-stage", "deploy-on-date"]]
ref.columns = ["whale_id", "sex", "life_stage", "deploy_on_date"]
t = t.merge(ref[["whale_id", "sex", "life_stage"]], on="whale_id", how="left")
t["date"] = t.timestamp.dt.date
t["year"] = t.timestamp.dt.year
t["month"] = t.timestamp.dt.month

cols = ["whale_id", "sex", "life_stage", "timestamp", "year", "month",
        "location-lat", "location-long"]
cols_out = ["whale_id", "sex", "life_stage", "timestamp", "year", "month", "latitude", "longitude"]

model = t[t.is_modelled].sort_values(["whale_id", "timestamp"])[cols]
model.columns = cols_out
raw = t[~t.is_modelled].sort_values(["whale_id", "timestamp"])[cols + ["argos:lc", "outlier"]]
raw.columns = cols_out + ["argos_location_class", "outlier"]

g = t[t.is_modelled].groupby("whale_id")
summ = g.agg(n_modelled_points=("timestamp", "size"), first_fix=("timestamp", "min"),
             last_fix=("timestamp", "max"), lat_min=("location-lat", "min"),
             lat_max=("location-lat", "max")).reset_index()
summ["days_tracked"] = (summ.last_fix - summ.first_fix).dt.days
summ = ref.merge(summ, on="whale_id", how="left")
summ = summ.merge(r.assign(whale_id=r["animal-id"].astype(str))[["whale_id", "deploy-on-latitude", "deploy-on-longitude", "tag-model"]]
                  .rename(columns={"deploy-on-latitude": "deploy_lat", "deploy-on-longitude": "deploy_lon", "tag-model": "tag_model"}),
                  on="whale_id")

readme = pd.DataFrame({"Satellite tag tracks: east Australian humpback whales": [
    "Source: Andrews-Goff V, Gales N, Childerhouse S, Laverick S, Polanowski A, Double M. 2023. Movebank Data Repository, doi:10.5441/001/1.294 (CC0).",
    "Paper: Andrews-Goff et al. 2023, Biodiversity Data Journal 11:e114729, https://doi.org/10.3897/BDJ.11.e114729",
    "50 whales of breeding stock E1, implantable satellite tags (Wildlife Computers), deployed 2008 to 2010; locations from Oct 2008 to Jul 2011.",
    "Sheet 'Tracks_modelled': state-space-model location estimates (regularised, smoothed positions). Use these for movement analysis.",
    "Sheet 'Tracks_raw_argos': raw Argos fixes. 'outlier' = TRUE rows were flagged by Movebank's filter; argos_location_class is Argos quality (3 best ... B worst).",
    "Sheet 'Whales': one row per whale with sex, life stage, deployment site and a summary of the modelled track.",
    "Sex: m male, f female, u unknown. Coordinates are decimal degrees (WGS84), timestamps are UTC.",
    "Please contact the data authors and cite the dataset if re-using."]})

out = "data/humpback_satellite_tags.xlsx"
with pd.ExcelWriter(out, engine="openpyxl", datetime_format="yyyy-mm-dd hh:mm") as w:
    readme.to_excel(w, sheet_name="Read_me", index=False)
    summ.to_excel(w, sheet_name="Whales", index=False)
    model.to_excel(w, sheet_name="Tracks_modelled", index=False)
    raw.to_excel(w, sheet_name="Tracks_raw_argos", index=False)
    for ws in w.book.worksheets:
        ws.freeze_panes = "A2"
        for i, c in enumerate(ws[1], 1):
            ws.column_dimensions[c.column_letter].width = 110 if ws.title == "Read_me" else max(12, min(24, len(str(c.value)) + 3))
print(len(summ), len(model), len(raw), summ.sex.value_counts().to_dict())
print(summ[["days_tracked", "lat_min", "lat_max"]].describe().round(1).to_string())
