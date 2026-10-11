"""How the tagged east Australian humpbacks' migration differs between tagging years
(2008 to 2011), with IMOS SST / fronts / currents along the track.

Needs: data/tag_points_env.csv (06_tag_imos_sst.py); optional data/tag_points_currents.csv
(07_tag_imos_currents.py). Writes data/tag_migration_metrics.csv and analysis/tag_*.png.
Run from repo root:  python analysis/05_tag_migration_by_year.py
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"   # validated categorical slots 1-3
plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                     "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6})

p = pd.read_csv("data/tag_points_env.csv", parse_dates=["timestamp"])
p["whale_id"] = p.whale_id.astype(str)
p.loc[p.sst_C < -1.9, ["sst_C", "sst_grad_C_per_100km"]] = np.nan      # ice / retrieval artefacts
if os.path.exists("data/tag_points_currents.csv"):
    c = pd.read_csv("data/tag_points_currents.csv", parse_dates=["timestamp"]); c["whale_id"] = c.whale_id.astype(str)
    p = p.merge(c, on=["whale_id", "timestamp"], how="left")
else:
    for k in ["ucur_ms", "vcur_ms", "gsla_m"]: p[k] = np.nan

ref = pd.read_csv("data-raw/tags/reference.csv"); ref["whale_id"] = ref["animal-id"].astype(str)
ref["dep"] = pd.to_datetime(ref["deploy-on-date"])
def cohort(d):
    if d.year == 2008: return "2008 south"      # off Eden, Oct-Nov 2008, southbound
    if d.year == 2009: return "2009 north"      # off Queensland, Jun 2009, northbound
    return "2010 feed" if d.month <= 3 else "2010 south"   # Antarctic feeding Feb-Mar; Qld Oct 2010
ref["cohort"] = ref.dep.map(cohort)
p = p.merge(ref[["whale_id", "cohort", "dep"]], on="whale_id").sort_values(["whale_id", "timestamp"])
p["lon360"] = p.longitude % 360
p["days"] = (p.timestamp - p.dep).dt.total_seconds() / 86400
p["doy_oct1"] = (p.timestamp - pd.to_datetime(p.timestamp.dt.year.where(p.timestamp.dt.month >= 6, p.timestamp.dt.year - 1).astype(str) + "-10-01")).dt.days

def hav(la1, lo1, la2, lo2):
    q = np.pi / 180
    a = np.sin((la2 - la1) * q / 2) ** 2 + np.cos(la1 * q) * np.cos(la2 * q) * np.sin((lo2 - lo1) * q / 2) ** 2
    return 2 * 6371 * np.arcsin(np.sqrt(a))

# ---- per-whale southbound metrics -------------------------------------------------
rows = []
for w, g in p[p.cohort.isin(["2008 south", "2010 south"])].groupby("whale_id"):
    g = g.sort_values("timestamp").reset_index(drop=True)
    r = dict(whale_id=w, cohort=g.cohort[0], deploy=g.dep[0].date(), n_points=len(g),
             min_lat=g.latitude.min(), days_tracked=(g.timestamp.iloc[-1] - g.timestamp.iloc[0]).days)
    for L in (-40, -50, -60):
        hit = g[(g.latitude <= L)]
        if len(hit) and g.latitude.iloc[0] > L:
            h = hit.iloc[0]
            r[f"date_cross_{abs(L)}S"] = h.timestamp.date(); r[f"days_to_{abs(L)}S"] = h.days
            r[f"lon_at_{abs(L)}S"] = h.lon360
            r[f"sst_at_{abs(L)}S"] = h.sst_C
    if "days_to_40S" in r and "days_to_60S" in r:
        a = g[(g.days >= r["days_to_40S"]) & (g.days <= r["days_to_60S"])]
        if len(a) > 1:
            km = hav(a.latitude.values[:-1], a.lon360.values[:-1], a.latitude.values[1:], a.lon360.values[1:]).sum()
            r["km_40_to_60S"] = km
            r["days_40_to_60S"] = r["days_to_60S"] - r["days_to_40S"]
            r["km_per_day_40_to_60S"] = km / max(r["days_40_to_60S"], 0.5)
            r["mean_sst_leg"] = a.sst_C.mean(); r["mean_front_leg"] = a.sst_grad_C_per_100km.mean()
    rows.append(r)
m = pd.DataFrame(rows); m.to_csv("data/tag_migration_metrics.csv", index=False)

cols = ["n_points", "min_lat", "days_to_40S", "days_to_50S", "days_to_60S", "lon_at_60S",
        "km_per_day_40_to_60S", "sst_at_40S", "mean_sst_leg", "mean_front_leg"]
print("whales per southbound cohort:", m.cohort.value_counts().to_dict())
print(m.groupby("cohort")[cols].median().round(1).T.to_string())
print("n whales with each metric:\n", m.groupby("cohort")[cols].count().T.to_string())
for k in ("date_cross_40S", "date_cross_60S"):
    d = pd.to_datetime(m[k]); print(k, d.groupby(m.cohort).median().dt.date.to_dict())

# ---- fig 1: tracks, one panel per tagging cohort (colour = days since tagging) ------------
names = {"2008 south": "Oct-Nov 2008\nsouthbound from Eden", "2009 north": "Jun 2009\nnorthbound from Qld",
         "2010 feed": "Feb-Mar 2010\ntagged in Antarctic waters", "2010 south": "Oct 2010\nsouthbound from Qld"}
fig, axs = plt.subplots(1, 4, figsize=(16, 5.6), sharey=True)
for ax, k in zip(axs, ["2008 south", "2009 north", "2010 feed", "2010 south"]):
    g = p[p.cohort == k]
    for _, gg in g.groupby("whale_id"):
        ax.plot(gg.lon360, gg.latitude, color="#9aa4b2", lw=0.5, alpha=0.5, zorder=1)
    s = ax.scatter(g.lon360, g.latitude, c=g.days, s=3, cmap="Blues", vmin=-30, vmax=g.days.max(), zorder=2)
    ax.set_title(f"{names[k]}\n{g.whale_id.nunique()} whales", fontsize=9.5, loc="left", color=INK, pad=8)
    ax.set_xlabel("Longitude (°E)"); ax.set_xlim(100, 190)
axs[0].set_ylabel("Latitude (°)"); axs[0].set_ylim(-72, -12)
fig.colorbar(s, ax=axs, shrink=0.7, pad=0.01, label="Days since tagging")
fig.suptitle("Satellite-tag tracks by tagging campaign (smoothed positions)", x=0.01, y=1.06, ha="left", fontsize=12)
fig.savefig("analysis/tag_tracks_by_cohort.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ---- fig 2: southbound timing, 2008 vs 2010 ---------------------------------------------
fig, ax = plt.subplots(figsize=(8.5, 5))
col = {"2008 south": BLUE, "2010 south": ORANGE}
for k, g in p[p.cohort.isin(col)].groupby("cohort"):
    g = g[g.doy_oct1.between(0, 200)]
    for _, gg in g.groupby("whale_id"):
        ax.plot(gg.doy_oct1, gg.latitude, color=col[k], lw=0.7, alpha=0.35)
    gb = g.assign(bin=(g.doy_oct1 // 5) * 5 + 2.5).groupby("bin")
    b = gb.latitude.median()[gb.whale_id.nunique() >= 4]     # medians only where >=4 whales are tracked
    ax.plot(b.index, b.values, color=col[k], lw=2.2)
    ax.text(b.index[-1] + 2, b.values[-1], f"{k.split()[0]} tagging\n(n={g.whale_id.nunique()})", color=INK, va="center", fontsize=9)
for L in (-40, -60): ax.axhline(L, color=INK2, lw=0.6, ls=":")
ax.set_xlim(0, 215); ax.set_xlabel("Days since 1 October"); ax.set_ylabel("Latitude (°)")
ax.set_title("Southbound migration timing, 2008 vs 2010 tagging\nthin = each whale, thick = median where 4+ whales are tracked", loc="left", fontsize=10.5)
fig.savefig("analysis/tag_southbound_timing_2008_vs_2010.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ---- fig 3: environment met along the way (binned by latitude) --------------------------
env = [("sst_C", "SST (°C)"), ("sst_grad_C_per_100km", "SST front strength (°C per 100 km)")]
if p.ucur_ms.notna().any():
    p["vcur_south_ms"] = -p.vcur_ms; env.append(("vcur_south_ms", "Southward current (m/s), IMOS GSLA"))
fig, axs = plt.subplots(1, len(env), figsize=(5 * len(env), 4.6))
for ax, (v, lab) in zip(np.atleast_1d(axs), env):
    for k, c in col.items():
        g = p[(p.cohort == k) & p[v].notna() & p.latitude.between(-62, -20)]
        gb = g.assign(bin=(g.latitude // 3) * 3 + 1.5).groupby("bin")
        b = gb[v].agg(["median", "size"]); b = b[(b["size"] >= 30) & (gb.whale_id.nunique() >= 3)]
        ax.plot(b["median"], b.index, color=c, lw=2, marker="o", ms=4, label=f"{k.split()[0]} tagging")
    ax.set_xlabel(lab); ax.set_ylabel("Latitude (°)")
np.atleast_1d(axs)[0].legend(frameon=False, loc="lower right")
fig.suptitle("Conditions at the whales' positions by latitude (median, 3° bins)", x=0.01, ha="left", fontsize=11)
fig.savefig("analysis/tag_environment_along_track.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("saved figures")
