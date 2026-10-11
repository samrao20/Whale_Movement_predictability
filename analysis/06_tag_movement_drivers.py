"""What drives humpback movement, and does it differ between groups?

Daily positions of the satellite-tagged whales are matched to SST, SST fronts, chlorophyll-a and
IMOS currents. For each group (tagging campaign x zone) we fit, with whale as the grouping unit:
  (1) log daily travel speed        ~ drivers   (linear mixed model, whale random intercept)
  (2) tortuous / area-restricted     ~ drivers   (logistic GEE, whale clusters)
Drivers are z-scored within group, so coefficients are per +1 SD and comparable across drivers.
SST and chlorophyll enter as anomalies relative to the latitude/month norm, because raw SST is
almost a function of latitude. Latitude itself is kept as a control.

Needs: data/tag_points_env.csv, tag_points_currents.csv, tag_points_chl.csv (scripts 06-08).
Writes data/tag_driver_models.csv, analysis/tag_driver_heatmap.png, analysis/tag_driver_binned.png
Run from repo root:  python analysis/06_tag_movement_drivers.py
"""
import warnings
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
warnings.filterwarnings("ignore")

INK, INK2, GRID, SURF, MID = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb", "#f0efec"
BLUE, RED, ORANGE, AQUA = "#2a78d6", "#e34948", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "text.color": INK, "axes.facecolor": SURF, "figure.facecolor": SURF})

# ---- load + merge -----------------------------------------------------------------------------------
p = pd.read_csv("data/tag_points_env.csv", parse_dates=["timestamp"]); p["whale_id"] = p.whale_id.astype(str)
for f in ("currents", "chl"):
    d = pd.read_csv(f"data/tag_points_{f}.csv", parse_dates=["timestamp"]); d["whale_id"] = d.whale_id.astype(str)
    p = p.merge(d.drop_duplicates(["whale_id", "timestamp"]), on=["whale_id", "timestamp"], how="left")
p.loc[p.sst_C < -1.9, ["sst_C", "sst_grad_C_per_100km"]] = np.nan
p["vsouth"] = -p.vcur_ms

ref = pd.read_csv("data-raw/tags/reference.csv"); ref["whale_id"] = ref["animal-id"].astype(str); ref["dep"] = pd.to_datetime(ref["deploy-on-date"])
ref["cohort"] = ref.dep.map(lambda d: "2008 south" if d.year == 2008 else "2009 north" if d.year == 2009
                            else "2010 feed" if d.month <= 3 else "2010 south")
p = p.merge(ref[["whale_id", "cohort", "animal-sex"]].rename(columns={"animal-sex": "sex"}), on="whale_id")
p["lon360"] = p.longitude % 360; p["date"] = p.timestamp.dt.normalize()

# ---- one position per whale-day --------------------------------------------------------------------
cols = ["latitude", "lon360", "sst_C", "sst_grad_C_per_100km", "chl_log10", "vsouth", "gsla_m"]
dly = p.groupby(["whale_id", "cohort", "sex", "date"], as_index=False)[cols].mean().sort_values(["whale_id", "date"])
dly["month"] = dly.date.dt.month

def hav(la1, lo1, la2, lo2):
    q = np.pi / 180
    a = np.sin((la2 - la1) * q / 2) ** 2 + np.cos(la1 * q) * np.cos(la2 * q) * np.sin((lo2 - lo1) * q / 2) ** 2
    return 2 * 6371 * np.arcsin(np.sqrt(a))

g = dly.groupby("whale_id")
dly["dt"] = g.date.diff().dt.days
dly["dist"] = hav(g.latitude.shift(), g.lon360.shift(), dly.latitude, dly.lon360)
dly["speed"] = dly.dist / dly.dt
# straightness over 4 steps (net displacement / path length) using only regular (<=2 d gaps) windows
reg = (dly.dt <= 2).astype(float)
path = dly.dist.where(dly.dt <= 2).groupby(dly.whale_id).transform(lambda s: s.rolling(4, min_periods=4).sum())
net = pd.Series(hav(g.latitude.shift(4), g.lon360.shift(4), dly.latitude, dly.lon360), index=dly.index)
dly["straight"] = (net / path).clip(upper=1)
dly.loc[reg.groupby(dly.whale_id).transform(lambda s: s.rolling(4, min_periods=4).min()) < 1, "straight"] = np.nan

s = dly[(dly.dt >= 1) & (dly.dt <= 2) & dly.speed.notna()].copy()
s = s[s.speed < 250]                               # drop implausible jumps (state-space model artefacts)
s["log_speed"] = np.log(s.speed + 1)
s["ars"] = (s.straight < 0.5).astype(float).where(s.straight.notna())
s["zone"] = pd.cut(s.latitude, [-90, -60, -45, 0], labels=["Antarctic (<60S)", "Southern Ocean (45-60S)", "Subtropical (>45S)"])
s["zone"] = s.zone.cat.rename_categories({"Subtropical (>45S)": "Subtropical/EAC (N of 45S)"})

# ---- anomalies relative to the latitude/month norm ----------------------------------------------------
for v, name in (("sst_C", "sst_anom"), ("chl_log10", "chl_anom")):
    ok = s[v].notna()
    fit = smf.ols(f"{v} ~ latitude + I(latitude**2) + I(latitude**3) + C(month)", data=s[ok]).fit()
    s.loc[ok, name] = fit.resid
    print(f"{v}: {fit.rsquared:.2f} of variance is latitude/season; anomaly kept as driver")
s = s.rename(columns={"sst_grad_C_per_100km": "front"})

DRIVERS = {"sst_anom": "Warmer than usual (SST anomaly)", "chl_anom": "Higher chlorophyll (anomaly)",
           "front": "Stronger SST front", "vsouth": "Southward current", "latitude": "Latitude (control)"}

def z(x): return (x - x.mean()) / x.std()

rows = []
for (coh, zone), d in s.groupby(["cohort", "zone"], observed=True):
    if d.whale_id.nunique() < 3 or len(d) < 60: continue
    use = [k for k in DRIVERS if d[k].notna().mean() >= 0.6 and d[k].std() > 0]
    dd = d.dropna(subset=use + ["log_speed"]).copy()
    if dd.whale_id.nunique() < 3 or len(dd) < 50: continue
    for k in use: dd[k] = z(dd[k])
    form = "log_speed ~ " + " + ".join(use)
    try:
        m = smf.mixedlm(form, dd, groups=dd.whale_id).fit(reml=True)
        for k in use:
            b, se = m.params[k], m.bse[k]
            rows.append(dict(cohort=coh, zone=str(zone), response="log speed", driver=k, coef=b, lo=b - 1.96 * se, hi=b + 1.96 * se,
                             p=m.pvalues[k], n_days=len(dd), n_whales=dd.whale_id.nunique()))
    except Exception as e:
        print("speed model failed", coh, zone, e)
    da = dd.dropna(subset=["ars"])
    if da.ars.nunique() == 2 and da.whale_id.nunique() >= 3:
        try:
            m = smf.gee("ars ~ " + " + ".join(use), "whale_id", da, family=sm.families.Binomial(),
                        cov_struct=sm.cov_struct.Exchangeable()).fit()
            for k in use:
                b, se = m.params[k], m.bse[k]
                rows.append(dict(cohort=coh, zone=str(zone), response="tortuous (log-odds)", driver=k, coef=b, lo=b - 1.96 * se,
                                 hi=b + 1.96 * se, p=m.pvalues[k], n_days=len(da), n_whales=da.whale_id.nunique()))
        except Exception as e:
            print("ars model failed", coh, zone, e)

res = pd.DataFrame(rows); res.to_csv("data/tag_driver_models.csv", index=False)
s.to_csv("data/tag_daily_steps.csv", index=False)
print(res.groupby(["cohort", "zone", "response"]).agg(n_days=("n_days", "first"), n_whales=("n_whales", "first")).to_string())

# ---- heatmaps: standardized effect by group -----------------------------------------------------------------
cmap = LinearSegmentedColormap.from_list("div", [BLUE, MID, RED])
SHORT = {"Subtropical/EAC (N of 45S)": "N of 45S", "Southern Ocean (45-60S)": "45-60S", "Antarctic (<60S)": "S of 60S"}
res["group"] = res.cohort + "\n" + res.zone.map(SHORT)
order = ["2008 south", "2010 south", "2009 north", "2010 feed"]
groups = [f"{c}\n{z}" for c in order for z in SHORT.values()]
groups = [g_ for g_ in groups if g_ in set(res.group)]
fig, axs = plt.subplots(1, 2, figsize=(17, 4.9), gridspec_kw=dict(wspace=0.62))
for ax, resp, title in zip(axs, ["log speed", "tortuous (log-odds)"],
                           ["Travel speed (change in log speed per +1 SD)", "Tortuous path (log-odds per +1 SD; exploratory, unstable)"]):
    r = res[res.response == resp]
    lim = max(0.15, np.nanpercentile(np.abs(r.coef), 95))
    for i, d_ in enumerate(DRIVERS):
        for j, g_ in enumerate(groups):
            c = r[(r.driver == d_) & (r.group == g_)]
            if c.empty:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, color=GRID, alpha=.35, lw=0)); continue
            c = c.iloc[0]; sig = (c.lo > 0) or (c.hi < 0)
            ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, color=cmap((np.clip(c.coef, -lim, lim) / lim + 1) / 2), lw=1, ec=SURF))
            ax.text(j, i, f"{c.coef:+.2f}" + ("*" if sig else ""), ha="center", va="center", fontsize=8.5,
                    color=INK, fontweight="bold" if sig else "normal")
    ax.set_xlim(-.5, len(groups) - .5); ax.set_ylim(len(DRIVERS) - .5, -.5)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels(groups, fontsize=8)
    ax.set_yticks(range(len(DRIVERS))); ax.set_yticklabels(list(DRIVERS.values()), fontsize=9)
    ax.tick_params(length=0); [sp.set_visible(False) for sp in ax.spines.values()]
    ax.set_title(title, loc="left", fontsize=10.5)
fig.suptitle("What is associated with whale movement, by group  (* = 95% interval excludes 0; grey = not enough data)",
             x=0.01, y=1.02, ha="left", fontsize=11.5)
fig.savefig("analysis/tag_driver_heatmap.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ---- binned look at the raw relationships (speed vs chlorophyll anomaly / SST anomaly) ----------------------
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
colz = {"2008 south": BLUE, "2010 south": ORANGE, "2009 north": AQUA}
for ax, v, lab in zip(axs, ["chl_anom", "sst_anom", "front"],
                      ["Chlorophyll anomaly (log10, vs latitude/month norm)", "SST anomaly (°C vs latitude/month norm)", "SST front (°C per 100 km)"]):
    for coh, c in colz.items():
        d = s[(s.cohort == coh) & s[v].notna()]
        if d.whale_id.nunique() < 3 or len(d) < 80: continue
        d = d.assign(bin=pd.qcut(d[v], 6, duplicates="drop"))
        b = d.groupby("bin", observed=True).agg(x=(v, "median"), y=("speed", "median"), n=("speed", "size"), w=("whale_id", "nunique"))
        ax.plot(b.x, b.y, color=c, lw=2, marker="o", ms=4, label=f"{coh} (n={d.whale_id.nunique()} whales)")
    ax.set_xlabel(lab); ax.grid(color=GRID, lw=.6); ax.spines[["top", "right"]].set_visible(False)
axs[0].set_ylabel("Median daily travel speed (km/day)"); axs[0].legend(frameon=False, fontsize=8.5, loc="upper right")
fig.suptitle("Travel speed against each driver (quantile bins, all latitudes)", x=0.01, ha="left", fontsize=11.5)
fig.savefig("analysis/tag_driver_binned.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("done")
