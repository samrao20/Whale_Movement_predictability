# Whale Movement Predictability — Model Status, Known Issues & Improvement Plan
**Date:** 2026-10-19 · **Repo:** github.com/samrao20/Whale_Movement_predictability
**Scope:** 1-D and 2-D advection–reaction–diffusion models for east-Australia humpback migration under ocean warming.

---

## 1. What the model does now

| Component | Script | Status |
|---|---|---|
| 1-D PDE (advection–diffusion, mass-conserving) | `05_pde_1d.R` | ✅ validated |
| Tag-derived movement params (v_max, D0, monthly drift) | `06_track_parameters.R` | ✅ |
| IMOS GSLA currents in advection | `07_pde_currents.R` | ✅ (second-order effect) |
| Interannual hindcast 2002–2015 | `08_hindcast_validation.R` | ✅ |
| Extended domain −29°..−65° | `09_pde_extended.R` | ✅ **best 1-D model: corr 0.84** |
| Demographic term R(C,T) | `10_pde_reaction.R` | ✅ corr 0.86 |
| Thermal-cue departure (climate-responsive timing) | `11_thermal_departure.R` | ✅ scenario device |
| Seasonal sightings density map | `12_sightings_map.R` | ✅ |
| 2-D drift field from tag fixes | `13_grid_2d_drift.R` | ✅ (95% fallback) |
| 2-D PDE (land mask, corridor mask, lat-profile validation) | `14_pde_2d.R` | ✅ Jun–Aug spearman 0.56–0.67 |

**Headline capabilities:**
- Northbound migration (May–Aug) is reproduced well in 1-D *and* 2-D (Jun centroid err +1.6°, Jul +0.1°; 2-D lat-profile spearman up to 0.67).
- The model conserves mass exactly, stays non-negative, and is demographically self-consistent (self-sustaining annual cycle).
- Climate response: +2 °C advances migration timing ~2 weeks (shoulder months), via a thermal-cue phase shift (β = 0.75 mo/°C, capped ±0.5 mo).

---

## 2. Known issues — model structure

### 2.1 Southbound leg fails (THE major gap)
- **Symptom:** 1-D centroid err Oct −24°, Nov −20°; 2-D lat-profile Oct spearman **−0.62** (anticorrelated).
- **Diagnosis:** tags show Oct southbound meridional speed hugely spread (−4.6 to −0.06 km/h, 10th–90th pct) — some whales sprint, many linger to feed en route. A single advection speed cannot represent this; the model has no foraging-stop behaviour.
- **Tested and ruled out:** latitude-stratified drift (corr got worse, 0.72); diffusion scaling ×5–×20 (no effect — symmetric diffusion doesn't move the centroid); 1-month drift phase lag (helps Oct but breaks Jun); corridor mask (no effect on Oct).
- **Fix required:** mechanistic foraging-stop behaviour (e.g. Chl-triggered residency: advection speed down-weighted where P = β·Chl·exp(−γ(T−T_upw)) is high, per the framework PDF's Holling-II formulation). This is a *process* addition, not tuning.

### 2.2 2-D drift field is 95% fallback
- Only 84 of 1,680 bin-months have ≥10 tag steps. Outside the corridor, advection is the corridor-wide monthly curve broadcast over longitude.
- Quiver maps (script 13) honestly flag fitted (red) vs fallback (grey) bins.
- **Implication:** 2-D zonal structure is trustworthy only along the migration corridor. Tasman/southern Coral Sea deployments would change this.

### 2.3 Thermal-cue β is not fitted
- Tags span 2008–2011 only (4 seasons) — an interannual timing–SST relationship cannot be estimated. β = 0.75 mo/°C is a literature-style sensitivity parameter; sweep over 0–1.0 reported. Interannual validation shows the cue is negligible at present-day anomaly scales (±0.4 °C) — it is a *scenario device*, not a validated mechanism.
- Cue regions are proxies (grid starts at −29°, forcing at −45°; true breeding-ground −25..−15 and feeding-ground −65..−55 regions off-grid).

### 2.4 Summer pile-up at southern boundary
- Model summer centroid sits at −65° (domain wall); observed −58..−62°. Density has nowhere to go and R(C,T) is smooth, so it accumulates at the edge. Wall-loss mechanism exists in script 10 (`USE_WALL_LOSS`) but is off by default (drains the population to ~2% mass in 3 yr if uncompensated).

### 2.5 Effort bias in sightings (data, not model)
- Presence-only sightings + whale-watch effort: September observed centroid (−28.4°) sits *north* of August (−26.6°) — Hervey Bay boats, not biology. Validation correlations are understated; never treat sighting density as abundance.
- `individual_count` missing ~43% of records → presence-only throughout.

### 2.6 Currents are second-order
- IMOS GSLA u/v in advection shifts monthly centroids <1° (EAC opposes southbound at only ~1–1.5 km/h). Keep for completeness; expect no skill from it.

### 2.7 Chl-a as prey proxy
- Chl proxies phytoplankton, not krill/forage fish — the dominant uncertainty in the ARS/diffusion and any foraging-stop term. Also 23% missing; entirely absent south of −55° (median fallback). Mysticetes barely feed on migration anyway — flag in writeup.

---

## 3. Known issues — numerics & code (fixed, but don't regress)

- **Boundary mass leak (fixed in 05/07/08+):** boundary cells must own their one interior face; freezing dC at boundaries leaked ~2/3 of mass over 3 years.
- **`aggregate(cbind(x,y) ~ g)` na.omit trap:** base R drops rows with NA in *any* response column before `na.rm` applies — one all-NA Chl month wiped everything ("no rows to aggregate" on non-empty data). Aggregate SST and Chl separately.
- **Floating-point month drift:** accumulated `t <- t + dt/hpm` gives 1.00000000001571; use `round()` before month grouping.
- **`ifelse()`/`pmax()` strip matrix dims** (script 14 bugs): rebuild with `matrix(..., nla, nlo)`.
- **Variable-name collision:** never name a loop month index `m0` when a mortality parameter `m0` exists (silently set mortality = 3, population extinct).
- **R `read.csv` reads logicals as character** ("True"), so `tr$visible == TRUE` silently drops all rows (script 13 bug).
- **RStudio session hygiene:** stale editor buffers override disk edits; masked base functions (`aggregate` shadowed by a workspace object) and hand-run functions (`getSrcFilename(fn) == ""`) caused hours of phantom errors. First response to any impossible error: `rm(list=ls())` + Restart R + re-source.

---

## 4. Improvement plan (priority order)

1. **Foraging-stop behaviour (the big one).** Modulate advection speed by prey field: v_eff = v_drift / (1 + P/h_v) with P the Holling-II prey proxy. Mechanistically motivated (matches framework PDF), directly targets the Oct/Nov gap, and it's the term that makes the model answer "what happens when prey distribution shifts with climate."
2. **Write up the 1-D findings** as the methods backbone — the story is complete: validated northbound mechanism, climate-cue timing response, honest southbound gap.
3. **More tag data for the 2-D field** (Tasman, southern Coral Sea, southbound leg). Alternatively ingest AATAMS/IMOS animal-tracking DB beyond the E1 set.
4. **Fit β (thermal cue)** once ≥10 seasons of tags exist, or borrow from other baleen species' published phenology–SST relationships.
5. **Southern boundary:** extend domain to full Southern Ocean or implement a stock-exchange boundary condition (whales leaving eastward/westward), replacing both the pile-up and the crude wall loss.
6. **Effort correction for validation:** model whale-watch effort (distance from ports, season) or restrict validation to systematic surveys; would raise apparent skill without touching the model.
7. **Krill proxy:** replace/augment Chl-a with a krill-relevant product (e.g. IMOS acoustic, or copepod/krill climatologies) for the feeding-ground dynamics.
8. **Numerics:** when scaling beyond the explicit prototype, move to `ReacTran` + `deSolve` (method-of-lines); positivity-preserving and far faster.

---

## 5. Reproduce

```r
# In RStudio, working directory = repo root (Files pane -> More -> Set As Working Directory)
source("analysis/06_track_parameters.R")   # needs data-raw/tags/tracks.csv
source("analysis/05_pde_1d.R")             # 1-D baseline
source("analysis/09_pde_extended.R")       # best 1-D model (corr 0.84)
source("analysis/10_pde_reaction.R")       # + R(C,T) (corr 0.86)
source("analysis/11_thermal_departure.R")  # climate scenario (several minutes)
source("analysis/13_grid_2d_drift.R")      # 2-D drift field (run before 14)
source("analysis/14_pde_2d.R")             # 2-D model
source("analysis/12_sightings_map.R")      # sightings map (needs install.packages("maps"))
```
Key fitted parameters: v_max ≈ 4.9 km/h (p90 directed), D0 ≈ 1 km²/h, h = 0.4 mg/m³, r_max = 0.25/mo, m0 = 0.05/mo (mismatch mortality), β = 0.75 mo/°C (sensitivity param). Fitted thermal envelopes: breeding T_opt 21.3 °C σ 4.6; feeding cold-plateau (≈1 below 12 °C).
