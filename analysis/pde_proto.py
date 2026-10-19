"""Prototype 1-D along-coast whale PDE (design validation before R port).

dC/dt = d/dx(D(x,t) dC/dx) - d/dx(v(x,t) C)

x = latitude, grid = 2-degree water cell centres (odd integers, -29..-55).
v = chemotactic drift on suitability S(T): v = v_max * tanh(alpha * dS/dx)
D = ARS diffusion: D0 / (1 + P/h), P = Chl-based prey index.
"""

import numpy as np
import pandas as pd
import pyreadr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"C:\Users\valpa\Whale_Movement_predictability\data"

# ---------------- grid & forcing ----------------
grid = np.arange(-29.0, -55.0 - 2.0, -2.0)   # decreasing: index 0 = QLD end
N = len(grid)
dx_km = 2.0 * 111.0

water = list(pyreadr.read_r(BASE + r"\water_extended.rds").values())[0]
strip = water[(water.cell_centre_lon >= 145) & (water.cell_centre_lon <= 156)]
clim = (strip.groupby(["cell_centre_lat", "month"], as_index=False)
        .agg(sst=("sea_surface_temp_C", "mean"),
             chl=("chlorophyll_mg_per_m3", "mean")))
chl_fill = clim.chl.median()  # southern cells lack chl; use median as background
sst = np.zeros((N, 12)); chl = np.zeros((N, 12))
for i, lat in enumerate(grid):
    sub = clim[clim.cell_centre_lat == lat].set_index("month")
    for m in range(1, 13):
        sst[i, m - 1] = sub.sst.get(m, np.nan)
        v = sub.chl.get(m, np.nan)
        chl[i, m - 1] = chl_fill if (v is None or np.isnan(v)) else v

# ---------------- suitability ----------------
def S_breeding(T):
    return np.exp(-(T - 21.32) ** 2 / (2 * 4.64 ** 2))

def S_feeding(T):  # saturating cold-plateau fit to GAM envelope
    return 1.0 / (1.0 + np.exp((T - 12.0) / 1.5))

S_b = S_breeding(sst)          # (N, 12)
S_f = S_feeding(sst)

# prey index: Holling-II term uses P = chl (already monthly per cell)
P = chl

# ---------------- parameters ----------------
v_max = 2.5          # km/h  (~60 km/day migration speed)
alpha = 3000.0       # tanh gain on dS/dx (per km; dS/dx ~ 1e-3 on a 2-degree grid)
D0 = 1000.0          # km^2/h  nominal diffusivity (ARS-reduced by prey)
h = 0.4              # Holling half-saturation on chl (mg/m3)
# NOTE: velocity field must be sign-consistent along the corridor in each season
# (uniform seasonal drift direction, magnitude modulated by |dS/dx|). A naive
# v ~ tanh(grad S) field tears the wave apart where dS/dx changes sign.
hours_per_month = 365.0 / 12.0 * 24.0
dt = 0.25            # h  (CFL: v*dt/dx ~ 1.5*0.25/222 << 1)
nsteps_month = int(round(hours_per_month / dt))

# ---------------- flux operators (upwind advection + centred diffusion) ----------------
def dSdx(S):
    """dS/dx per km, with x pointing poleward (grid is decreasing lat)."""
    g = np.zeros_like(S)
    # poleward direction = increasing index; physical dx is negative in lat
    # but we only need the signed gradient along the corridor: S[next]-S[prev]
    g[1:-1] = (S[2:] - S[:-2]) / (2 * dx_km)
    g[0] = (S[1] - S[0]) / dx_km
    g[-1] = (S[-1] - S[-2]) / dx_km
    return g

def step(C, v, D):
    """One explicit step. v positive = poleward (toward higher index)."""
    # faces between cells: face i+0.5
    vf = 0.5 * (v[:-1] + v[1:])
    up = np.where(vf > 0, C[:-1], C[1:])          # upwind value at face
    adv_flux = vf * up                             # poleward flux per face
    Df = 0.5 * (D[:-1] + D[1:])
    dif_flux = -Df * (C[1:] - C[:-1]) / dx_km      # diffusive flux (poleward)
    net = adv_flux + dif_flux
    dC = np.zeros_like(C)
    dC[1:-1] = -(net[1:] - net[:-1]) / dx_km
    # zero-flux boundaries
    return C + dt * dC

# ---------------- run 3 years (monthly forcing, linear interp) ----------------
def month_weights(t_months):
    """Return (m0, m1, w) with t_months fractional month index starting Jan=0."""
    m = t_months % 12.0
    m0 = int(np.floor(m)) % 12
    m1 = (m0 + 1) % 12
    return m0, m1, m - np.floor(m)

C = np.zeros(N)
C[0:4] = 1.0                       # initial pulse off SE QLD / N NSW
snapshots = []                     # (time_months, C)
t_months = 0.0
total_months = 36
steps_done = 0
for step_i in range(int(total_months * hours_per_month / dt)):
    m0, m1, w = month_weights(t_months)
    # migration calendar from fitted seasonal envelopes:
    #   northbound (equatorward, toward QLD)  : May-Aug  -> use S_breeding
    #   southbound (poleward, toward Antarctic): Oct-Feb -> use S_feeding
    #   transition months use whichever envelope is active; v sign flips with dS/dx.
    mo = (t_months % 12.0) + 1.0  # 1..12
    if 5 <= mo <= 8:      # austral autumn/winter: breed in the warm north
        S = (1 - w) * S_b[:, m0] + w * S_b[:, m1]
        g = dSdx(S)
        v = v_max * np.tanh(alpha * np.abs(g)) * -1.0  # drive equatorward
    elif mo >= 10 or mo <= 2:  # austral summer: feed in the cold south
        S = (1 - w) * S_f[:, m0] + w * S_f[:, m1]
        g = dSdx(S)
        v = v_max * np.tanh(alpha * np.abs(g))         # drive poleward
    else:                    # shoulder months: no directed drift, diffuse only
        S = np.ones(N)
        v = np.zeros(N)
    Pp = (1 - w) * P[:, m0] + w * P[:, m1]
    D = D0 / (1.0 + Pp / h)

    C = step(C, v, D)
    if step_i % (nsteps_month // 2) == 0:
        snapshots.append((t_months, C.copy(), S.copy(), v.copy()))
    t_months += dt / hours_per_month

# ---------------- plots ----------------
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
ax = axes[0, 0]
for (tm, Cs, Ss, vs) in snapshots:
    if tm < 12 and abs(tm - round(tm)) < 0.01 and int(round(tm)) % 1 == 0:
        ax.plot(grid, Cs, label=f"m{int(round(tm))}")
ax.set_title("Year 1: density vs latitude (monthly)"); ax.legend(fontsize=6); ax.set_xlabel("lat")

ax = axes[0, 1]
T = np.array([tm for tm, *_ in snapshots]); Z = np.array([Cs for _, Cs, *_ in snapshots]).T
pc = ax.pcolormesh(T, grid, Z, shading="auto", cmap="viridis")
ax.set_title("Hovmoller: density (3 yr)"); ax.set_xlabel("months"); plt.colorbar(pc, ax=ax)

ax = axes[1, 0]
Zs = np.array([Ss for _, _, Ss, _ in snapshots]).T
pc = ax.pcolormesh(T, grid, Zs, shading="auto", cmap="RdYlGn")
ax.set_title("Suitability S(x,t)"); ax.set_xlabel("months"); plt.colorbar(pc, ax=ax)

ax = axes[1, 1]
Zv = np.array([vs for _, _, _, vs in snapshots]).T
pc = ax.pcolormesh(T, grid, Zv, shading="auto", cmap="RdBu",
                   vmin=-v_max, vmax=v_max)
ax.set_title("Advection v(x,t) km/h (+ = poleward)"); ax.set_xlabel("months"); plt.colorbar(pc, ax=ax)

plt.tight_layout()
plt.savefig(r"C:\Users\valpa\Whale_Movement_predictability\analysis\pde_proto_check.png", dpi=110)
print("final total density:", C.sum())
print("min density:", min(s[1].min() for s in snapshots))
# monthly centroid per year
for tm, Cs, _, _ in snapshots:
    if abs(tm - round(tm)) < 0.01 and round(tm) % 3 == 0 and Cs.sum() > 0:
        print(f"month {int(tm):3d} centroid lat {np.sum(grid*Cs)/Cs.sum():6.1f}")
