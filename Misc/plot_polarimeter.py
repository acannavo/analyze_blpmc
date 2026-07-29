"""
BNL 200 MeV Polarimeter — Wire-scanner analysis
------------------------------------------------
Reads Sheet1.csv from the same folder as this script.

Produces three figures:
  Fig 1 · polarimeter_scan.png      — raw scan + per-peak Gaussian fits
  Fig 2 · rl_profiles.png           — RL-reconstructed beam profiles + convergence
  Fig 3 · beam_reconstruction.png   — σ_beam: quadrature subtraction vs RL

Run from VS Code:  Ctrl+F5
Dependencies: numpy, pandas, matplotlib, scipy
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit

# ── Paths ─────────────────────────────────────────────────────────────────────
HERE = Path(__file__).parent
DATA_FILE = HERE / "Sheet1.csv"

# ── 1. Load data ──────────────────────────────────────────────────────────────
df   = pd.read_csv(DATA_FILE, header=None,
                   names=["step", "mm", "VL", "VD"],
                   skiprows=2, dtype=float)
step = df["step"].values
mm   = df["mm"].values
VL   = df["VL"].values
VD   = df["VD"].values
avg  = (VL + VD) / 2.0    # used for all fits

# ── 2. Step ↔ mm conversion ───────────────────────────────────────────────────
slope     = (mm[-1] - mm[0]) / (step[-1] - step[0])
intercept = mm[0] - slope * step[0]

def step_to_mm(s): return slope * np.asarray(s) + intercept
def mm_to_step(m): return (np.asarray(m) - intercept) / slope

# ── 3. Peak definitions ────────────────────────────────────────────────────────
peak_info = {
    270_000: "270 k  =  4 mm foil",
    390_000: "390 k  =  2 mm foil",
    510_000: "510 k  =  Ø 0.3 mm wire",
    630_000: "630 k  =  Ø 0.5 mm wire",
    755_000: "755 k  =  Ø 0.5 mm wire",
}
peak_colours = ["#e63946", "#f4a261", "#457b9d", "#2a9d8f", "#6a4c93"]
peak_windows = {
    270_000: (248_000, 290_000),
    390_000: (372_000, 408_000),
    510_000: (494_000, 530_000),
    630_000: (616_000, 656_000),
    755_000: (736_000, 760_000),    # data ends at 760 k
}
# Effective geometric σ of each target (mm)
#   Rectangular foil of width w  →  σ_target = w / √12
#   Circular wire of diameter D  →  σ_target = D / 4
target_sigma_geom = {
    270_000: 4.0 / np.sqrt(12),    # 1.155 mm
    390_000: 2.0 / np.sqrt(12),    # 0.577 mm
    510_000: 0.3 / 4,              # 0.075 mm
    630_000: 0.5 / 4,              # 0.125 mm
    755_000: 0.5 / 4,              # 0.125 mm
}
# Target shape for kernel construction
target_specs = {
    270_000: ('foil', 4.0),
    390_000: ('foil', 2.0),
    510_000: ('wire', 0.3),
    630_000: ('wire', 0.5),
    755_000: ('wire', 0.5),
}

def gaussian(x, A, mu, sig):
    return A * np.exp(-0.5 * ((x - mu) / sig) ** 2)

# ── 4. Gaussian fits ──────────────────────────────────────────────────────────
fits = {}
print("Gaussian fit results")
print("=" * 72)
for sc, (sl, sh) in peak_windows.items():
    ml, mh = step_to_mm(sl), step_to_mm(sh)
    mask = (mm >= ml) & (mm <= mh)
    xw, yw = mm[mask], avg[mask]
    p0 = [yw.max(), xw[np.argmax(yw)], (mh - ml) / 5]
    popt, pcov = curve_fit(gaussian, xw, yw, p0=p0,
                           bounds=([0, ml, 0.1], [np.inf, mh, mh - ml]))
    perr = np.sqrt(np.diag(pcov))
    xd   = np.linspace(ml, mh, 600)
    fits[sc] = (popt, perr, xd, gaussian(xd, *popt))
    A, mu, sig = popt; dA, dmu, dsig = perr
    print(f"  Peak {sc//1000:3d} k  |  A={A:6.1f}±{dA:4.1f}  "
          f"μ={mu:.3f}±{dmu:.3f} mm  σ={sig:.3f}±{dsig:.3f} mm  "
          f"FWHM={2.355*sig:.3f} mm")
print("=" * 72)

# ── 5. Quadrature subtraction  σ_beam = √(σ_meas² − σ_target²) ───────────────
recon = {}
print("\nQuadrature subtraction  —  beam reconstruction")
print("=" * 72)
for sc, sig_t in target_sigma_geom.items():
    if sc not in fits: continue
    popt, perr, _, _ = fits[sc]
    sig_m, dsig_m = popt[2], perr[2]
    disc = sig_m**2 - sig_t**2
    if disc <= 0:
        print(f"  Peak {sc//1000:3d} k  |  unphysical (σ_meas ≤ σ_target)"); continue
    sig_b  = np.sqrt(disc)
    dsig_b = (sig_m / sig_b) * dsig_m
    recon[sc] = (sig_b, dsig_b, sig_m, sig_t)
    print(f"  Peak {sc//1000:3d} k  |  σ_meas={sig_m:.3f}  σ_target={sig_t:.3f}  "
          f"→  σ_beam={sig_b:.3f} ± {dsig_b:.3f} mm  (err amplif. ×{sig_m/sig_b:.2f})")

sig_vals   = np.array([v[0] for v in recon.values()])
dsig_vals  = np.array([v[1] for v in recon.values()])
weights    = 1.0 / dsig_vals**2
sigma_wmean  = np.sum(weights * sig_vals) / np.sum(weights)
dsigma_wmean = 1.0 / np.sqrt(np.sum(weights))
print(f"\n  Weighted mean :  σ_beam = {sigma_wmean:.3f} ± {dsigma_wmean:.3f} mm")
print(f"  FWHM (beam)   :  {2.355*sigma_wmean:.3f} ± {2.355*dsigma_wmean:.3f} mm")
print("=" * 72)

# ── 6. Richardson-Lucy deconvolution ─────────────────────────────────────────
# Settings
UPSAMPLE = 5      # fine-grid factor relative to native 2 000-step spacing
N_ITER   = 150    # fixed number of RL iterations
EPS      = 1e-10  # numerical floor to avoid log(0) / division by zero

def build_kernel(t_type, param, dx):
    """Normalized PSF kernel on the fine grid."""
    if t_type == 'foil':
        n  = max(3, int(np.ceil(param / dx)) * 2 + 1)   # always odd
        xk = (np.arange(n) - n // 2) * dx
        T  = np.where(np.abs(xk) <= param / 2, 1.0, 0.0)
    else:   # circular wire
        R  = param / 2
        n  = max(3, int(np.ceil(2 * R / dx)) * 2 + 1)
        xk = (np.arange(n) - n // 2) * dx
        T  = np.where(np.abs(xk) < R,
                      2 * np.sqrt(np.maximum(R**2 - xk**2, 0)) / (np.pi * R**2), 0.0)
    s = T.sum()
    return T / s if s > 0 else np.array([1.0])

rl_results = {}
print("\nRichardson-Lucy deconvolution")
print("=" * 72)

for sc, (sl, sh) in peak_windows.items():
    ml, mh = step_to_mm(sl), step_to_mm(sh)
    mask = (mm >= ml) & (mm <= mh)
    xw, Mw = mm[mask], avg[mask]

    # Upsample via cubic spline
    n_fine = (len(xw) - 1) * UPSAMPLE + 1
    xf     = np.linspace(xw[0], xw[-1], n_fine)
    dx_f   = xf[1] - xf[0]
    Mf     = np.maximum(CubicSpline(xw, Mw)(xf), 0.0)

    # Build target kernel
    t_type, param = target_specs[sc]
    T  = build_kernel(t_type, param, dx_f)
    Tf = T[::-1]    # flipped kernel for back-projection (T is symmetric, kept for clarity)

    # RL iterations — B initialized to measured profile
    B      = np.maximum(Mf.copy(), EPS)
    ll_hist = []
    for _ in range(N_ITER):
        fwd    = np.maximum(np.convolve(B, T, mode='same'), EPS)
        B      = np.maximum(B * np.convolve(Mf / fwd, Tf, mode='same'), EPS)
        ll     = float(np.sum(np.where(Mf > 0, Mf * np.log(fwd) - fwd, -fwd)))
        ll_hist.append(ll)

    # Fit Gaussian to RL result for σ_beam,RL
    mu_fit = fits[sc][0][1]
    xc     = xf - mu_fit     # centered coordinates
    try:
        po, pc = curve_fit(gaussian, xc, B,
                           p0=[B.max(), 0.0, 1.0],
                           bounds=([0, -3, 0.05], [np.inf, 3, 8]))
        sig_rl  = abs(po[2])
        dsig_rl = np.sqrt(max(pc[2, 2], 0))
        xd_rl   = np.linspace(xc[0], xc[-1], 500)
        Bfit    = gaussian(xd_rl, *po)
    except Exception as e:
        sig_rl, dsig_rl = np.nan, np.nan
        xd_rl, Bfit = None, None

    rl_results[sc] = (xc, B, np.array(ll_hist), sig_rl, dsig_rl, xd_rl, Bfit)

    sig_q = recon[sc][0] if sc in recon else np.nan
    print(f"  Peak {sc//1000:3d} k  |  "
          f"σ_quad={sig_q:.3f}   σ_RL={sig_rl:.3f} ± {dsig_rl:.3f} mm")

print("=" * 72)

# ══════════════════════════════════════════════════════════════════════════════
# ── Figure 1: raw scan + Gaussian fits ────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
fig1, ax = plt.subplots(figsize=(15, 6.5))
l1, = ax.plot(mm, VL, color="steelblue",  ls="", marker="o", ms=3, alpha=0.55,
              label="VL  (left  detector)")
l2, = ax.plot(mm, VD, color="darkorange", ls="", marker="o", ms=3, alpha=0.55,
              label="VD  (right detector)")
ax.set_xlabel("Position  (mm)", fontsize=13)
ax.set_ylabel("Counts / 2 000-step bin", fontsize=13)
ax.set_title("BNL 200 MeV Polarimeter — Wire-scanner profile  ·  Gaussian fits to (VL+VD)/2",
             fontsize=13, pad=14)
ax.grid(True, linestyle="--", alpha=0.35)
ax.set_ylim(bottom=0, top=avg.max() * 1.55)
ax_s = ax.secondary_xaxis("top", functions=(mm_to_step, step_to_mm))
ax_s.set_xlabel("Position  (step)", fontsize=13)
ax_s.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
glines = []
for (sc, calib), col in zip(peak_info.items(), peak_colours):
    ax.axvline(step_to_mm(sc), color=col, lw=0.8, linestyle=":", alpha=0.5)
    popt, perr, xd, yd = fits[sc]; A, mu, sig = popt; fwhm = 2.355 * sig
    lbl = f"{calib}\n  μ={mu:.2f} mm   σ={sig:.2f} mm   FWHM={fwhm:.2f} mm"
    gl, = ax.plot(xd, yd, color=col, lw=2.4, label=lbl); glines.append(gl)
ax.legend(handles=[l1, l2] + glines, loc="upper left", fontsize=8.5,
          framealpha=0.93, edgecolor="grey", handlelength=1.6,
          borderpad=0.8, labelspacing=0.55)
plt.tight_layout()
plt.savefig(HERE / "polarimeter_scan.png", dpi=150, bbox_inches="tight")

# ══════════════════════════════════════════════════════════════════════════════
# ── Figure 2: RL reconstructed profiles + convergence ─────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
fig2, (ax_prof, ax_ll) = plt.subplots(1, 2, figsize=(15, 5.5))

for (sc, label), col in zip(peak_info.items(), peak_colours):
    xc, B, ll_hist, sig_rl, dsig_rl, xd_rl, Bfit = rl_results[sc]
    Bnorm = B / B.max()
    sig_q = recon[sc][0] if sc in recon else np.nan
    leg   = f"{label}\n  σ_RL={sig_rl:.3f} mm   σ_quad={sig_q:.3f} mm"
    ax_prof.plot(xc, Bnorm, color=col, lw=1.8, alpha=0.85, label=leg)
    if xd_rl is not None:
        ax_prof.plot(xd_rl, Bfit / B.max(), color=col, lw=1.2,
                     linestyle="--", alpha=0.6)
    ax_ll.plot(np.arange(1, N_ITER + 1), ll_hist, color=col, lw=1.6,
               label=label.split("=")[0].strip())

ax_prof.set_xlabel("Position relative to peak centre  (mm)", fontsize=12)
ax_prof.set_ylabel("Normalised beam intensity  B(x) / max", fontsize=12)
ax_prof.set_title("RL-reconstructed beam profiles  (centred)\n"
                  "Dashed = Gaussian fit to RL result", fontsize=11, pad=10)
ax_prof.set_xlim(-6, 6)
ax_prof.grid(True, linestyle="--", alpha=0.35)
ax_prof.legend(fontsize=7.5, framealpha=0.92, loc="upper right", labelspacing=0.6)

ax_ll.set_xlabel("RL iteration", fontsize=12)
ax_ll.set_ylabel("Poisson log-likelihood", fontsize=12)
ax_ll.set_title("RL convergence per target", fontsize=11, pad=10)
ax_ll.grid(True, linestyle="--", alpha=0.35)
ax_ll.legend(fontsize=9, framealpha=0.92, loc="lower right")

plt.tight_layout()
plt.savefig(HERE / "rl_profiles.png", dpi=150, bbox_inches="tight")

# ══════════════════════════════════════════════════════════════════════════════
# ── Figure 3: σ_beam comparison — quadrature subtraction vs RL ───────────────
# ══════════════════════════════════════════════════════════════════════════════
target_labels = {
    270_000: "4 mm foil\n(270 k)",
    390_000: "2 mm foil\n(390 k)",
    510_000: "Ø 0.3 mm wire\n(510 k)",
    630_000: "Ø 0.5 mm wire\n(630 k)",
    755_000: "Ø 0.5 mm wire\n(755 k)",
}
pt_col = {270_000:"#e63946", 390_000:"#f4a261", 510_000:"#457b9d",
          630_000:"#2a9d8f",  755_000:"#6a4c93"}

keys   = list(recon.keys())
xi     = np.arange(len(keys))
dx_off = 0.13

fig3, ax3 = plt.subplots(figsize=(10, 5.2))
for i, sc in enumerate(keys):
    sig_q, dsig_q, _, _ = recon[sc]
    _, _, _, sig_rl, dsig_rl, _, _ = rl_results[sc]
    col = pt_col[sc]
    ax3.errorbar(i - dx_off, sig_q,  yerr=dsig_q,  fmt="o", ms=9,
                 capsize=5, capthick=1.8, lw=1.8, color=col)
    ax3.errorbar(i + dx_off, sig_rl, yerr=dsig_rl, fmt="s", ms=9,
                 capsize=5, capthick=1.8, lw=1.8, color=col, mfc="white")
    ax3.annotate(f"Q:{sig_q:.3f}",  xy=(i - dx_off, sig_q),
                 xytext=(-20, 10), textcoords="offset points",
                 fontsize=7.5, color=col)
    ax3.annotate(f"RL:{sig_rl:.3f}", xy=(i + dx_off, sig_rl),
                 xytext=(4, 10), textcoords="offset points",
                 fontsize=7.5, color=col)

ax3.axhline(sigma_wmean, color="black", lw=1.5, ls="--",
            label=f"Quad. weighted mean = {sigma_wmean:.3f} mm")
ax3.axhspan(sigma_wmean - dsigma_wmean, sigma_wmean + dsigma_wmean,
            alpha=0.10, color="black", label=f"±1σ = {dsigma_wmean:.3f} mm")

legend_handles = [
    Line2D([0],[0], marker='o', color='k', ms=8, ls='none',
           label='Quadrature subtraction  ●'),
    Line2D([0],[0], marker='s', color='k', ms=8, ls='none', mfc='white',
           label='Richardson-Lucy  □'),
    Line2D([0],[0], color='k', lw=1.5, ls='--',
           label=f"Quad. weighted mean = {sigma_wmean:.3f} mm"),
]
ax3.legend(handles=legend_handles, fontsize=9, framealpha=0.92, loc="upper right")
ax3.set_xticks(xi)
ax3.set_xticklabels([target_labels[sc] for sc in keys], fontsize=9)
ax3.set_ylabel("σ_beam  (mm)", fontsize=12)
ax3.set_title("BNL 200 MeV Polarimeter — σ_beam reconstruction\n"
              "Quadrature subtraction  vs  Richardson-Lucy",
              fontsize=11, pad=12)
ax3.set_ylim(0, max(v[0] for v in recon.values()) * 1.55)
ax3.grid(True, axis="y", linestyle="--", alpha=0.4)
ax3_fwhm = ax3.secondary_yaxis("right",
    functions=(lambda s: 2.355*s, lambda f: f/2.355))
ax3_fwhm.set_ylabel("FWHM  (mm)", fontsize=12)
plt.tight_layout()
plt.savefig(HERE / "beam_reconstruction.png", dpi=150, bbox_inches="tight")

plt.show()
print(f"\nFig 1 → polarimeter_scan.png")
print(f"Fig 2 → rl_profiles.png")
print(f"Fig 3 → beam_reconstruction.png")