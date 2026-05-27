# blpmc Analysis Framework

ROOT-based analysis framework for proton-carbon elastic scattering polarimetry. Extracts the analyzing power **Aᵧ** from Geant4 Monte Carlo simulations as a function of a scanned geometry parameter (e.g. target transverse position, TRANSX).

---

## Typical Usage

```bash
# 1. Go to the Analysis folder
cd ~/Analysis

# 2. Start ROOT (rootlogon.C compiles the classes automatically)
root

# 3. Load the macro and run the full parameter scan
.L analyze_blpmc.C
analyze_blpmc_scan()

# 4. Draw results
draw_scan_p12()       # Ay vs scan parameter for 12.0 deg detector
draw_scan_p16()       # Ay vs scan parameter for 16.2 deg detector
draw_scan_both()      # both detectors

# 5. Draw raw scintillator hit counts vs scan parameter
draw_count_vs_param() # counts per scintillator (left/right, spin-up/down)

# 6. Inspect a specific (angle, param) point
inspect("16p2", 0.0)        # all spins combined
inspect("12p0", 10.0, 1)    # spin-up only
inspect("16p2", 5.0, -1)    # spin-down only
insp_draw_all()              # all plots for the inspected point
```

---

## Repository Structure

```
analyze_blpmc/
├── analyze_blpmc.C   ← main analysis macro
├── PolAnalysis.h     ← asymmetry and Ay calculations
├── PolHistograms.h   ← histogram filling and drawing
├── rootlogon.C       ← auto-compiles classes on ROOT startup
└── README.md
```

Data lives outside the repository and is linked in via a symlink:

```bash
ln -s /path/to/Data ~/Analysis/Data
```

---

## Scan Parameter

The scan parameter is set by two constants at the top of `analyze_blpmc.C`:

```cpp
const TString kScanParam = "TRANSX";   // folder name under Data/
const TString kScanUnit  = "mm";       // axis label unit
```

Supported values:

| kScanParam | kScanUnit | Description |
|------------|-----------|-------------|
| `TRANSX`   | `mm`      | Target transverse position X |
| `TRANSY`   | `mm`      | Target transverse position Y |
| `TRANSZ`   | `mm`      | Target longitudinal position |
| `ROTX`     | `deg`     | Target rotation around X |
| `ROTY`     | `deg`     | Target rotation around Y |
| `ROTZ`     | `deg`     | Target rotation around Z |

All folder discovery, graph titles, axis labels, canvas names, and inspection paths follow automatically from these two constants — no other changes needed.

---

## Data Structure

Simulations are organized by reaction type, spin state, detector angle, and parameter value:

```
Data/{kScanParam}/
└── pC_{Elas|Inel443}_200MeV_MT_P80_{SpinUp|SpinDown}_{12p0|16p2}_MEYER/
    └── {0p0 .. 20p0}/
        ├── blpmc_YYYY-MM-DD_HH-MM-SS.root   ← simulation output
        ├── config.conf                        ← per-run configuration
        ├── pluto.root -> ...                  ← symlink to PLUTO input
        └── run.log                            ← Geant4 output log
```

Each ROOT file contains a single TTree named `tree` with branches for energy deposition in each scintillator of both detector arms (P12: 3 scintillators, P16: 4 scintillators). Both detector branches are present in every file regardless of which detector angle was simulated:

```
eDepP12LS0, eDepP12LS1, eDepP12LS2   ← P12 Left arm
eDepP12RS0, eDepP12RS1, eDepP12RS2   ← P12 Right arm
eDepP16LS0, ..., eDepP16LS3          ← P16 Left arm
eDepP16RS0, ..., eDepP16RS3          ← P16 Right arm
```

---

## Running Simulations — scan.sh

`scan.sh` loops over all 8 PLUTO files (2 reactions × 2 spin states × 2 angles) and all parameter values, launching Geant4 jobs with a configurable parallelism limit.

Key parameters at the top of the script:

```bash
START=0.0      # first parameter value
END=20.0       # last parameter value
STEP=1.0       # step size
MAX_JOBS=4     # max simultaneous Geant4 processes
```

The script automatically skips folders where a valid ROOT file (>10 KB) already exists, making it safe to re-run after partial failures.

```bash
cd ~/blpmc/build
./scripts/scan.sh
```

> **Note:** Exit code 139 (segfault on exit) is expected and harmless — it occurs during Geant4/PLUTO destructor cleanup after the ROOT file has already been written and closed. Data integrity is not affected.

---

## Analysis Logic

### Class Overview

**`PolAnalysis`** handles all physics calculations for one detector:
- Stores left/right arm event counts for spin-up and spin-down separately
- Computes simple asymmetry `ε = (nL - nR) / (nL + nR)` per spin state
- Computes cross-ratio asymmetry `ε = (√(nL↑·nR↓) - √(nR↑·nL↓)) / (√(nL↑·nR↓) + √(nR↑·nL↓))`, which cancels detector acceptance asymmetries
- Computes analyzing power `Aᵧ = ε / P_beam` with full error propagation
- Compares to theoretical Aᵧ and reports statistical significance

**`PolHistograms`** handles histogram management:
- 1D energy deposition histograms per scintillator: raw, cut-filtered (both spins), and spin-separated
- 2D correlation histograms between adjacent scintillators
- Drawing methods with cut line overlays and spin comparison

### Inelastic Normalization

For each (angle, parameter) point, elastic and inelastic events are combined into one `PolAnalysis` object. The inelastic tree is truncated to a fraction of its total entries controlled by `kInelFactor = 0.1077`, applied independently to spin-up and spin-down:

```
maxInel = round(kInelFactor × nInel_total)
```

### Scan Loop (`analyze_blpmc_scan`)

For each parameter value and each detector angle the function:

1. Discovers ROOT files by scanning the folder tree — no hardcoded filenames
2. Opens the four trees (elastic spin-up/down, inelastic spin-up/down)
3. Counts events passing scintillator energy cuts via `TTree::Draw`
4. Computes Aᵧ via cross-ratio method
5. Stores the result as one point in a `TGraphErrors`

### Selection Cuts

| Detector | Scintillator | Min [MeV] | Max [MeV] |
|----------|-------------|-----------|-----------|
| P12 (12°) | S0 | 1.5 | 2.5 |
| P12 (12°) | S1 | 3.5 | 5.7 |
| P12 (12°) | S2 | 6.5 | 9.0 |
| P16 (16.2°) | S0 | 1.5 | 2.5 |
| P16 (16.2°) | S1 | 10.0 | 14.0 |
| P16 (16.2°) | S2 | 11.0 | 21.0 |
| P16 (16.2°) | S3 | 15.0 | 35.0 |

Cuts are defined as constants at the top of `analyze_blpmc.C` and can be tuned without touching the class code.

---

## Inspection Mode

After running the scan, any single (angle, parameter) point can be inspected interactively:

```
inspect("16p2", 0.0)      // all spins combined (default)
inspect("16p2", 0.0,  1)  // SpinUp only
inspect("16p2", 0.0, -1)  // SpinDown only
```

Available draw functions after `inspect()`:

| Function | Description |
|----------|-------------|
| `insp_draw_1d_raw()` | Raw energy deposition spectra for all scintillators |
| `insp_draw_1d_comparison()` | Raw vs cut-filtered overlay with cut lines |
| `insp_draw_1d_spin()` | Spin-up (blue) vs spin-down (red) after cuts |
| `insp_draw_2d()` | 2D correlation plots between adjacent scintillators |
| `insp_draw_all()` | All of the above |

The analysis printout (event counts, asymmetries, Aᵧ, comparison to theory) always uses both spin states regardless of the spin filter argument — the filter only affects which histograms are filled.

### Spin Comparison Plot (`insp_draw_1d_spin`)

Shows the energy deposition spectra after cuts separately for spin-up (blue) and spin-down (red) in each scintillator pad. Dashed gray lines mark the cut boundaries for each scintillator. This plot makes the left/right asymmetry between spin states directly visible — the relative height difference between the two curves is proportional to the analyzing power.

---

## Raw Count Summary Plot (`draw_count_vs_param`)

```
draw_count_vs_param()
```

Produces one canvas per detector (P12, P16) with one pad per scintillator, each showing four graphs of raw hit counts (energy deposition > 0) vs the scan parameter:

| Marker | Description |
|--------|-------------|
| Blue solid circle | Left arm, Spin↑ |
| Blue open circle | Left arm, Spin↓ |
| Red solid square | Right arm, Spin↑ |
| Red open square | Right arm, Spin↓ |

Only elastic trees are used (inelastic normalization would distort a raw-count overview). Error bars are Poisson `√N`. This plot is useful for quickly checking detector acceptance and spin asymmetry as a function of the scan parameter before running the full Aᵧ extraction.

Can be called stand-alone at any time — does not require `analyze_blpmc_scan()` to have been run first.

---

## Configuration Reference

All physics constants are defined at the top of `analyze_blpmc.C`:

```cpp
const TString  kScanParam   = "TRANSX"; // scan parameter folder name
const TString  kScanUnit    = "mm";     // axis label unit
const Double_t kInelFactor  = 0.1077;  // inelastic normalization fraction
const Double_t kBeamPol     = 0.80;    // beam polarization
const Double_t kP16_TheoAy  = 0.993;   // theoretical Ay for P16
const Double_t kP12_TheoAy  = 0.786;   // theoretical Ay for P12
```

---

## Known Issues and Fixes

### ROOT 6.32 — `inspect()` segfault on P12 detector

On ROOT 6.32, calling `inspect("12p0", ...)` causes a segmentation fault inside `libCling.so` during the `CountEvents` call. The crash originates in the Cling JIT symbol lookup triggered by the `TTree::Draw` histogram accumulation syntax (`>>+`), which is unstable under heavy interpreter load in ROOT 6.32.

**Fix (applied in `PolAnalysis.h`):** The `>>+` accumulation operator is replaced with a unique histogram name generated per call via a static counter (`_hCntL_0`, `_hCntL_1`, ...). This avoids the JIT symbol lookup entirely while keeping the fast `TTree::Draw` code path. Temporary histograms are deleted immediately after use to prevent memory buildup during long scans.

This fix is transparent — behaviour and results are identical on all ROOT versions.

### Stale compiled objects after header changes

If `PolHistograms.h` or `PolAnalysis.h` are updated, the compiled `.so` files must be removed before starting ROOT, otherwise the old binary layout will be used and the session will crash:

```bash
cd ~/Analysis
rm -f *.so *.pcm *.d
root
```