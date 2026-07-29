# blpmc Analysis Framework

ROOT-based analysis framework for the BNL 200 MeV proton-carbon (pC)
polarimeter. Extracts the analyzing power **Aᵧ** from Geant4 Monte Carlo
simulations of elastic *and* inelastic pC scattering, at two detector
angles (P12 = 12.0°, P16 = 16.2°), as a function of a scanned experimental
parameter (e.g. beam rotation/translation, `ROTY`/`TRANSX`). Elastic and
inelastic statistics are always read in full (no subsampling) and combined
into one cross-section-weighted `A_y` per parameter point — see *Analysis
Logic* below for why pooling (rather than averaging two separate
asymmetries) is the physically correct combination.

---

## Typical Usage

```bash
# 1. Go to the Analysis folder
cd ~/Analysis

# 2. Start ROOT (rootlogon.C enables multithreading and compiles the classes)
root

# 3. Load the macro and run the full parameter scan for a given scan
#    directory and beam kinetic energy [MeV]
.L analyze_blpmc.C
analyze_blpmc_scan("ROTY_m0p2_p0p1_0p01_Pin_P80", 199.5)

# 4. Draw results (defaults to the cross-section-weighted A_y)
draw_scan_p16()
draw_scan_p12()
draw_scan_both()
draw_scan_p16("elas")          # elastic-only A_y
draw_scan_p16("inel")          # inelastic-only A_y
draw_scan_p16("all")           # elastic + inelastic + weighted, overlaid

# 5. Draw raw scintillator hit counts vs scan parameter (uses the cache
#    filled by the scan above — no files are reopened)
draw_count_vs_param()
draw_count_vs_param("inel", true)

# 6. Inspect a specific (angle, param) point
inspect("16p2", 0.0, 0, 199.5)     # all spins combined
inspect("12p0", 10.0, 1, 199.5)    # spin-up only
inspect("16p2", 5.0, -1, 199.5)    # spin-down only
insp_draw_all()                     # all inspection plots for that point
```

---

## Repository Structure

```
analyze_blpmc/
├── analyze_blpmc.C          ← main analysis macro (scan, draw, inspect)
├── rootlogon.C               ← enables ROOT multithreading, compiles classes
├── Include/
│   ├── PolAnalysis.h/.C      ← event counting (RDataFrame) + A_y calculations
│   ├── PolHistograms.h       ← histogram filling and drawing
│   ├── Kinematics.C/.h       ← CM↔lab angle conversion, t-Mandelstam
│   ├── MeyerScattering.C/.h  ← Meyer elastic/inelastic cross-section splines
│   ├── DetectorConfig.C/.h   ← detector acceptance geometry
│   ├── PhysicalConstants.h   ← masses, excitation energy, etc.
│   └── *Spline_*.C           ← tabulated cross-section data used by
│                                 MeyerScattering
├── CODE_OVERVIEW.md          ← condensed developer-facing reference
│                                 (file map, physics formulas, function list)
├── RESTART_PROMPT.md         ← running log of historical bugs/decisions;
│                                 read before changing counting/normalization
│                                 logic — several subtle bugs here have
│                                 already been found and fixed once
└── README.md
```

Data lives outside the repository and is linked in via a symlink:

```bash
ln -s /path/to/Data ~/Analysis/Data
```

---

## Multithreading

`PolAnalysis::CountEvents` (the scan's main bottleneck — it opens every
elastic and inelastic ROOT file and counts events passing the scintillator
coincidence cuts) is implemented with `ROOT::RDataFrame` instead of
`TTree::Draw`. RDataFrame automatically parallelizes its `Filter`/`Count`
actions across worker threads once ROOT's implicit multithreading is
enabled — which `rootlogon.C` now does once, globally, at startup:

```cpp
ROOT::EnableImplicitMT();   // uses std::thread::hardware_concurrency() by default
```

No code changes are needed elsewhere — every file opened during a scan or
an `inspect()` call benefits automatically. To control the thread count
explicitly (e.g. on a shared machine), set `BLPMC_NTHREADS` before starting
ROOT:

```bash
BLPMC_NTHREADS=8 root
```

The same RDataFrame-based counting is used for the raw per-scintillator
occupancy counts (`countBranchRaw`) shown by `draw_count_vs_param()`.

---

## Scan Directory Naming

`analyze_blpmc_scan()` and `inspect()` take a scan-directory name (a folder
under `Data/`) whose **first underscore-delimited token** is the scanned
parameter name, e.g.:

| Scan directory example                    | Parameter | Unit  |
|--------------------------------------------|-----------|-------|
| `TRANSX_0_5_0p1_MT_P80_dTh5_dPh5_MEYER`     | `TRANSX`  | mm    |
| `ROTY_m0p2_p0p1_0p01_Pin_P80`               | `ROTY`    | deg   |

The unit is inferred automatically: parameters starting with `TRANS` are
`mm`, everything else is `deg`. This name is also used to decide the beam
model for the cross-section normalization — see below.

Set it explicitly via the first argument, or globally:

```cpp
analyze_blpmc_scan("ROTY_m0p2_p0p1_0p01_Pin_P80", 199.5);
// or
gScanDir = "TRANSX_0_5_0p1_MT_P80_dTh5_dPh5_MEYER";
analyze_blpmc_scan();
```

---

## Data Structure

Simulations are organized by reaction type, spin state, detector angle, and
parameter value:

```
Data/{scanDir}/
└── pC_{Elas|Inel443}_{energy}_{beam}_{SpinUp|SpinDown}_{12p0|16p2}_MEYER/
    └── {paramKey}/
        ├── blpmc_YYYY-MM-DD_HH-MM-SS.root   ← simulation output
        ├── config.conf                        ← per-run configuration
        ├── pluto.root -> ...                  ← symlink to PLUTO input
        └── run.log                             ← Geant4 output log
```

`paramKey` (the parameter-value sub-folder name) is auto-detected in any of
three formats, so old and new datasets coexist without renaming:

| Format                | Example         | Meaning |
|------------------------|-----------------|---------|
| Sign-prefixed          | `p1p5`, `m0p2`  | `p`/`m` sign, `p` as decimal point → `+1.5`, `-0.2` |
| Literal-minus, decimal | `-0p2`, `0p1`   | literal `-`, `p` as decimal point → `-0.2`, `0.1` |
| Bare integer (legacy)  | `0`, `5`        | whole numbers only |

Each ROOT file contains a single `TTree` named `tree` with branches for
energy deposition in each scintillator of both detector arms (P12: 3
scintillators, P16: 4 scintillators). Both detector branches are present in
every file regardless of which detector angle was simulated:

```
eDepP12LS0, eDepP12LS1, eDepP12LS2   ← P12 Left arm
eDepP12RS0, eDepP12RS1, eDepP12RS2   ← P12 Right arm
eDepP16LS0, ..., eDepP16LS3          ← P16 Left arm
eDepP16RS0, ..., eDepP16RS3          ← P16 Right arm
```

---

## Analysis Logic

### Class Overview

**`PolAnalysis`** handles all physics calculations for one detector, one
reaction sample (elastic *or* inelastic) at a time:
- Counts left/right-arm coincidences within the scintillator energy cuts,
  separately for spin-up and spin-down, via `ROOT::RDataFrame` (§ Multithreading)
- Computes the cross-ratio asymmetry
  `ε = (√(nL↑·nR↓) − √(nR↑·nL↓)) / (√(nL↑·nR↓) + √(nR↑·nL↓))`,
  which cancels detector acceptance asymmetries
- Computes analyzing power `Aᵧ = ε / P_beam` with full error propagation
- Compares to theoretical Aᵧ and reports statistical significance

**`PolHistograms`** handles histogram management:
- 1D energy deposition histograms per scintillator: raw, cut-filtered
  (both spins), and spin-separated
- 2D correlation histograms between adjacent scintillators
- Drawing methods with cut-line overlays, spin comparison, and an
  elastic-vs-inelastic overlay that can rescale the inelastic histogram
  for display (see below)

### Full-Statistics Elastic + Inelastic, Combined by Pooled Counts

Earlier versions of this framework went through two prior methods before
arriving at the current one — both discussed here since the reasoning that
ruled each one out is part of why the current approach is correct:

1. **Subsampling + pooled cross-ratio.** Read only a fraction of the
   inelastic tree (chosen to match the cross-section ratio) and pooled
   the resulting counts into one cross-ratio calculation. Problems:
   discarded most of the (expensive to generate) inelastic statistics,
   and pooling *before* rescaling for the cross-section ratio silently
   assumed the ratio didn't need correcting for cut efficiency.
2. **Full-statistics elastic/inelastic, weighted-average of two
   independently-computed cross-ratios**, using an efficiency-corrected
   effective weight `R_eff = R0·(eff_inel/eff_elas)`. This fixed the
   subsampling problem, but introduced a different one: computing
   `ε_elas` and `ε_inel` *separately* and only combining them afterward
   is not mathematically the same as pooling the raw counts and computing
   one cross-ratio — the cross-ratio is a nonlinear (square-root) function
   of the counts, so averaging two nonlinear outputs isn't equivalent to
   evaluating the nonlinear function on pooled inputs, whenever
   `ε_elas ≠ ε_inel`.

**The current approach is pooling — but done correctly.** The key physical
fact that makes pooling valid here: elastic and inelastic events are
**indistinguishable once they've survived the coincidence cuts** — a real
detector just registers "a hit passed," with no per-event tag of which
reaction produced it. So the combined asymmetry genuinely *should* be
computed from pooled counts, exactly like a real detector would see them
— the earlier "weighted-average-of-two-cross-ratios" approach was the
artificial one, relying on a separation that only exists in simulation.

1. Read **every** elastic event and **every** inelastic event — no
   subsampling. `ε_elas`, `ε_inel` are still computed independently and
   plotted as separate diagnostic curves (`draw_scan_p16("elas")` /
   `("inel")`), but they no longer feed the weighted combination directly.
2. Before pooling, rescale the inelastic counts by `R0` to correct for the
   fact that elastic and inelastic samples were generated with the *same
   arbitrary* statistics (~10M events each in this setup) regardless of
   their different physical cross sections:

   ```
   Lu = Lu_elas + R0·Lu_inel      (and similarly Ru, Ld, Rd)

   A = √(Lu·Rd),  B = √(Ru·Ld),  S = A+B
   ε_total = (A−B)/S
   A_y     = ε_total/P_beam
   ```
   `Lu/Ru/Ld/Rd` are the four raw left/right-arm, spin-up/down counts —
   `R0` scales the *already cut* inelastic counts directly, which is
   sufficient on its own: since passing the cuts is a per-event property
   independent of how many total events were generated, scaling the final
   surviving count by `R0` gives the same expectation as generating the
   inelastic sample at its true relative rate in the first place, then
   cutting. No separate cut-efficiency correction is needed.

3. Error propagation follows directly, treating `R0` as an exact constant
   (no associated uncertainty) and each raw count as Poisson-distributed:

   ```
   Var(Lu) = Lu_elas + R0²·Lu_inel     (and similarly Var(Ru), Var(Ld), Var(Rd))

   δε_total = (A·B/S²) · √[ Var(Lu)/Lu² + Var(Ru)/Ru² + Var(Ld)/Ld² + Var(Rd)/Rd² ]
   ```
   This is the direct generalization of the plain cross-ratio error
   (`PolAnalysis::GetCrossRatioError()`): setting `R0=1` (no elastic/
   inelastic split) makes `Var(Lu)=Lu` exactly, and the formula reduces
   to the original cross-ratio error formula algebraically. See
   `ComputeWeightedAy()` in `analyze_blpmc.C` for the full implementation.

**Limiting behavior**, both confirmed against real scan output:
   - If this geometry rejects essentially all inelastic events
     (`Lu_inel ≈ Ru_inel ≈ Ld_inel ≈ Rd_inel ≈ 0`), the pooled counts
     reduce to the elastic-only counts and `ε_total → ε_elas` exactly —
     as physically expected when the detector cleanly separates the two
     channels.
   - If inelastic genuinely contributes, its pooled contribution scales
     with `R0` and the real cut-survival counts at that specific
     geometry — not a fixed nominal assumption independent of the scan
     parameter.

`R0` is computed from the Meyer elastic/inelastic cross-section splines,
integrated over the detector acceptance (see `IntegrateXS_local` /
`GetInelFactor` in `analyze_blpmc.C`), and depends on:
- the detector angle (12.0° or 16.2°),
- the beam kinetic energy, which must be supplied explicitly as the
  second argument to `analyze_blpmc_scan()` / `inspect()`, and
- the beam model — `Pin` (perfectly-collimated beam, uses the differential
  cross section at the detector-centre angle) vs `MT`/realistic beam
  (integrates over an acceptance window whose half-opening, in mrad, is
  read from a `_dThNN_` token in the scan directory name).

`R0` can be overridden directly at runtime:

```cpp
gInelFactor_P16 = 0.1077;   // force a fixed value for P16
gInelFactor_P12 = 0.0;      // disable the inelastic contribution for P12
```

### Beam Polarization Auto-Detection

Beam polarization `P_beam` is **no longer a hardcoded constant**. It is
auto-parsed from a `_PXX_` token in the scan directory name, mirroring the
`_dThNN_` window-parsing pattern used for `R0` above:

- `ParseBeamPolarization(scanDir)` scans for `"_P"` immediately followed by
  a digit — this naturally skips `"_Pin_"` (`'i'` is not a digit) and any
  `"_dPhNN_"` token (no underscore directly precedes that `'P'`) — and
  returns the percentage as a fraction, e.g. `..._P80_...` → `0.80`.
- `GetBeamPolarization(scanDir)` caches the resolved value in the global
  `gBeamPol` (reset to `-1.0` at the top of `analyze_blpmc_scan()`, the same
  convention as `gInelFactor_P16`/`gInelFactor_P12`), and falls back to
  `kBeamPol_Fallback = 0.80` **with a loud printed warning** if no such
  token is found in the folder name.
- The resolved value is threaded through explicitly: `RunCounting()` takes
  `beamPol` as an argument, and both `analyze_blpmc_scan()` and `inspect()`
  resolve it once via `GetBeamPolarization(gScanDir)` before passing it
  down — no code path reads a bare constant anymore.

Just like `R0`, it can be overridden directly at runtime:

```cpp
gBeamPol = 0.75;   // force a fixed polarization, bypassing folder-name parsing
```

> **Not yet propagated:** `beamPolarizationError` is still hardcoded to
> `0.0` everywhere — polarization uncertainty is not currently folded into
> any `A_y` error bar. This is a known, deliberately deferred limitation
> (see *Known Limitations* below), not a bug.

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

Cuts are constants at the top of `analyze_blpmc.C` and can be tuned
without touching the class code.

### Scan Loop (`analyze_blpmc_scan`)

```cpp
analyze_blpmc_scan(userScanDir = "", beamEnergy = 200.0,
                    startVal = -1e30, stopVal = 1e30)
```

`startVal`/`stopVal` restrict which parameter sub-folders are processed
(inclusive bounds on the parsed numeric parameter value). For each
parameter value and each detector angle the function:

1. Discovers ROOT files by scanning the folder tree — no hardcoded filenames
2. Runs two independent counting passes (elastic, inelastic), each reading
   its sample in full, multithreaded via RDataFrame
3. Computes `ε_elas`, `ε_inel`, and the cross-section-weighted `ε_w` → `A_y`
4. Stores each as one point in its own `TGraphErrors`
   (`gGraph_{P16,P12}_{elas,inel,weighted}`)
5. Also caches raw and full-cut scintillator counts so
   `draw_count_vs_param()` can plot them without reopening any files

---

## Inspection Mode

After running the scan (or standalone), any single (angle, parameter)
point can be inspected interactively:

```cpp
inspect("16p2", 0.0, 0, 199.5)      // all spins combined (default)
inspect("16p2", 0.0, 1, 199.5)      // SpinUp only
inspect("16p2", 0.0, -1, 199.5)     // SpinDown only
```

Elastic and inelastic events are read in full (same as the scan) and
tracked in **separate** `PolAnalysis` objects, so the printed summary
reports the elastic-only, inelastic-only, and cross-section-weighted `A_y`
side by side — consistent with what the scan computes.

Available draw functions after `inspect()`:

| Function | Description |
|----------|-------------|
| `insp_draw_1d_raw()` | Raw energy deposition spectra for all scintillators |
| `insp_draw_1d_comparison()` | Raw vs cut-filtered overlay with cut lines |
| `insp_draw_1d_spin()` | Spin-up (blue) vs spin-down (red) after cuts |
| `insp_draw_2d()` | 2D correlation plots between adjacent scintillators |
| `insp_draw_elas_inel()` | Elastic vs inelastic RAW spectra, inelastic scaled to R |
| `insp_draw_elas_inel_cuts()` | Same, but post-cut spectra |
| `insp_draw_all()` | All of the above |

### Elastic vs Inelastic Overlay (`insp_draw_elas_inel[_cuts]`)

Since both samples are now read in full, their raw event counts no longer
reflect the physical cross-section ratio (elastic and inelastic MC samples
aren't necessarily generated in physically-proportional numbers). For
**display purposes only**, the inelastic histogram in this overlay is
independently rescaled per scintillator so its integral relative to the
elastic one equals `R` — the same ratio used in the `A_y` combination. The
underlying stored histograms (and all counts used for physics results) are
never modified by this — only the plotted clone is scaled.

### Spin Comparison Plot (`insp_draw_1d_spin`)

Shows the energy deposition spectra after cuts separately for spin-up
(blue) and spin-down (red) in each scintillator pad. Dashed gray lines mark
the cut boundaries. The relative height difference between the two curves
is proportional to the analyzing power.

---

## Raw Count Summary Plot (`draw_count_vs_param`)

```cpp
draw_count_vs_param()                  // elastic raw only (default)
draw_count_vs_param("elas", true)      // elastic raw + elastic full-cut, overlaid
draw_count_vs_param("inel")            // inelastic raw only (full statistics)
draw_count_vs_param("inel", true)      // inelastic raw + inelastic full-cut, overlaid
draw_count_vs_param("both", true)      // all four curve groups, overlaid (16 curves/pad)
```

Produces one canvas per detector (P12, P16) with one pad per scintillator.
`mode` selects which reaction(s) to show — `"elas"`, `"inel"`, or `"both"` —
and controls **both** the raw and (if `showCuts=true`) the full-cut curves
for that reaction, so raw and cut can be compared directly for the same
reaction on the same pad:

- **Raw** (branch > 0, per scintillator): shows how many hits each
  individual scintillator registers on its own — no coincidence
  requirement.
- **Full-cut** (`showCuts=true`): the same full multi-scintillator
  coincidence count that actually feeds the A_y calculation. This curve
  is identical across every scintillator pad for a given reaction — a
  coincidence is a single per-event yes/no, not a per-scintillator
  quantity — it's shown on each pad only for side-by-side reference
  against that scintillator's raw occupancy.

Comparing raw vs full-cut for the *same* scintillator and reaction shows
how much of that scintillator's raw occupancy actually survives the full
coincidence cut. Comparing elastic vs inelastic (`mode="both"`) shows
whether inelastic events leak into the elastic selection window. All
inelastic curves (raw and full-cut) reflect the full inelastic statistics
— no subsampling.

Reads entirely from the cache filled by `analyze_blpmc_scan()` — no files
are reopened. Run the scan at least once first.

---

## Configuration Reference

All physics constants are defined at the top of `analyze_blpmc.C`:

```cpp
const Double_t kBeamPol_Fallback = 0.80;   // used ONLY if "_PXX_" can't be
                                            // parsed from the scan dir name
                                            // (see Beam Polarization Auto-
                                            // Detection above) — never read
                                            // directly by any counting code
const Double_t kP16_TheoAy = 0.993;        // theoretical Ay for P16
const Double_t kP12_TheoAy = 0.786;        // theoretical Ay for P12
// + per-scintillator energy cuts for P12 and P16 (see Selection Cuts above)
```

---

## Known Limitations (not yet fixed)

1. **Beam-type auto-detection silently defaults on unrecognised folder
   names.** Whether a scan directory is treated as a `Pin` (perfectly
   collimated) or `MT`/realistic beam is inferred purely from substrings in
   its name (`_Pin_` vs `_dThNN_`). If a folder has **neither** token, the
   code falls back to `MT` + 5 mrad with only a printed warning — no hard
   error. This has been hit once in practice (a folder named only after its
   scan parameter, with no beam tokens); the fix so far has been to rename
   the folder consistently rather than change the code. Revisit if it
   recurs.
2. **`beamPolarizationError` is hardcoded to `0.0` throughout.** Beam
   polarization is now auto-parsed (see *Beam Polarization Auto-Detection*
   above), but its uncertainty is not yet propagated into any `A_y` error
   bar. Deliberately deferred — a separate task if/when the polarization
   uncertainty needs to be included in the final error budget.

---

## Environment / Workflow

- Remote analysis machine: `ssh -Y -p7777 acannavo@210.119.41.60`
  (ROOT 6.32.12), working directory `~/Analysis`, data linked in at
  `~/Analysis/Data` (symlink to the real data location — see *Repository
  Structure* above).
- Local edits are **not** automatically visible on the remote machine —
  sync explicitly via `rsync`/`scp` (see `sync_analysis.sh`), excluding
  compiled ACLiC artifacts (`*.so`, `*.d`, `*_ACLiC_dict*`, `*_rdict*`) so
  stale binaries are never carried across.
- `rootlogon.C` compiles/loads everything in a specific dependency order —
  keep this order when adding new headers:
  1. `PhysicalConstants.h+`
  2. `DetectorConfig.h+` / `DetectorConfig.C+`
  3. `Kinematics.C` (interpreted, **no** `+` — ACLiC has trouble with it)
  4. `MeyerScattering.C` (interpreted, **no** `+` — it `#include`s spline
     data files that ACLiC can't resolve)
  5. `PolHistograms.h+`
  6. `PolAnalysis.h+`
  7. `analyze_blpmc.C` (interpreted)

---

## Known Issues and Fixes

### Stale compiled objects after header changes

If `PolHistograms.h` or `PolAnalysis.h` are updated, the compiled `.so`
files must be removed before starting ROOT, otherwise the old binary
layout will be used and the session will crash:

```bash
cd ~/Analysis
rm -f *.so *.pcm *.d
root
```

### ROOT 6.32 `inspect()` segfault on P12 — resolved

Earlier versions counted events with `TTree::Draw`'s `>>+` histogram
accumulation syntax, which triggered a Cling JIT symbol lookup that
crashed on ROOT 6.32 under heavy interpreter load. `PolAnalysis::CountEvents`
no longer uses `TTree::Draw` at all — it counts via `ROOT::RDataFrame` —
so this class of crash no longer applies.