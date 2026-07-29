# Mathematical Reference: Elastic/Inelastic Weighted Aᵧ

This document collects the exact formulas, error propagation, safeguard
conditions, and code locations for the four core quantities behind the
cross-section-weighted analyzing power calculation in `analyze_blpmc`.

**Notation.** For one reaction sample (elastic or inelastic), the four raw
counts are `nLu, nRu, nLd, nRd` (Left/Right arm, spin Up/Down). `N_cut`
denotes their sum (total events surviving the coincidence cut); `N_gen`
denotes the sample's total generated (pre-cut) event count.

---

## 1. ε_elas and ε_inel — the cross-ratio asymmetry

Computed identically for both samples — `ε_elas` from the elastic
`PolAnalysis` object, `ε_inel` from the inelastic one — using only that
sample's own four counts. No counts are shared between the two.

### Formula

```
A = √(nLu · nRd)
B = √(nRu · nLd)

ε = (A − B) / (A + B)
```

### Error formula

Derived by propagating Poisson statistics (δnᵢ² = nᵢ) directly through all
four counts via the chain rule:

```
δε = √( B²·(nLu + nRd) + A²·(nRu + nLd) )  /  (A + B)²
```

*(Equivalent closed form: δε² = A²B²/(A+B)⁴ · (1/nLu + 1/nRu + 1/nLd + 1/nRd)
— algebraically identical, the code uses the first form.)*

### Safeguard condition

```
if (nLu == 0 || nRu == 0 || nLd == 0 || nRd == 0):
    ε = 0.0
    δε = 0.0
    print "Warning: Zero counts detected, cannot compute cross-ratio!"
```

### Code location

| Quantity | Function | File : Line |
|---|---|---|
| ε | `PolAnalysis::GetCrossRatioAsymmetry()` | `Include/PolAnalysis.h:160` |
| δε | `PolAnalysis::GetCrossRatioError()` | `Include/PolAnalysis.h:182` |

### ⚠️ Potential source of error

The guard fires whenever **any single one** of the four counts is zero —
even if the *other three* have excellent statistics. Mathematically, the
honest result in that situation is **`ε = ±1` with `δε → ∞`** (the formula
genuinely diverges — e.g. if `nLu=0` then `A=0` and `ε=(0−B)/(0+B)=−1`
exactly). Returning `0 ± 0` instead:

- silently discards signal from the three populated bins,
- displays as a clean, confident-looking measurement rather than an
  undetermined one,
- and — critically — this can happen even when the *sample total*
  (`N_cut = nLu+nRu+nLd+nRd`) is large, since the guard checks individual
  bins, not the sum.

This is the most direct way a spurious value can enter `ε_w` (§4).

---

## 2. R_eff — the effective cross-section weight

### Formula

```
eff_elas = N_cut_elas / N_gen_elas
eff_inel = N_cut_inel / N_gen_inel

R_eff = R0 · (eff_inel / eff_elas)
```
`R0` (from `GetInelFactor`) and both `N_gen` values are treated as **exact
constants** — `R0` is a theory input (Meyer cross-section ratio); `N_gen` is
a deterministic row-count of the file, not a statistical measurement. Only
`N_cut_elas` and `N_cut_inel` carry statistical uncertainty.

### Error formula

Since `R_eff ∝ N_cut_inel / N_cut_elas` (a ratio of two independent Poisson
counts, all other factors fixed), relative variances add:

```
δR_eff = R_eff · √( 1/N_cut_inel + 1/N_cut_elas )
```

### Safeguard conditions

```
if eff_elas == 0 (i.e. N_gen_elas == 0):     R_eff = 0.0            [config guard, shouldn't occur in practice]
if eff_elas_value == 0 (N_cut_elas == 0):    R_eff = 0.0            [see effElas > 0.0 check below]
if N_cut_inel == 0  OR  N_cut_elas == 0:     δR_eff = 0.0  (forced)
```
Precisely, from the code: `R = (effElas > 0.0) ? R0·effInel/effElas : 0.0`,
and `dR` is only computed `if (nCutI > 0 && nCutE > 0)`, otherwise `0.0`.

### Code location

| Quantity | Function | File : Line |
|---|---|---|
| R_eff | `ComputeWeightedAy()` | `analyze_blpmc.C:605` |
| δR_eff | `ComputeWeightedAy()` | `analyze_blpmc.C:616–619` |
| N_cut (via `GetCountTotal`) | `PolAnalysis::GetCountTotal()` | `Include/PolAnalysis.h:248` |
| N_gen (via accessor) | `PolAnalysis::GetGeneratedTotal()` | `Include/PolAnalysis.h:68` |

### ⚠️ Potential sources of error

- **`N_cut_inel = 0`**: `R_eff = 0` is the *correct* point estimate (no
  measured inelastic contamination). But `δR_eff = 0` is a **deliberate
  simplification**, not a rigorous result — the mathematically proper
  treatment of "zero observed counts" is a Poisson upper limit, not a flat
  zero uncertainty. Not implemented.
- **`N_cut_elas = 0`**: `R_eff` falls back to `0` via the `effElas > 0.0`
  guard — but this is coincidental, not meaningful, since `ε_elas` is
  independently also `0±0` in that case (§1's guard), so `ε_w` ends up `0`
  regardless of what `R_eff` does. This is really a "no data at all" state,
  not "no inelastic contamination."
- **Low but nonzero `N_cut_inel`** (e.g. 1–2 events): `R_eff` and `δR_eff`
  are both well-defined, but `δR_eff` can be comparable to or exceed
  `R_eff` — this is *expected*, correct low-statistics behavior, not a bug.
- **Independence assumption**: `R_eff` is built from the same underlying
  `N_cut_elas`/`N_cut_inel` totals that `ε_elas`/`ε_inel` are built from
  (just summed differently — individual bins vs. their sum). Treating
  `R_eff`, `ε_elas`, `ε_inel` as statistically independent in §4's error
  propagation ignores this shared-sample covariance. A structural
  approximation, not corrected for.

---

## 3. Why the combination has the form it does

Before the formula: `ε_w` is a **weighted average**, `x_w = (w₁x₁+w₂x₂)/(w₁+w₂)`,
with `w_elas = 1` and `w_inel = R_eff`. The denominator being the *sum of
weights* is not a stylistic choice — it's required so that `x_w` reduces to
the common value when `x₁ = x₂`, and so it correctly limits to `ε_elas` as
`R_eff → 0` and to `ε_inel` as `R_eff → ∞`.

---

## 4. ε_w — the combined weighted asymmetry

### Formula

```
ε_w = (ε_elas + R_eff · ε_inel) / (1 + R_eff)

A_y = ε_w / P_beam
```

### Error formula

Three independent inputs (`P_beam`'s own uncertainty is **not** included —
see below): `ε_elas`, `ε_inel`, `R_eff`. Partial derivatives:

```
∂ε_w/∂ε_elas = 1/(1+R_eff)
∂ε_w/∂ε_inel = R_eff/(1+R_eff)
∂ε_w/∂R_eff  = (ε_inel − ε_elas)/(1+R_eff)²
```

Combined in quadrature:

```
δε_w = √[ δε_elas²/(1+R_eff)²  +  (R_eff·δε_inel/(1+R_eff))²  +  ((ε_inel−ε_elas)·δR_eff/(1+R_eff)²)² ]

δA_y = δε_w / P_beam
```

### Safeguard condition

```
if P_beam ≤ 0:  entire result (Ay, dAy, R_eff, effElas, effInel) = 0
```
No other guard is applied at this level — `ε_w`/`δε_w` simply inherit
whatever `ε_elas`, `ε_inel`, `δε_elas`, `δε_inel`, `R_eff`, `δR_eff` already
are, including any of their own safeguard fallbacks (§1, §2).

### Code location

| Quantity | Function | File : Line |
|---|---|---|
| ε_w, A_y | `ComputeWeightedAy()` | `analyze_blpmc.C:613–614, 631` |
| δε_w, δA_y | `ComputeWeightedAy()` | `analyze_blpmc.C:621–632` |
| Called from scan loop | `analyze_blpmc_scan()` | `analyze_blpmc.C:909` |
| Called from inspection | `inspect()` | `analyze_blpmc.C:1462` |

### ⚠️ Potential sources of error

- **Inherited from §1**: if `ε_elas` or `ε_inel` individually hits the
  per-bin zero-count fallback (`0±0`) while `R_eff` (built from *totals*,
  a different check) is still a normal nonzero value, that spurious `0`
  is blended into `ε_w` with a real, nonzero weight — with no compensating
  mechanism to suppress it (unlike the clean `R_eff=0` case, where
  `δR_eff=0` naturally kills the inelastic contribution).
- **Missing `P_beam` uncertainty**: `GetAnalyzingPowerError_CrossRatio()`
  (used elsewhere in the code) has a branch that folds in
  `beamPolarizationError`; `ComputeWeightedAy()` does not — it divides by
  `beamPol` with no associated error term. Currently harmless since
  `beamPolarizationError` is set to `0.0` everywhere in this pipeline, but
  would silently under-report `δA_y` if that ever changes.
- **Correctness check available**: when `R_eff = 0` exactly, `ε_w` and
  `δε_w` must equal `ε_elas` and `δε_elas` to full numerical precision —
  confirmed against real scan output (ROTY = −1, P16: elastic-only and
  weighted both print `A_y = 0.968026 ± 0.140089`).

---

## Summary: known open items

| # | Item | Where it enters |
|---|---|---|
| 1 | Per-bin zero-count guard returns `0±0` instead of the true `ε=±1, δε→∞` divergence | §1 → propagates into §4 |
| 2 | `δR_eff` forced to `0` at `N_cut_inel=0` instead of a proper Poisson upper limit | §2 |
| 3 | `P_beam` uncertainty not propagated into `δA_y` | §4 |
| 4 | `R_eff` uses total-count guard; `ε` uses per-bin guard — the two can disagree on "no data" | §1 vs §2 |
| 5 | `R_eff` and `ε_elas`/`ε_inel` are not statistically independent (shared counts), but treated as such in the quadrature sum | §4 |