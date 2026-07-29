// =========================================================
// analyze_blpmc.C — Parameter Scan
//
// Automatically scans Data/{gScanDir}/ folder structure:
//   {gScanDir}/
//     pC_{Elas|Inel443}_{energy}_{beam}_{SpinUp|SpinDown}_{12p0|16p2}_MEYER/
//       {p0p0 .. pNpN}/          ← new: p/m prefix for ±, 'p' as decimal point
//         blpmc_*-*-*_*-*-*.root
//
// Filename convention for parameter values (as of new naming scheme):
//   p0p1  = +0.1    (p = plus,  'p' = decimal point)
//   m0p2  = -0.2    (m = minus, 'p' = decimal point)
//   p1p5  = +1.5
//   m10p0 = -10.0
//
// Usage:
//   analyze_blpmc_scan("ROTY_m0p2_p0p1_0p01_Pin_P80", 199.5)
//   analyze_blpmc_scan("", 199.5, -0.1, 0.1)       // only process points in [-0.1, 0.1]
//   inspect("12p0", -0.1, 0, 199.5)
//   draw_scan_p16()               // weighted A_y (default, physically correct)
//   draw_scan_p16("elas")         // elastic-only A_y
//   draw_scan_p16("inel")         // inelastic-only A_y
//   draw_scan_p16("weighted")     // weighted combination (same as default)
//   draw_scan_p16("all")          // all three modes on one canvas
//   draw_count_vs_param()                  // elastic raw only (cached, no re-reading)
//   draw_count_vs_param("inel", true)      // inelastic raw + inelastic full-cut (cached)
//
// IMPORTANT: draw_count_vs_param() reads a cache filled by
// analyze_blpmc_scan() — run the scan at least once first.
//
// Multithreading:
//   PolAnalysis::CountEvents uses ROOT::RDataFrame, which is
//   automatically parallelized once ROOT::EnableImplicitMT() has
//   run — this happens once, globally, in rootlogon.C. No action
//   is needed here; every scan point benefits automatically.
//
// Statistics and normalization:
//   Elastic and inelastic samples are ALWAYS read in full — no
//   event subsampling. Each population's analyzing power (A_y) is
//   measured independently at full statistical precision, then the
//   two are combined via a cross-section-weighted average:
//
//     R    = inel/elas weight (physics cross-section ratio, see below)
//     ε_w  = (ε_elas + R·ε_inel) / (1 + R)
//     A_y  = ε_w / P_beam
//
//   The inel/elas weight R is computed from the Meyer cross sections
//   integrated over the detector acceptance (same physics as
//   ComputeNormalization.C). It depends on the detector angle
//   (12.0 or 16.2 deg) and the beam kinetic energy, which you must
//   supply explicitly as the second argument:
//     analyze_blpmc_scan("ROTY_m0p2_p0p1_0p01", 199.5)
//     inspect("16p2", 0.0, 0, 199.5)
//   You can also override the factors directly at runtime:
//     gInelFactor_P16 = 0.1077;   // force a fixed value for P16
//     gInelFactor_P12 = 0.0;      // disable inelastic for P12
// =========================================================

#include <map>
#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"

// =========================================================
// Configuration — physics constants
// =========================================================

// P16 (16.2 deg, 4 scintillators) energy cuts [MeV]
const Double_t kP16_S0_min = 1.5,  kP16_S0_max = 2.5;
const Double_t kP16_S1_min = 10.0, kP16_S1_max = 14.0;
const Double_t kP16_S2_min = 11.0, kP16_S2_max = 21.0;
const Double_t kP16_S3_min = 15.0, kP16_S3_max = 35.0;

// P12 (12.0 deg, 3 scintillators) energy cuts [MeV]
const Double_t kP12_S0_min = 1.5,  kP12_S0_max = 2.5;
const Double_t kP12_S1_min = 3.5,  kP12_S1_max = 5.7;
const Double_t kP12_S2_min = 6.5,  kP12_S2_max = 9.0;

// Beam polarization and theoretical Ay
// kBeamPol_Fallback is used ONLY if a scan directory name has no
// parseable "_PXX_" token (see ParseBeamPolarization/GetBeamPolarization
// below) — this mirrors the ParseThetaWindow "5 mrad fallback" convention.
// Do not use this constant directly anywhere else; always go through
// GetBeamPolarization(scanDir) so the value is auto-parsed per scan.
const Double_t kBeamPol_Fallback = 0.80;
const Double_t kP16_TheoAy = 0.993;
const Double_t kP12_TheoAy = 0.786;

// =========================================================
// Runtime globals
// =========================================================
TString gScanDir = "";

// Physics-based inel/elas weight (computed per detector angle).
// Negative value = "not yet computed, use IntegrateXS".
// Set to a positive value to override.
Double_t gInelFactor_P16 = -1.0;
Double_t gInelFactor_P12 = -1.0;

// Beam polarization, auto-parsed from a "_PXX_" token in the scan
// directory name (see ParseBeamPolarization/GetBeamPolarization).
// Negative value = "not yet resolved for this scan".
// Set to a positive value (0-1) to override, exactly like
// gInelFactor_P16/P12 can be overridden.
Double_t gBeamPol = -1.0;

// =========================================================
// Cached scan results — filled by analyze_blpmc_scan,
// consumed by draw_scan_* without re-reading any files.
//
// Three TGraphErrors per detector: elas-only, inel-only, weighted
// (the physically-correct cross-section-weighted combination).
// =========================================================
TGraphErrors *gGraph_P16_elas     = nullptr;
TGraphErrors *gGraph_P16_inel     = nullptr;
TGraphErrors *gGraph_P16_weighted = nullptr;   // cross-section-weighted combination

TGraphErrors *gGraph_P12_elas     = nullptr;
TGraphErrors *gGraph_P12_inel     = nullptr;
TGraphErrors *gGraph_P12_weighted = nullptr;   // cross-section-weighted combination

// =========================================================
// Cached COUNT results — filled by analyze_blpmc_scan,
// consumed by draw_count_vs_param without re-reading any files.
//
// Layout per detector: vector<TGraphErrors*> of size nScint * 16.
// For scintillator `is`, slot `ig` (0-15) lives at index is*16+ig:
//   0-3   : raw      elas  L-Up, L-Down, R-Up, R-Down   (per-scintillator)
//   4-7   : full-cut elas  L-Up, L-Down, R-Up, R-Down   (same for all is)
//   8-11  : raw      inel  L-Up, L-Down, R-Up, R-Down   (per-scintillator, full statistics)
//   12-15 : full-cut inel  L-Up, L-Down, R-Up, R-Down   (same for all is, full statistics)
// =========================================================
std::vector<TGraphErrors*> gCountGraphs_P16;
std::vector<TGraphErrors*> gCountGraphs_P12;
Bool_t gCountsValid_P16 = false;
Bool_t gCountsValid_P12 = false;

// =========================================================
// Global inspection objects (reused across inspect() calls)
// =========================================================
PolHistograms *gInsp_hist      = nullptr;
PolAnalysis   *gInsp_ana       = nullptr;   // elastic counts + cuts (also used by draw fns for cut lines)
PolAnalysis   *gInsp_ana_inel  = nullptr;   // inelastic counts, tracked separately (see inspect())
TString        gInsp_label = "";
Double_t       gInsp_R     = -1.0;          // inel/elas cross-section ratio from the last inspect() call

// Separate elastic-only / inelastic-only histogram sets, used by
// insp_draw_elas_inel() to visually distinguish which population
// contributes inside the cut window — a single combined histogram
// (gInsp_hist above) cannot show this since elastic and inelastic
// are filled into the same bins.
PolHistograms *gInsp_hist_elas = nullptr;
PolHistograms *gInsp_hist_inel = nullptr;


// =========================================================
// ─── SECTION 1: VALUE PARSING ────────────────────────────
//
// Handles the OLD format  (no sign prefix):
//   "0p1"  → 0.1    "1p5" → 1.5    "2" → 2.0
//
// AND the NEW format (p = plus, m = minus):
//   "p0p1" → +0.1   "m0p2" → -0.2
//   "p1p5" → +1.5   "m10p0" → -10.0
//
// Note: in new format the FIRST character is p/m (sign),
//       then the 'p' inside the number is the decimal point.
// =========================================================
Double_t ParsePValue(const TString &s) {
    if (s.Length() == 0) return 0.0;

    char first = s[0];

    // ── New format with letter sign: 'p' (plus) or 'm' (minus) prefix ──
    // e.g. "p0p1" = +0.1, "m0p2" = -0.2, "p5p0" = +5.0
    if (first == 'p' || first == 'm') {
        Double_t sign = (first == 'm') ? -1.0 : +1.0;
        TString  body = s(1, s.Length() - 1);
        body.ReplaceAll("p", ".");
        return sign * body.Atof();
    }

    // ── Literal-minus format: '-' prefix, 'p' as decimal separator ──
    // e.g. "-0p1" = -0.1, "-1" = -1.0, "-1p5" = -1.5
    if (first == '-') {
        TString body = s(1, s.Length() - 1);
        body.ReplaceAll("p", ".");
        return -body.Atof();
    }

    // ── Old format (no sign character): 'p' as decimal separator ──
    // e.g. "0p1" = 0.1, "5" = 5.0, "0p0" = 0.0
    TString tmp = s;
    tmp.ReplaceAll("p", ".");
    return tmp.Atof();
}

// =========================================================
// Format a Double_t back to the NEW folder-name convention
//   -0.2  → "m0p2"
//   +1.5  → "p1p5"
//   +5.0  → "p5p0"   (or "5" in old format — handled by scan)
// Used by inspect() to build the folder path when the scan
// uses the new naming scheme.
// =========================================================
TString FormatPValue(Double_t v) {
    TString sign = (v < 0) ? "m" : "p";
    TString body = Form("%.10g", TMath::Abs(v));
    body.ReplaceAll(".", "p");
    return sign + body;
}

// =========================================================
// ─── SECTION 2: DETECTOR CONFIG ──────────────────────────
// =========================================================
void ParseScanDir(const TString &dir, TString &param, TString &unit) {
    Ssiz_t pos = dir.First('_');
    param = (pos == kNPOS) ? dir : dir(0, pos);
    unit  = param.BeginsWith("TRANS") ? "mm" : "deg";
}

Bool_t ResolveScanDir(TString &scanDir, TString &param, TString &unit,
                      const TString &userArg = "") {
    if (userArg != "") gScanDir = userArg;
    if (gScanDir == "") {
        cout << "[Error] No scan directory set." << endl;
        cout << "  Pass it as argument:  analyze_blpmc_scan(\"ROTY_m0p2_p0p1_0p01\")" << endl;
        cout << "  Or set globally:      gScanDir = \"TRANSX_0_5_0p1\"" << endl;
        return false;
    }
    ParseScanDir(gScanDir, param, unit);
    scanDir = "Data/" + gScanDir;
    return true;
}

// NOTE: Beam energy is supplied explicitly by the user as a function
// argument to analyze_blpmc_scan() and inspect().  See usage above.


// =========================================================
// ─── SECTION 3: PHYSICS-BASED NORMALIZATION ──────────────
//
// Two modes, selected automatically from the scan directory name:
//
//  PIN BEAM  (folder contains "_Pin_"):
//    ΔΘ = ΔΦ = 0 — all events hit exactly the detector centre.
//    Ratio = dσ_inel(θ_centre) / dσ_elas(θ_centre)  [point ratio].
//
//  MT / REALISTIC BEAM  (no "_Pin_" in folder name):
//    Ratio = ∫ dσ_inel dΩ / ∫ dσ_elas dΩ
//    integrated over [θ_centre ± DETECTOR_THETA_WINDOW].
//    Identical to ComputeNormalization.C.
//
// Requires PEG functions loaded by rootlogon.C:
//   ComputeCMAngleRange, ConvertThetaCMtoT,
//   MeyerXS_Elastic, MeyerXS_Inelastic
// =========================================================

namespace NormConfig {
    const Double_t kThetaP16_rad = 16.2 * TMath::DegToRad();
    const Double_t kThetaP12_rad = 12.0 * TMath::DegToRad();
}

Bool_t IsPinBeam(const TString &scanDir) {
    return scanDir.Contains("_Pin_") || scanDir.Contains("_pin_");
}

// Parse the theta window [rad] from the folder name.
// Looks for "_dThXX_" token where XX is the half-opening in mrad.
// Examples:
//   TRANSY_0_5_0p1_MT_P80_dTh5_dPh5_MEYER   → 5e-3 rad
//   TRANSY_0_5_0p1_MT_P80_dTh20_dPh20_MEYER  → 20e-3 rad
// Returns -1.0 if the token is not found (caller should handle).
Double_t ParseThetaWindow(const TString &scanDir) {
    // Find "_dTh" token
    Ssiz_t pos = scanDir.Index("_dTh");
    if (pos == kNPOS) {
        // Try lowercase
        pos = scanDir.Index("_dth");
        if (pos == kNPOS) return -1.0;
    }
    // Skip past "_dTh"
    Ssiz_t start = pos + 4;
    // Read digits
    TString numStr = "";
    while (start < scanDir.Length() && isdigit(scanDir[start])) {
        numStr += scanDir[start];
        start++;
    }
    if (numStr.Length() == 0) return -1.0;
    return numStr.Atof() * 1e-3;   // mrad → rad
}

// Parse the beam polarization [fraction 0-1] from a "_PXX_" token in
// the scan directory name (XX = percent). Mirrors ParseThetaWindow's
// "_dThXX_" parsing exactly, one token, one parser.
// Examples:
//   ROTY_m0p2_p0p1_0p01_Pin_P80           → 0.80
//   TRANSY_0_5_0p1_MT_P70_dTh5_dPh5_MEYER → 0.70
// Scans for "_P" followed immediately by a digit — this naturally
// skips "_Pin_" ('i' is not a digit) and any "_dPh.." token (no
// underscore directly precedes the 'P' there). Returns -1.0 if no
// such token is found anywhere in the name (caller should handle,
// same convention as ParseThetaWindow's -1.0 sentinel).
Double_t ParseBeamPolarization(const TString &scanDir) {
    Ssiz_t pos = 0;
    while (true) {
        pos = scanDir.Index("_P", pos);
        if (pos == kNPOS) return -1.0;

        Ssiz_t start = pos + 2;
        if (start < scanDir.Length() && isdigit(scanDir[start])) {
            TString numStr = "";
            Ssiz_t i = start;
            while (i < scanDir.Length() && isdigit(scanDir[i])) {
                numStr += scanDir[i];
                i++;
            }
            return numStr.Atof() / 100.0;   // percent → fraction
        }
        // Not a match at this position (e.g. "_Pin_") — keep looking.
        pos += 2;
    }
}

// Resolve and cache the beam polarization for the CURRENT scan.
// Cache (gBeamPol) is reset to -1.0 at the top of analyze_blpmc_scan()
// so it is re-parsed per scan, exactly like gInelFactor_P16/P12.
// Falls back to kBeamPol_Fallback with a loud warning if the scan
// directory name has no parseable "_PXX_" token — this is the same
// silent-danger case flagged in RESTART_PROMPT.md/CODE_OVERVIEW.md
// §6.1, now made loud instead of silent.
Double_t GetBeamPolarization(const TString &scanDir) {
    if (gBeamPol > 0.0) {
        cout << Form("  [Norm] Using beam polarization P = %.4f (cached/override)", gBeamPol) << endl;
        return gBeamPol;
    }

    Double_t parsed = ParseBeamPolarization(scanDir);
    if (parsed <= 0.0) {
        cout << Form("[Warning] Could not parse beam polarization (\"_PXX_\") from \"%s\""
                     " — using %.2f fallback. Check the folder name!",
                     scanDir.Data(), kBeamPol_Fallback) << endl;
        gBeamPol = kBeamPol_Fallback;
    } else {
        cout << Form("  [Norm] Beam polarization parsed from folder name: P = %.4f", parsed) << endl;
        gBeamPol = parsed;
    }
    return gBeamPol;
}

// Trapezoidal integration of dσ/dΩ × sin(θ_CM) over CM angle range [deg].
Double_t IntegrateXS_local(Double_t theta_cm_min_deg,
                            Double_t theta_cm_max_deg,
                            Double_t ekin,
                            Bool_t   inelastic,
                            Int_t    n_steps = 10000) {
    Double_t dtheta = (theta_cm_max_deg - theta_cm_min_deg) / n_steps;
    Double_t sum    = 0.;
    for (Int_t i = 0; i <= n_steps; i++) {
        Double_t theta_deg = theta_cm_min_deg + i * dtheta;
        Double_t theta_rad = theta_deg * TMath::DegToRad();
        Double_t t         = ConvertThetaCMtoT(theta_deg, ekin);
        Double_t xs        = inelastic ? MeyerXS_Inelastic(t, ekin)
                                       : MeyerXS_Elastic  (t, ekin);
        Double_t w = (i == 0 || i == n_steps) ? 0.5 : 1.0;
        sum += w * xs * TMath::Sin(theta_rad) * dtheta * TMath::DegToRad();
    }
    return 2. * TMath::Pi() * sum;   // [mb]
}

// Compute and cache the inel/elas weight for one detector.
// scanDir is used to detect Pin vs MT beam.
Double_t GetInelFactor(Bool_t isP16, Double_t ekin, const TString &scanDir) {
    Double_t &cache = isP16 ? gInelFactor_P16 : gInelFactor_P12;

    if (cache > 0.0) {
        cout << Form("  [Norm] Using %s inel factor = %.6f (cached/override)",
                     isP16 ? "P16" : "P12", cache) << endl;
        return cache;
    }

    Double_t theta_center_rad = isP16 ? NormConfig::kThetaP16_rad
                                      : NormConfig::kThetaP12_rad;
    Double_t theta_center_deg = theta_center_rad * TMath::RadToDeg();

    if (IsPinBeam(scanDir)) {
        // ── Pin beam: point ratio at detector centre ──────────────────
        // Use ComputeCMAngleRange with zero window to get the CM angle
        // corresponding to the lab centre angle.
        Double_t theta_cm_ctr = 0., dummy = 0.;
        {
            std::streambuf *old = std::cout.rdbuf();
            std::ostringstream devnull;
            std::cout.rdbuf(devnull.rdbuf());
            ComputeCMAngleRange(ekin, theta_center_rad, 0.0, theta_cm_ctr, dummy);
            std::cout.rdbuf(old);
        }
        Double_t t_ctr   = ConvertThetaCMtoT(theta_cm_ctr, ekin);
        Double_t xs_elas = MeyerXS_Elastic  (t_ctr, ekin);
        Double_t xs_inel = MeyerXS_Inelastic(t_ctr, ekin);

        if (xs_elas <= 0.) {
            cout << "[Warning] Zero elastic XS at centre — using legacy factor 0.1077" << endl;
            cache = 0.1077;
            return cache;
        }
        cache = xs_inel / xs_elas;
        cout << Form("  [Norm] %s  E=%.2f MeV  Pin beam  θ_lab=%.1f°  θ_CM=%.3f°",
                     isP16 ? "P16" : "P12", ekin, theta_center_deg, theta_cm_ctr) << endl;
        cout << Form("         dσ_elas=%.4f mb/sr  dσ_inel=%.4f mb/sr  R=%.6f",
                     xs_elas, xs_inel, cache) << endl;

    } else {
        // ── MT beam: integrate over acceptance window read from folder name ──
        Double_t thetaWindow = ParseThetaWindow(scanDir);
        if (thetaWindow < 0.) {
            cout << Form("[Warning] Could not parse dTh from \"%s\" — using 5 mrad fallback",
                         scanDir.Data()) << endl;
            thetaWindow = 5e-3;   // safe fallback; warn loudly so user notices
        }
        Double_t theta_cm_min = 0., theta_cm_max = 0.;
        {
            std::streambuf *old = std::cout.rdbuf();
            std::ostringstream devnull;
            std::cout.rdbuf(devnull.rdbuf());
            ComputeCMAngleRange(ekin, theta_center_rad, thetaWindow,
                                theta_cm_min, theta_cm_max);
            std::cout.rdbuf(old);
        }
        Double_t sigma_elas = IntegrateXS_local(theta_cm_min, theta_cm_max, ekin, false);
        Double_t sigma_inel = IntegrateXS_local(theta_cm_min, theta_cm_max, ekin, true);

        if (sigma_elas <= 0.) {
            cout << "[Warning] IntegrateXS returned zero elastic XS — using legacy factor 0.1077" << endl;
            cache = 0.1077;
            return cache;
        }
        cache = sigma_inel / sigma_elas;
        cout << Form("  [Norm] %s  E=%.2f MeV  MT beam  ΔΘ=%.0f mrad  θ_CM=[%.3f,%.3f]deg",
                     isP16 ? "P16" : "P12", ekin, thetaWindow*1e3, theta_cm_min, theta_cm_max) << endl;
        cout << Form("         σ_elas=%.6f mb  σ_inel=%.6f mb  R=%.6f",
                     sigma_elas, sigma_inel, cache) << endl;
    }
    return cache;
}


// =========================================================
// ─── SECTION 4: FILE UTILITIES ───────────────────────────
// =========================================================
TString FindOutputRoot(const TString &runDir) {
    TSystemDirectory dir(runDir.Data(), runDir.Data());
    TList *files = dir.GetListOfFiles();
    if (!files) return "";
    TSystemFile *f;
    TIter next(files);
    while ((f = (TSystemFile*)next())) {
        TString fname = f->GetName();
        if (fname.EndsWith(".root") && fname.BeginsWith("blpmc_"))
            return runDir + "/" + fname;
    }
    return "";
}

TTree* OpenTree(const TString &rootFile, TFile* &fileHandle) {
    fileHandle = TFile::Open(rootFile);
    if (!fileHandle || fileHandle->IsZombie()) {
        cout << "[Error] Cannot open: " << rootFile << endl;
        fileHandle = nullptr;
        return nullptr;
    }
    TTree *t = (TTree*)fileHandle->Get("tree");
    if (!t) {
        cout << "[Error] TTree 'tree' not found in: " << rootFile << endl;
        fileHandle->Close();
        fileHandle = nullptr;
        return nullptr;
    }
    cout << "  Opened: " << rootFile
         << "  (" << t->GetEntries() << " entries)" << endl;
    return t;
}

void SafeClose(TFile* &f) {
    if (f) { f->Close(); f = nullptr; }
}


// =========================================================
// ─── SECTION 5: FOLDER DISCOVERY ─────────────────────────
//
// Accepts BOTH old format ("0p1", "5") and new format
// ("p0p1", "m0p2") for the parameter sub-folder names.
// =========================================================
Bool_t IsParamFolder(const TString &name) {
    if (name.Length() == 0) return false;
    char first = name[0];
    // p/m prefix format:  "p0p1", "m0p2"
    if ((first == 'p' || first == 'm') && name.Length() > 1 && isdigit(name[1]))
        return true;
    // Literal-minus format: "-0p1", "-1", "-1p5"
    if (first == '-' && name.Length() > 1 && (isdigit(name[1]) || name[1] == '0'))
        return true;
    // Old/no-sign format:  "0p1", "0p0", "5", "1p5"
    if (isdigit(first))
        return true;
    return false;
}

// =========================================================
// ─── SECTION 6: ANALYSIS HELPER ──────────────────────────
//
// Runs CountEvents for one (angle, paramKey, mode) combination.
// mode: "elas" = elastic only     (full statistics)
//       "inel" = inelastic only   (full statistics)
//
// The two populations are always read in full; the physical
// inel/elas cross-section ratio is applied afterwards as a weight
// between the two independently-measured asymmetries (see the
// weighted combination in analyze_blpmc_scan()), not by
// subsampling events here.
//
// Returns a new PolAnalysis (caller owns it).
// =========================================================
struct RoleMap { TString eu, ed, iu, id; };   // paths to the four ROOT dirs

PolAnalysis* RunCounting(const RoleMap &rm,
                          Bool_t isP16,
                          const TString &mode,
                          Double_t beamPol,
                          Bool_t silent = false) {
    TString detName = isP16 ? "P16" : "P12";
    Int_t   nScint  = isP16 ? 4     : 3;

    PolAnalysis *ana = new PolAnalysis(detName, nScint);
    if (isP16) {
        ana->SetCuts(0, kP16_S0_min, kP16_S0_max);
        ana->SetCuts(1, kP16_S1_min, kP16_S1_max);
        ana->SetCuts(2, kP16_S2_min, kP16_S2_max);
        ana->SetCuts(3, kP16_S3_min, kP16_S3_max);
        ana->SetTheoreticalAy(kP16_TheoAy);
    } else {
        ana->SetCuts(0, kP12_S0_min, kP12_S0_max);
        ana->SetCuts(1, kP12_S1_min, kP12_S1_max);
        ana->SetCuts(2, kP12_S2_min, kP12_S2_max);
        ana->SetTheoreticalAy(kP12_TheoAy);
    }
    ana->SetBeamPolarization(beamPol);
    ana->SetBeamPolarizationError(0.0);

    Bool_t doElas = (mode == "elas");
    Bool_t doInel = (mode == "inel");

    // For silent mode, redirect cout so noisy OpenTree / CountEvents
    // messages are suppressed.
    std::streambuf *savedBuf = nullptr;
    std::ostringstream devnull;
    if (silent) {
        savedBuf = std::cout.rdbuf();
        std::cout.rdbuf(devnull.rdbuf());
    }

    TFile *hEU=nullptr, *hED=nullptr, *hIU=nullptr, *hID=nullptr;
    TTree *tEU=nullptr, *tED=nullptr, *tIU=nullptr, *tID=nullptr;

    if (doElas) {
        tEU = OpenTree(rm.eu, hEU);
        tED = OpenTree(rm.ed, hED);
        if (!tEU || !tED) {
            if (!silent) cout << "[Warning] Could not open elastic trees — skipping." << endl;
            doElas = false;
        }
    }
    if (doInel) {
        tIU = OpenTree(rm.iu, hIU);
        tID = OpenTree(rm.id, hID);
        if (!tIU || !tID) {
            if (!silent) cout << "[Warning] Could not open inelastic trees — skipping inel." << endl;
            doInel = false;
        }
    }

    if (doElas) {
        if (!silent) cout << "    Counting elastic SpinUp..."   << endl;
        ana->CountEvents(tEU, true);
        if (!silent) cout << "    Counting elastic SpinDown..." << endl;
        ana->CountEvents(tED, false);
    }
    if (doInel) {
        if (!silent) cout << Form("    Counting inel SpinUp   (full statistics, %lld events)...",
                                  tIU->GetEntries()) << endl;
        ana->CountEvents(tIU, true);
        if (!silent) cout << Form("    Counting inel SpinDown (full statistics, %lld events)...",
                                  tID->GetEntries()) << endl;
        ana->CountEvents(tID, false);
    }

    if (silent) std::cout.rdbuf(savedBuf);   // restore output

    SafeClose(hEU); SafeClose(hED);
    SafeClose(hIU); SafeClose(hID);
    return ana;
}

// =========================================================
// ComputeWeightedAy
//
// METHOD (revised): elastic and inelastic events are physically
// INDISTINGUISHABLE once they've survived the coincidence cuts — a real
// detector just sees "a hit passed the cut," with no per-event tag of
// which reaction produced it. So the combined asymmetry must be computed
// on POOLED counts, not as an average of two independently-computed
// cross-ratios (that approach — computing ε_elas and ε_inel separately
// and only combining them afterward — was found to disagree with pooling
// whenever ε_elas ≠ ε_inel, because the cross-ratio is a nonlinear
// function of the counts; see the derivation this replaces).
//
// The elastic and inelastic samples were generated with the SAME
// arbitrary statistics (~10M events each) regardless of their different
// physical cross sections, so before pooling, the inelastic counts must
// be rescaled by R0 = σ_inel/σ_elas (from GetInelFactor — exact, no
// associated uncertainty):
//
//   Lu = Lu_elas + R0·Lu_inel      (and similarly Ru, Ld, Rd)
//
//   A = √(Lu·Rd),  B = √(Ru·Ld),  S = A+B
//   ε_total = (A−B)/S
//   A_y     = ε_total / P_beam
//
// Error propagation: Lu, Ru, Ld, Rd are independent (built from disjoint
// elastic/inelastic count pairs), each with variance from Poisson
// statistics on the raw counts (R0 treated as an exact constant, so no
// R0² cross term beyond simple rescaling):
//
//   Var(Lu) = Lu_elas + R0²·Lu_inel     (and similarly Var(Ru), Var(Ld), Var(Rd))
//
//   δε_total = (A·B/S²) · √[ Var(Lu)/Lu² + Var(Ru)/Ru² + Var(Ld)/Ld² + Var(Rd)/Rd² ]
//
// This is the direct generalization of PolAnalysis::GetCrossRatioError():
// setting R0=1 (no elastic/inelastic split) reduces Var(Lu)=Lu exactly,
// and the formula above becomes algebraically identical to the original.
// =========================================================
struct WeightedAyResult {
    Double_t Ay, dAy;
    Double_t R0;                          // cross-section ratio used (exact, no error)
    Double_t Lu, Ru, Ld, Rd;              // pooled (R0-scaled) counts, for diagnostics
};

WeightedAyResult ComputeWeightedAy(PolAnalysis *anaE, PolAnalysis *anaI,
                                    Double_t R0, Double_t beamPol) {
    WeightedAyResult res = {0., 0., R0, 0., 0., 0., 0.};
    if (beamPol <= 0) return res;

    Double_t Lu_e = (Double_t)anaE->GetCountLeft(true),   Lu_i = (Double_t)anaI->GetCountLeft(true);
    Double_t Ru_e = (Double_t)anaE->GetCountRight(true),  Ru_i = (Double_t)anaI->GetCountRight(true);
    Double_t Ld_e = (Double_t)anaE->GetCountLeft(false),  Ld_i = (Double_t)anaI->GetCountLeft(false);
    Double_t Rd_e = (Double_t)anaE->GetCountRight(false), Rd_i = (Double_t)anaI->GetCountRight(false);

    Double_t Lu = Lu_e + R0 * Lu_i;
    Double_t Ru = Ru_e + R0 * Ru_i;
    Double_t Ld = Ld_e + R0 * Ld_i;
    Double_t Rd = Rd_e + R0 * Rd_i;
    res.Lu = Lu; res.Ru = Ru; res.Ld = Ld; res.Rd = Rd;

    // Need all four pooled totals strictly positive for a well-defined
    // cross-ratio — same convention as PolAnalysis::GetCrossRatioAsymmetry's
    // zero-count guard (returns Ay=0, dAy=0 rather than dividing by zero).
    if (Lu <= 0.0 || Ru <= 0.0 || Ld <= 0.0 || Rd <= 0.0) return res;

    Double_t A = TMath::Sqrt(Lu * Rd);
    Double_t B = TMath::Sqrt(Ru * Ld);
    Double_t S = A + B;

    Double_t epsTotal = (A - B) / S;

    Double_t varLu = Lu_e + R0*R0*Lu_i;
    Double_t varRu = Ru_e + R0*R0*Ru_i;
    Double_t varLd = Ld_e + R0*R0*Ld_i;
    Double_t varRd = Rd_e + R0*R0*Rd_i;

    Double_t AB_S2 = (A * B) / (S * S);
    Double_t depsTotal = AB_S2 * TMath::Sqrt(
        varLu/(Lu*Lu) + varRu/(Ru*Ru) + varLd/(Ld*Ld) + varRd/(Rd*Rd)
    );

    res.Ay  = epsTotal  / beamPol;
    res.dAy = depsTotal / beamPol;
    return res;
}


// =========================================================
// Helper — build a styled TGraphErrors
// =========================================================
TGraphErrors* MakeGraph(Int_t nPts, const TString &name, const TString &title,
                        Int_t markerStyle, Int_t color) {
    TGraphErrors *g = new TGraphErrors(nPts);
    g->SetName(name.Data());
    g->SetTitle(title.Data());
    g->SetMarkerStyle(markerStyle);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    return g;
}


// =========================================================
// ─── SECTION 7: MAIN SCAN FUNCTION ───────────────────────
//
// Usage:
//   analyze_blpmc_scan("ROTY_m0p2_p0p1_0p01_Pin_P80", 199.5)
//   analyze_blpmc_scan()                          // uses existing gScanDir, default energy
//   analyze_blpmc_scan("", 199.5, -0.1, 0.1)       // only process points in [-0.1, 0.1]
//
// startVal/stopVal restrict which parameter sub-folders are processed.
// Default (-1e30, +1e30) processes ALL discovered points (unchanged
// behaviour). Bounds are INCLUSIVE and compare against the parsed
// numeric parameter value (e.g. ROTY = -0.05), not the folder string.
//
// In addition to the A_y graphs (gGraph_P16/P12_elas/inel/weighted),
// this function ALSO computes and caches scintillator counts — both
// raw per-PMT occupancy (branch > 0) and full-coincidence-cut counts
// for elastic and full-statistics inelastic — so that draw_count_vs_param()
// can plot them WITHOUT re-opening any files. This adds extra
// RDataFrame count calls per parameter point (raw occupancy is new
// work; full-cut counts are already computed by RunCounting and were
// previously discarded).
// =========================================================
void analyze_blpmc_scan(TString userScanDir = "", Double_t beamEnergy = 200.0,
                         Double_t startVal = -1e30, Double_t stopVal = 1e30) {

    TString scanDir, kScanParam, kScanUnit;
    if (!ResolveScanDir(scanDir, kScanParam, kScanUnit, userScanDir)) return;

    cout << "\n=== Scan directory : " << gScanDir   << " ===" << endl;
    cout << "=== Parameter      : " << kScanParam  << " [" << kScanUnit << "] ===" << endl;
    cout << Form("=== Beam energy    : %.2f MeV (user-supplied) ===", beamEnergy) << endl;
    if (startVal > -1e29 || stopVal < 1e29)
        cout << Form("=== Range filter   : [%.6g, %.6g] ===\n", startVal, stopVal) << endl;
    else
        cout << "=== Range filter   : none (all points) ===\n" << endl;

    // ─── 1. Discover folder structure ────────────────────────────────
    // data[angle][paramKey][role] = folder path
    std::map<TString, std::map<TString, std::map<TString,TString>>> data;

    TSystemDirectory topDir(scanDir.Data(), scanDir.Data());
    TList *plutoFolders = topDir.GetListOfFiles();
    if (!plutoFolders) {
        cout << "[Error] Cannot list: " << scanDir << endl;
        return;
    }

    TSystemFile *plutoEntry;
    TIter nextPluto(plutoFolders);
    while ((plutoEntry = (TSystemFile*)nextPluto())) {
        TString pname = plutoEntry->GetName();
        if (!plutoEntry->IsDirectory()) continue;
        if (pname == "." || pname == "..") continue;

        Bool_t isUp   = pname.Contains("SpinUp");
        Bool_t isDown = pname.Contains("SpinDown");
        Bool_t isElas = pname.Contains("Elas") && !pname.Contains("Inel");
        Bool_t isInel = pname.Contains("Inel");
        Bool_t is12   = pname.Contains("_12p0_") || pname.EndsWith("_12p0");
        Bool_t is16   = pname.Contains("_16p2_") || pname.EndsWith("_16p2");

        if ((!isUp && !isDown) || (!isElas && !isInel) || (!is12 && !is16)) {
            cout << "[Skip] Unrecognised folder: " << pname << endl;
            continue;
        }

        TString angle = is12 ? "12p0" : "16p2";
        TString role;
        if      (isElas && isUp)   role = "elas_up";
        else if (isElas && isDown) role = "elas_down";
        else if (isInel && isUp)   role = "inel_up";
        else if (isInel && isDown) role = "inel_down";

        TString plutoPath = scanDir + "/" + pname;
        TSystemDirectory plutoDir(plutoPath.Data(), plutoPath.Data());
        TList *paramFolders = plutoDir.GetListOfFiles();
        if (!paramFolders) continue;

        TSystemFile *paramEntry;
        TIter nextParam(paramFolders);
        while ((paramEntry = (TSystemFile*)nextParam())) {
            TString paramName = paramEntry->GetName();
            if (!paramEntry->IsDirectory()) continue;
            if (paramName == "." || paramName == "..") continue;
            if (!IsParamFolder(paramName)) continue;
            data[angle][paramName][role] = plutoPath + "/" + paramName;
        }
    }

    if (data.empty()) {
        cout << "[Error] No valid folders found under " << scanDir << endl;
        return;
    }

    // ─── 2. Sort keys, then apply range filter ───────────────────────
    std::vector<TString> angles;
    for (auto &a : data) angles.push_back(a.first);
    sort(angles.begin(), angles.end());

    std::vector<TString> paramKeysAll;
    for (auto &t : data[angles[0]]) paramKeysAll.push_back(t.first);
    sort(paramKeysAll.begin(), paramKeysAll.end(), [](const TString &a, const TString &b){
        return ParsePValue(a) < ParsePValue(b);
    });

    std::vector<TString> paramKeys;
    for (auto &k : paramKeysAll) {
        Double_t v = ParsePValue(k);
        if (v >= startVal && v <= stopVal) paramKeys.push_back(k);
    }

    if (paramKeys.empty()) {
        cout << "[Error] No parameter points fall inside range [" << startVal
             << ", " << stopVal << "]" << endl;
        return;
    }

    cout << "Angles found:  ";
    for (auto &a : angles) cout << ParsePValue(a) << "  ";
    cout << endl;
    cout << kScanParam << " values processed (" << paramKeys.size() << "/" << paramKeysAll.size() << "): ";
    for (auto &t : paramKeys) cout << ParsePValue(t) << "  ";
    cout << endl;

    // ─── 3. Compute normalisation factors ────────────────────────────
    cout << "\n--- Computing inel/elas normalization factors ---" << endl;
    // Reset cached values so they are recomputed for this energy
    gInelFactor_P16 = -1.0;
    gInelFactor_P12 = -1.0;
    gBeamPol        = -1.0;   // re-parse beam polarization for this scan too
    Double_t inelFactor_P16 = GetInelFactor(true,  beamEnergy, gScanDir);
    Double_t inelFactor_P12 = GetInelFactor(false, beamEnergy, gScanDir);
    Double_t beamPol        = GetBeamPolarization(gScanDir);
    cout << Form("  P16: R = %.6f", inelFactor_P16) << endl;
    cout << Form("  P12: R = %.6f", inelFactor_P12) << endl;
    cout << Form("  Beam polarization P = %.4f\n", beamPol) << endl;

    // ─── 4. Allocate result graphs ────────────────────────────────────
    Int_t nPoints = (Int_t)paramKeys.size();

    auto makeTitle = [&](const TString &det, const TString &mode) -> TString {
        TString angle = det.Contains("16") ? "16.2" : "12.0";
        return Form("A_{y} vs %s  (%s, %s#circ, %s);%s [%s];A_{y}",
                    kScanParam.Data(), det.Data(), angle.Data(),
                    mode.Data(), kScanParam.Data(), kScanUnit.Data());
    };

    // Delete old A_y graphs if they exist
    auto safeDelG = [](TGraphErrors* &g) { if (g) { delete g; g = nullptr; } };
    safeDelG(gGraph_P16_elas);     safeDelG(gGraph_P16_inel);     safeDelG(gGraph_P16_weighted);
    safeDelG(gGraph_P12_elas);     safeDelG(gGraph_P12_inel);     safeDelG(gGraph_P12_weighted);

    gGraph_P16_elas     = MakeGraph(nPoints, "gAy_P16_elas",     makeTitle("P16","elastic"),          21, kBlue+1);
    gGraph_P16_inel     = MakeGraph(nPoints, "gAy_P16_inel",     makeTitle("P16","inelastic"),         24, kBlue+1);
    gGraph_P16_weighted = MakeGraph(nPoints, "gAy_P16_weighted", makeTitle("P16","weighted"),          33, kViolet+1);

    gGraph_P12_elas     = MakeGraph(nPoints, "gAy_P12_elas",     makeTitle("P12","elastic"),          21, kRed+1);
    gGraph_P12_inel     = MakeGraph(nPoints, "gAy_P12_inel",     makeTitle("P12","inelastic"),         24, kRed+1);
    gGraph_P12_weighted = MakeGraph(nPoints, "gAy_P12_weighted", makeTitle("P12","weighted"),          33, kViolet+1);

    // Delete old count-cache graphs and allocate fresh ones.
    // Layout: index = is*16 + slot  (slot 0-15, see cache comment above globals)
    auto allocCountGraphs = [&](Int_t nScint) -> std::vector<TGraphErrors*> {
        std::vector<TGraphErrors*> v(nScint * 16, nullptr);
        for (auto &g : v) g = new TGraphErrors(nPoints);
        return v;
    };
    auto freeCountGraphs = [](std::vector<TGraphErrors*> &v) {
        for (auto &g : v) if (g) { delete g; g = nullptr; }
        v.clear();
    };
    freeCountGraphs(gCountGraphs_P16);
    freeCountGraphs(gCountGraphs_P12);
    gCountGraphs_P16 = allocCountGraphs(4);   // P16: 4 scintillators
    gCountGraphs_P12 = allocCountGraphs(3);   // P12: 3 scintillators
    gCountsValid_P16 = false;
    gCountsValid_P12 = false;

    // Raw single-scintillator threshold count (branch > 0), via RDataFrame
    // so it benefits from the same implicit-MT parallelization as
    // PolAnalysis::CountEvents.
    auto countBranchRaw = [](TTree *t, const char *branch) -> Long64_t {
        if (!t) return 0LL;
        ROOT::RDataFrame df(*t);
        return *df.Filter(Form("%s>0", branch)).Count();
    };

    std::vector<TString> branchL_P16 = {"eDepP16LS0","eDepP16LS1","eDepP16LS2","eDepP16LS3"};
    std::vector<TString> branchR_P16 = {"eDepP16RS0","eDepP16RS1","eDepP16RS2","eDepP16RS3"};
    std::vector<TString> branchL_P12 = {"eDepP12LS0","eDepP12LS1","eDepP12LS2"};
    std::vector<TString> branchR_P12 = {"eDepP12RS0","eDepP12RS1","eDepP12RS2"};

    // ─── 5. Loop over parameter values ───────────────────────────────
    for (Int_t ipt = 0; ipt < nPoints; ipt++) {
        TString  paramKey = paramKeys[ipt];
        Double_t paramVal = ParsePValue(paramKey);

        cout << "\n=============================" << endl;
        cout << " " << kScanParam << " = " << paramVal << endl;
        cout << "=============================" << endl;

        for (auto &angle : angles) {
            Bool_t isP16      = angle.Contains("16");
            Double_t inelFac  = isP16 ? inelFactor_P16 : inelFactor_P12;

            cout << "\n  --- Angle: " << ParsePValue(angle) << " deg ---" << endl;

            auto &roleMap = data[angle][paramKey];
            Bool_t ok = (roleMap.count("elas_up")  && roleMap.count("elas_down") &&
                         roleMap.count("inel_up")   && roleMap.count("inel_down"));
            if (!ok) {
                cout << "[Warning] Missing folders for angle=" << angle
                     << " " << kScanParam << "=" << paramKey << " — skipping." << endl;
                continue;
            }

            TString fEU = FindOutputRoot(roleMap["elas_up"]);
            TString fED = FindOutputRoot(roleMap["elas_down"]);
            TString fIU = FindOutputRoot(roleMap["inel_up"]);
            TString fID = FindOutputRoot(roleMap["inel_down"]);

            if (fEU.IsNull()||fED.IsNull()||fIU.IsNull()||fID.IsNull()) {
                cout << "[Warning] Missing ROOT file for angle=" << angle
                     << " " << kScanParam << "=" << paramKey << " — skipping." << endl;
                continue;
            }

            RoleMap rm = { fEU, fED, fIU, fID };

            // Two counting passes, each reading its sample in FULL —
            // no event subsampling. Files are opened/closed inside
            // RunCounting; the RDataFrame-based CountEvents call is
            // multithreaded automatically (see rootlogon.C).
            cout << "  [elastic]" << endl;
            PolAnalysis *anaE = RunCounting(rm, isP16, "elas", beamPol, true);

            cout << "  [inelastic]" << endl;
            PolAnalysis *anaI = RunCounting(rm, isP16, "inel", beamPol, true);

            anaE->PrintResults();
            anaI->PrintResults();

            // Store the independent A_y measurements
            TGraphErrors *gE = isP16 ? gGraph_P16_elas : gGraph_P12_elas;
            TGraphErrors *gI = isP16 ? gGraph_P16_inel : gGraph_P12_inel;

            gE->SetPoint(ipt, paramVal, anaE->GetAnalyzingPower_CrossRatio());
            gE->SetPointError(ipt, 0., anaE->GetAnalyzingPowerError_CrossRatio());
            gI->SetPoint(ipt, paramVal, anaI->GetAnalyzingPower_CrossRatio());
            gI->SetPointError(ipt, 0., anaI->GetAnalyzingPowerError_CrossRatio());

            // ── Cross-section-weighted combination (pooled counts) ──────
            // See ComputeWeightedAy() for the full derivation. In short:
            // elastic and inelastic events are indistinguishable once
            // they've survived the coincidence cuts, so the combined
            // asymmetry is computed on POOLED counts (inelastic rescaled
            // by R0 = inelFac to correct for the equal, arbitrary
            // generated statistics), not by averaging two independently
            // computed cross-ratios.
            {
                TGraphErrors *gW = isP16 ? gGraph_P16_weighted : gGraph_P12_weighted;
                WeightedAyResult wr = ComputeWeightedAy(anaE, anaI, inelFac, beamPol);

                gW->SetPoint(ipt, paramVal, wr.Ay);
                gW->SetPointError(ipt, 0., wr.dAy);
                cout << Form("  [weighted] A_y = %.6f ± %.6f  (R0=%.6f, pooled L↑=%.0f R↑=%.0f L↓=%.0f R↓=%.0f)",
                             wr.Ay, wr.dAy, wr.R0, wr.Lu, wr.Ru, wr.Ld, wr.Rd) << endl;
            }

            // ── Cache full-cut counts (already computed by RunCounting,
            //    just extract them instead of discarding) ──────────────
            std::vector<TGraphErrors*> &cg = isP16 ? gCountGraphs_P16 : gCountGraphs_P12;
            Int_t nScint = isP16 ? 4 : 3;

            Long64_t eLU = anaE->GetCountLeft(true),  eLD = anaE->GetCountLeft(false);
            Long64_t eRU = anaE->GetCountRight(true), eRD = anaE->GetCountRight(false);
            Long64_t iLU = anaI->GetCountLeft(true),  iLD = anaI->GetCountLeft(false);
            Long64_t iRU = anaI->GetCountRight(true), iRD = anaI->GetCountRight(false);

            for (Int_t is = 0; is < nScint; is++) {
                // slots 4-7: full-cut elastic (same value repeated across all is)
                cg[is*16+4]->SetPoint(ipt, paramVal, (Double_t)eLU);
                cg[is*16+4]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)eLU));
                cg[is*16+5]->SetPoint(ipt, paramVal, (Double_t)eLD);
                cg[is*16+5]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)eLD));
                cg[is*16+6]->SetPoint(ipt, paramVal, (Double_t)eRU);
                cg[is*16+6]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)eRU));
                cg[is*16+7]->SetPoint(ipt, paramVal, (Double_t)eRD);
                cg[is*16+7]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)eRD));

                // slots 12-15: full-cut inelastic (full statistics)
                cg[is*16+12]->SetPoint(ipt, paramVal, (Double_t)iLU);
                cg[is*16+12]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)iLU));
                cg[is*16+13]->SetPoint(ipt, paramVal, (Double_t)iLD);
                cg[is*16+13]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)iLD));
                cg[is*16+14]->SetPoint(ipt, paramVal, (Double_t)iRU);
                cg[is*16+14]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)iRU));
                cg[is*16+15]->SetPoint(ipt, paramVal, (Double_t)iRD);
                cg[is*16+15]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)iRD));
            }

            // ── Raw per-scintillator occupancy (slots 0-3 elastic, 8-11
            //    inelastic) — NEW work, requires re-opening all four files
            //    since RunCounting already closed them. Computed for BOTH
            //    reactions unconditionally, so draw_count_vs_param() can
            //    show a true raw-vs-cut comparison for either one. ───────
            std::vector<TString> &bL = isP16 ? branchL_P16 : branchL_P12;
            std::vector<TString> &bR = isP16 ? branchR_P16 : branchR_P12;

            auto fillRaw = [&](const TString &fileUp, const TString &fileDown, Int_t slotBase) {
                TFile *hU=nullptr, *hD=nullptr;
                TTree *tU = OpenTree(fileUp,   hU);
                TTree *tD = OpenTree(fileDown, hD);
                if (tU && tD) {
                    for (Int_t is = 0; is < nScint; is++) {
                        Long64_t luR = countBranchRaw(tU, bL[is].Data());
                        Long64_t ldR = countBranchRaw(tD, bL[is].Data());
                        Long64_t ruR = countBranchRaw(tU, bR[is].Data());
                        Long64_t rdR = countBranchRaw(tD, bR[is].Data());

                        cg[is*16+slotBase+0]->SetPoint(ipt, paramVal, (Double_t)luR);
                        cg[is*16+slotBase+0]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)luR));
                        cg[is*16+slotBase+1]->SetPoint(ipt, paramVal, (Double_t)ldR);
                        cg[is*16+slotBase+1]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)ldR));
                        cg[is*16+slotBase+2]->SetPoint(ipt, paramVal, (Double_t)ruR);
                        cg[is*16+slotBase+2]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)ruR));
                        cg[is*16+slotBase+3]->SetPoint(ipt, paramVal, (Double_t)rdR);
                        cg[is*16+slotBase+3]->SetPointError(ipt, 0., TMath::Sqrt((Double_t)rdR));
                    }
                } else {
                    cout << "[Warning] Could not reopen files for raw counts (slot base "
                         << slotBase << ")." << endl;
                }
                SafeClose(hU);
                SafeClose(hD);
            };

            fillRaw(fEU, fED, 0);    // elastic raw   -> slots 0-3
            fillRaw(fIU, fID, 8);    // inelastic raw -> slots 8-11

            if (isP16) gCountsValid_P16 = true; else gCountsValid_P12 = true;

            delete anaE; delete anaI;
        }
    }

    cout << "\n=== Scan complete ===" << endl;
    cout << "Call draw_scan_p16() / draw_scan_p12() / draw_scan_both() to view results." << endl;
    cout << "  mode = \"weighted\" (default, correct)  \"elas\"  \"inel\"  \"all\"" << endl;
    cout << "Call draw_count_vs_param() to view cached scintillator counts (no re-reading)." << endl;
    cout << "Call SaveScanCache() to persist this cache to disk (Cache/<gScanDir>.root)," << endl;
    cout << "  so a future ROOT session can LoadScanCache() and skip re-scanning entirely." << endl;
}


// =========================================================
// ─── SECTION 8: DRAW FUNCTIONS ───────────────────────────
//
// draw_scan_p16()              → weighted (default, correct)
// draw_scan_p16("elas")        → elastic-only
// draw_scan_p16("inel")        → inelastic-only
// draw_scan_p16("all")         → all three on the same canvas
//
// All draw functions use the CACHED graphs — no file re-reading.
// =========================================================

// Internal: compute a padded Y range covering all points (± error bars)
// across one or more graphs, skipping any non-finite (NaN) points — the
// true-singularity guard in PolAnalysis/ComputeWeightedAy can now produce
// occasional NaN points, and ROOT's own TGraph auto-range does not reliably
// skip them, which can otherwise corrupt the whole axis range.
void ComputeYRangeFinite(const std::vector<TGraphErrors*> &graphs,
                          Double_t &ylo, Double_t &yhi) {
    ylo = std::numeric_limits<Double_t>::max();
    yhi = -std::numeric_limits<Double_t>::max();
    for (auto *g : graphs) {
        if (!g) continue;
        Int_t n = g->GetN();
        Double_t *y  = g->GetY();
        Double_t *ey = g->GetEY();
        for (Int_t i = 0; i < n; i++) {
            if (!TMath::Finite(y[i])) continue;
            Double_t e = (ey && TMath::Finite(ey[i])) ? ey[i] : 0.0;
            Double_t lo = y[i] - e;
            Double_t hi = y[i] + e;
            if (lo < ylo) ylo = lo;
            if (hi > yhi) yhi = hi;
        }
    }
    if (ylo > yhi) { ylo = 0.0; yhi = 1.0; }   // fallback: every point was non-finite
    Double_t pad = 0.08 * (yhi - ylo);
    if (pad <= 0.0) pad = 0.1;
    ylo -= pad;
    yhi += pad;
}

void DrawScanGraph(TGraphErrors *g, const TString &param, const TString &unit,
                   Double_t theoAy, const TString &canvName, const TString &canvTitle) {
    TCanvas *c = new TCanvas(canvName.Data(), canvTitle.Data(), 900, 600);
    c->SetGrid();
    g->Draw("AP");
    g->GetXaxis()->SetTitle(Form("%s [%s]", param.Data(), unit.Data()));
    g->GetYaxis()->SetTitle("A_{y}");

    // Explicit NaN-safe range — don't rely on ROOT's own auto-range, which
    // does not reliably exclude non-finite points (the true-singularity
    // guard in PolAnalysis/ComputeWeightedAy can now produce occasional
    // NaN points).
    Double_t ylo, yhi;
    ComputeYRangeFinite({g}, ylo, yhi);
    Double_t ymin = TMath::Min(ylo, theoAy);
    Double_t ymax = TMath::Max(yhi, theoAy);
    g->GetYaxis()->SetRangeUser(ymin, ymax);

    // Draw theory line spanning graph x-range
    Double_t xlo = g->GetXaxis()->GetXmin();
    Double_t xhi = g->GetXaxis()->GetXmax();
    TLine *lTheo = new TLine(xlo, theoAy, xhi, theoAy);
    lTheo->SetLineStyle(2);
    lTheo->SetLineColor(kGray+2);
    lTheo->SetLineWidth(2);
    lTheo->Draw();

    TLegend *leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(g,     g->GetTitle(),                             "lep");
    leg->AddEntry(lTheo, Form("Theory A_{y} = %.3f", theoAy), "l");
    leg->Draw();
    c->Update();
}

// Internal: overlay elas / inel / weighted on one canvas
void DrawScanAll(TGraphErrors *gE, TGraphErrors *gI, TGraphErrors *gW,
                 const TString &param, const TString &unit,
                 Double_t theoAy, const TString &det) {

    TCanvas *c = new TCanvas(Form("cScan_%s_all", det.Data()),
                             Form("A_{y} vs %s — %s (all modes)", param.Data(), det.Data()),
                             900, 600);
    c->SetGrid();

    Double_t xlo = gW->GetXaxis()->GetXmin();
    Double_t xhi = gW->GetXaxis()->GetXmax();

    gE->SetMarkerStyle(21); gE->SetMarkerColor(kBlue+1);   gE->SetLineColor(kBlue+1);   gE->SetLineStyle(1);
    gI->SetMarkerStyle(24); gI->SetMarkerColor(kOrange+7); gI->SetLineColor(kOrange+7); gI->SetLineStyle(2);
    gW->SetMarkerStyle(33); gW->SetMarkerColor(kViolet+1); gW->SetLineColor(kViolet+1); gW->SetLineStyle(1);

    gE->Draw("AP");   // draw first to set axes
    gE->GetXaxis()->SetTitle(Form("%s [%s]", param.Data(), unit.Data()));
    gE->GetYaxis()->SetTitle("A_{y}");
    gI->Draw("P SAME");
    gW->Draw("P SAME");

    // Fix: the axis range set by gE->Draw("AP") only reflects gE's own
    // data — gI/gW were being drawn "P SAME" onto those pre-existing axes
    // with no range adjustment, so any inelastic point (or error bar)
    // outside elastic's range was silently clipped off the visible plot.
    // Recompute the range across all three graphs (± error bars, skipping
    // any NaN points) and apply it explicitly.
    Double_t ylo, yhi;
    ComputeYRangeFinite({gE, gI, gW}, ylo, yhi);
    Double_t ymin = TMath::Min(ylo, theoAy);
    Double_t ymax = TMath::Max(yhi, theoAy);
    gE->GetYaxis()->SetRangeUser(ymin, ymax);

    TLine *lTheo = new TLine(xlo, theoAy, xhi, theoAy);
    lTheo->SetLineStyle(2);
    lTheo->SetLineColor(kGray+2);
    lTheo->SetLineWidth(2);
    lTheo->Draw();

    TLegend *leg = new TLegend(0.55, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gE,    "Elastic only",                    "lep");
    leg->AddEntry(gI,    "Inelastic only",                  "lep");
    leg->AddEntry(gW,    "Weighted (cross-section-ratio)",  "lep");
    leg->AddEntry(lTheo, Form("Theory = %.3f", theoAy),     "l");
    leg->Draw();
    c->Update();
}

void draw_scan_p16(TString mode = "weighted") {
    if (!gGraph_P16_weighted) {
        cout << "[Error] No scan data. Run analyze_blpmc_scan() first!" << endl;
        return;
    }
    TString param, unit;
    ParseScanDir(gScanDir, param, unit);

    if (mode == "all") {
        DrawScanAll(gGraph_P16_elas, gGraph_P16_inel, gGraph_P16_weighted,
                    param, unit, kP16_TheoAy, "P16");
        return;
    }
    TGraphErrors *g = nullptr;
    if      (mode == "elas") g = gGraph_P16_elas;
    else if (mode == "inel") g = gGraph_P16_inel;
    else                     g = gGraph_P16_weighted;   // default

    DrawScanGraph(g, param, unit, kP16_TheoAy,
                  Form("cScan_P16_%s", mode.Data()),
                  Form("A_{y} vs %s — P16 16.2° (%s)", param.Data(), mode.Data()));
}

void draw_scan_p12(TString mode = "weighted") {
    if (!gGraph_P12_weighted) {
        cout << "[Error] No scan data. Run analyze_blpmc_scan() first!" << endl;
        return;
    }
    TString param, unit;
    ParseScanDir(gScanDir, param, unit);

    if (mode == "all") {
        DrawScanAll(gGraph_P12_elas, gGraph_P12_inel, gGraph_P12_weighted,
                    param, unit, kP12_TheoAy, "P12");
        return;
    }
    TGraphErrors *g = nullptr;
    if      (mode == "elas") g = gGraph_P12_elas;
    else if (mode == "inel") g = gGraph_P12_inel;
    else                     g = gGraph_P12_weighted;   // default

    DrawScanGraph(g, param, unit, kP12_TheoAy,
                  Form("cScan_P12_%s", mode.Data()),
                  Form("A_{y} vs %s — P12 12.0° (%s)", param.Data(), mode.Data()));
}

// Convenience: draw all modes for both detectors
void draw_scan_both(TString mode = "weighted") {
    draw_scan_p16(mode);
    draw_scan_p12(mode);
}

// Convenience: show all three modes for both detectors in one call
void draw_scan_all_modes() {
    draw_scan_p16("all");
    draw_scan_p12("all");
}


// =========================================================
// ─── SECTION 9: INSPECTION MODE ──────────────────────────
//
// Usage:
//   inspect("12p0", -0.1)      // uses gScanDir, both spins
//   inspect("16p2",  0.0, 1)   // SpinUp only
//   inspect("16p2", -0.2, -1)  // SpinDown only
//
// After calling inspect(), use:
//   insp_draw_1d_raw()
//   insp_draw_1d_comparison()
//   insp_draw_1d_spin()
//   insp_draw_2d()
//   insp_draw_all()
// =========================================================

// =========================================================
// InspFillTree
//
// Fills the inspection histograms from EVERY entry in the tree —
// elastic and inelastic alike. There is no subsampling; the
// cross-section-ratio weighting between elastic and inelastic is
// applied only in the summary A_y numbers (see inspect()), not by
// dropping events at fill time.
// =========================================================
void InspFillTree(TTree *tree, PolHistograms *hist, PolAnalysis *ana,
                  Int_t nScint, Bool_t isSpinUp) {
    if (!tree) return;

    Double_t eDepP16LS0, eDepP16LS1, eDepP16LS2, eDepP16LS3;
    Double_t eDepP16RS0, eDepP16RS1, eDepP16RS2, eDepP16RS3;
    Double_t eDepP12LS0, eDepP12LS1, eDepP12LS2;
    Double_t eDepP12RS0, eDepP12RS1, eDepP12RS2;

    tree->SetBranchAddress("eDepP16LS0", &eDepP16LS0);
    tree->SetBranchAddress("eDepP16LS1", &eDepP16LS1);
    tree->SetBranchAddress("eDepP16LS2", &eDepP16LS2);
    tree->SetBranchAddress("eDepP16LS3", &eDepP16LS3);
    tree->SetBranchAddress("eDepP16RS0", &eDepP16RS0);
    tree->SetBranchAddress("eDepP16RS1", &eDepP16RS1);
    tree->SetBranchAddress("eDepP16RS2", &eDepP16RS2);
    tree->SetBranchAddress("eDepP16RS3", &eDepP16RS3);
    tree->SetBranchAddress("eDepP12LS0", &eDepP12LS0);
    tree->SetBranchAddress("eDepP12LS1", &eDepP12LS1);
    tree->SetBranchAddress("eDepP12LS2", &eDepP12LS2);
    tree->SetBranchAddress("eDepP12RS0", &eDepP12RS0);
    tree->SetBranchAddress("eDepP12RS1", &eDepP12RS1);
    tree->SetBranchAddress("eDepP12RS2", &eDepP12RS2);

    Long64_t nTotal = tree->GetEntries();

    for (Long64_t i = 0; i < nTotal; i++) {
        tree->GetEntry(i);
        if (nScint == 4) {
            Double_t e[8] = {eDepP16LS0, eDepP16LS1, eDepP16LS2, eDepP16LS3,
                             eDepP16RS0, eDepP16RS1, eDepP16RS2, eDepP16RS3};
            for (Int_t j = 0; j < 8; j++) {
                hist->Fill(j, e);
                hist->FillWithCuts(j, e, ana);
                hist->FillWithCutsSpin(j, e, ana, isSpinUp);
            }
            hist->Fill2D(e);
        } else {
            Double_t e[6] = {eDepP12LS0, eDepP12LS1, eDepP12LS2,
                             eDepP12RS0, eDepP12RS1, eDepP12RS2};
            for (Int_t j = 0; j < 6; j++) {
                hist->Fill(j, e);
                hist->FillWithCuts(j, e, ana);
                hist->FillWithCutsSpin(j, e, ana, isSpinUp);
            }
            hist->Fill2D(e);
        }
    }
}

// =========================================================
// Find a reaction/spin/angle sub-folder by scanning
// (works for both MT and Pin naming conventions)
// =========================================================
TString FindSubDir(const TString &dataPath,
                   const TString &reac,
                   const TString &spinStr,
                   const TString &angle) {
    TSystemDirectory dir(dataPath.Data(), dataPath.Data());
    TList *entries = dir.GetListOfFiles();
    if (!entries) return "";
    TSystemFile *e;
    TIter next(entries);
    while ((e = (TSystemFile*)next())) {
        TString n = e->GetName();
        if (!e->IsDirectory() || n == "." || n == "..") continue;
        Bool_t hasAngle = n.Contains("_" + angle + "_") || n.EndsWith("_" + angle);
        if (n.Contains(reac) && n.Contains(spinStr) && hasAngle)
            return dataPath + "/" + n;
    }
    return "";
}

// =========================================================
// Resolve paramVal → folder name, trying multiple formats
// against what actually exists on disk:
//   1. New format:           "p0p1", "m0p2", "p0"  (sign + decimal-as-p)
//   2. Old format, decimal:  "0p0", "0p1", "5p0"    (decimal-as-p, no sign char)
//   3. Old format, integer:  "0", "5"                (bare integer, legacy)
//
// IMPORTANT: candidate 2 is generated unconditionally — even for whole
// numbers — because some scans (e.g. ones using ParsePValue's old-format
// branch) write "0p0" for zero, not bare "0". Relying on integer
// rounding to decide the format silently fails for exactly this case.
// =========================================================
TString ResolveParamKey(Double_t paramVal, const TString &baseDir) {
    std::vector<TString> candidates;

    // 1. p/m prefix format (new):       "p0p1", "m0p2", "p0"
    candidates.push_back(FormatPValue(paramVal));

    // 2. Literal-minus format (newest dataset): "-0p1", "-1p5", "-1", "0p1", "1"
    {
        TString litFmt = Form("%.10g", TMath::Abs(paramVal));
        if (!litFmt.Contains(".")) litFmt += ".0";
        litFmt.ReplaceAll(".", "p");
        if (paramVal < 0.0)
            candidates.push_back("-" + litFmt);  // "-0p1", "-1p5"
        else
            candidates.push_back(litFmt);         // "0p1", "1p5" (no sign)
    }

    // 3. Old decimal-as-p, unconditionally (covers "0p0", "0p1", "5p0")
    {
        TString oldFmtDec = Form("%.10g", TMath::Abs(paramVal));
        if (!oldFmtDec.Contains("."))  oldFmtDec += ".0";
        oldFmtDec.ReplaceAll(".", "p");
        candidates.push_back(oldFmtDec);          // "0p0", "0p1"
    }

    // 4. Old format, bare integer (legacy fallback: "0", "5")
    Double_t roundedVal = TMath::Nint(paramVal);
    if (TMath::Abs(paramVal - roundedVal) < 1e-9)
        candidates.push_back(Form("%lld", (Long64_t)roundedVal));

    // Deduplicate and try each candidate
    std::vector<TString> seen;
    for (auto &cand : candidates) {
        Bool_t dup = false;
        for (auto &s : seen) if (s == cand) { dup = true; break; }
        if (dup) continue;
        seen.push_back(cand);
        TString path = baseDir + "/" + cand;
        if (gSystem->AccessPathName(path.Data()) == kFALSE)
            return cand;
    }

    // None found — report all candidates tried
    cout << "[Warning] No folder matched paramVal=" << paramVal
         << " under " << baseDir << endl;
    cout << "          Tried: ";
    for (auto &cand : seen) cout << "\"" << cand << "\" ";
    cout << endl;
    return FormatPValue(paramVal);
}

void inspect(TString angleKey, Double_t paramVal, Int_t spin = 0,
             Double_t beamEnergy = 200.0) {
    if (gScanDir == "") {
        cout << "[Error] No scan directory set. Set: gScanDir = \"...\"" << endl;
        return;
    }

    TString kScanParam, kScanUnit;
    ParseScanDir(gScanDir, kScanParam, kScanUnit);

    Bool_t  isP16     = angleKey.Contains("16");
    TString detName   = isP16 ? "P16" : "P12";
    Int_t   nScint    = isP16 ? 4 : 3;
    TString angleFull = isP16 ? "16p2" : "12p0";

    TString dataPath = "Data/" + gScanDir;

    // Resolve parameter folder name (handles old and new conventions)
    TString subEU = FindSubDir(dataPath, "Elas",    "SpinUp",   angleFull);
    TString subED = FindSubDir(dataPath, "Elas",    "SpinDown", angleFull);
    TString subIU = FindSubDir(dataPath, "Inel443", "SpinUp",   angleFull);
    TString subID = FindSubDir(dataPath, "Inel443", "SpinDown", angleFull);

    if (subEU.IsNull()||subED.IsNull()||subIU.IsNull()||subID.IsNull()) {
        cout << "[Error] Could not find one or more subdirectories in: " << dataPath << endl;
        cout << "  EU: " << (subEU.IsNull() ? "NOT FOUND" : subEU) << endl;
        cout << "  ED: " << (subED.IsNull() ? "NOT FOUND" : subED) << endl;
        cout << "  IU: " << (subIU.IsNull() ? "NOT FOUND" : subIU) << endl;
        cout << "  ID: " << (subID.IsNull() ? "NOT FOUND" : subID) << endl;
        return;
    }

    // Resolve the parameter sub-folder name (handles new and old conventions)
    TString paramKey = ResolveParamKey(paramVal, subEU);
    cout << Form("  Resolved paramVal=%.6g → folder name \"%s\"", paramVal, paramKey.Data()) << endl;

    TString dirEU = subEU + "/" + paramKey;
    TString dirED = subED + "/" + paramKey;
    TString dirIU = subIU + "/" + paramKey;
    TString dirID = subID + "/" + paramKey;

    TString fEU = FindOutputRoot(dirEU);
    TString fED = FindOutputRoot(dirED);
    TString fIU = FindOutputRoot(dirIU);
    TString fID = FindOutputRoot(dirID);

    if (fEU.IsNull()||fED.IsNull()||fIU.IsNull()||fID.IsNull()) {
        cout << "[Error] Could not find ROOT files:" << endl;
        cout << "  " << dirEU << "  → " << (fEU.IsNull() ? "NOT FOUND" : fEU) << endl;
        cout << "  " << dirED << "  → " << (fED.IsNull() ? "NOT FOUND" : fED) << endl;
        cout << "  " << dirIU << "  → " << (fIU.IsNull() ? "NOT FOUND" : fIU) << endl;
        cout << "  " << dirID << "  → " << (fID.IsNull() ? "NOT FOUND" : fID) << endl;
        return;
    }

    TString spinTag  = (spin == 1) ? "_SpinUp" : (spin == -1) ? "_SpinDown" : "_AllSpin";
    gInsp_label = Form("%s_%s%g%s", detName.Data(), kScanParam.Data(), paramVal, spinTag.Data());

    cout << "\n=== Inspection: " << detName
         << "  " << kScanParam << "=" << paramVal
         << "  [" << gScanDir << "] ===" << endl;

    // Compute normalization for this detector using user-supplied beam energy
    cout << Form("  Beam energy    = %.2f MeV (user-supplied)", beamEnergy) << endl;
    Double_t inelFac = GetInelFactor(isP16, beamEnergy, gScanDir);
    Double_t beamPol = GetBeamPolarization(gScanDir);
    gInsp_R = inelFac;

    TFile *hEU=nullptr, *hED=nullptr, *hIU=nullptr, *hID=nullptr;
    TTree *tEU = OpenTree(fEU, hEU);
    TTree *tED = OpenTree(fED, hED);
    TTree *tIU = OpenTree(fIU, hIU);
    TTree *tID = OpenTree(fID, hID);

    if (!tEU||!tED||!tIU||!tID) {
        cout << "[Error] Failed to open one or more trees." << endl;
        SafeClose(hEU); SafeClose(hED);
        SafeClose(hIU); SafeClose(hID);
        return;
    }

    cout << Form("  Inel factor R = %.6f (cross-section ratio; applied only when"
                 " combining A_y, not to event counts)", inelFac) << endl;
    cout << "  Inel SpinUp:   " << tIU->GetEntries() << " total (read in full)" << endl;
    cout << "  Inel SpinDown: " << tID->GetEntries() << " total (read in full)" << endl;

    // Rebuild inspection objects
    if (gInsp_hist) { delete gInsp_hist; gInsp_hist = nullptr; }
    if (gInsp_ana)  { delete gInsp_ana;  gInsp_ana  = nullptr; }
    if (gInsp_ana_inel) { delete gInsp_ana_inel; gInsp_ana_inel = nullptr; }
    if (gInsp_hist_elas) { delete gInsp_hist_elas; gInsp_hist_elas = nullptr; }
    if (gInsp_hist_inel) { delete gInsp_hist_inel; gInsp_hist_inel = nullptr; }

    TH1::AddDirectory(kFALSE);
    gInsp_hist = new PolHistograms(detName, detName, nScint, 500, 0., 50.);
    gInsp_hist->Create2DHistograms();
    // Separate elastic-only and inelastic-only sets. Create2DHistograms()
    // is required here even though insp_draw_elas_inel only uses 1D raw
    // histograms — InspFillTree unconditionally calls Fill2D(), which
    // indexes into h2D[] and would be out-of-bounds on an empty vector.
    gInsp_hist_elas = new PolHistograms(detName + "_elas", detName + " Elastic", nScint, 500, 0., 50.);
    gInsp_hist_elas->Create2DHistograms();
    gInsp_hist_inel = new PolHistograms(detName + "_inel", detName + " Inelastic", nScint, 500, 0., 50.);
    gInsp_hist_inel->Create2DHistograms();
    TH1::AddDirectory(kTRUE);

    auto configureCuts = [&](PolAnalysis *ana) {
        if (isP16) {
            ana->SetCuts(0, kP16_S0_min, kP16_S0_max);
            ana->SetCuts(1, kP16_S1_min, kP16_S1_max);
            ana->SetCuts(2, kP16_S2_min, kP16_S2_max);
            ana->SetCuts(3, kP16_S3_min, kP16_S3_max);
            ana->SetTheoreticalAy(kP16_TheoAy);
        } else {
            ana->SetCuts(0, kP12_S0_min, kP12_S0_max);
            ana->SetCuts(1, kP12_S1_min, kP12_S1_max);
            ana->SetCuts(2, kP12_S2_min, kP12_S2_max);
            ana->SetTheoreticalAy(kP12_TheoAy);
        }
        ana->SetBeamPolarization(beamPol);
        ana->SetBeamPolarizationError(0.0);
    };

    // gInsp_ana tracks ELASTIC counts (and is also read by the draw
    // functions for cut boundaries). gInsp_ana_inel tracks INELASTIC
    // counts separately — pooling the two into one object would repeat
    // the same "combined(pooled)" bias the scan moved away from.
    gInsp_ana      = new PolAnalysis(detName, nScint);
    gInsp_ana_inel = new PolAnalysis(detName, nScint);
    configureCuts(gInsp_ana);
    configureCuts(gInsp_ana_inel);

    if (spin >= 0) {
        cout << "  Filling SpinUp elastic (full statistics)..." << endl;
        InspFillTree(tEU, gInsp_hist,      gInsp_ana, nScint, true);
        InspFillTree(tEU, gInsp_hist_elas, gInsp_ana, nScint, true);
        cout << "  Filling SpinUp inelastic (full statistics)..." << endl;
        InspFillTree(tIU, gInsp_hist,      gInsp_ana_inel, nScint, true);
        InspFillTree(tIU, gInsp_hist_inel, gInsp_ana_inel, nScint, true);
    }
    if (spin <= 0) {
        cout << "  Filling SpinDown elastic (full statistics)..." << endl;
        InspFillTree(tED, gInsp_hist,      gInsp_ana, nScint, false);
        InspFillTree(tED, gInsp_hist_elas, gInsp_ana, nScint, false);
        cout << "  Filling SpinDown inelastic (full statistics)..." << endl;
        InspFillTree(tID, gInsp_hist,      gInsp_ana_inel, nScint, false);
        InspFillTree(tID, gInsp_hist_inel, gInsp_ana_inel, nScint, false);
    }

    // Recount (multithreaded via RDataFrame) for the printed summary
    gInsp_ana->ResetCounts();
    gInsp_ana_inel->ResetCounts();
    if (spin >= 0) { gInsp_ana->CountEvents(tEU, true);  gInsp_ana_inel->CountEvents(tIU, true);  }
    if (spin <= 0) { gInsp_ana->CountEvents(tED, false); gInsp_ana_inel->CountEvents(tID, false); }

    cout << "\n########## ELASTIC ##########" << endl;
    gInsp_ana->PrintResults();
    cout << "########## INELASTIC (full statistics) ##########" << endl;
    gInsp_ana_inel->PrintResults();

    // Cross-section-weighted (pooled-count) combination — same helper as analyze_blpmc_scan()
    {
        WeightedAyResult wr = ComputeWeightedAy(gInsp_ana, gInsp_ana_inel, inelFac, beamPol);
        Double_t theoAy = isP16 ? kP16_TheoAy : kP12_TheoAy;
        cout << "\n########## WEIGHTED (pooled elastic + R0·inelastic counts) ##########" << endl;
        cout << Form("  R0 = %.6f", wr.R0) << endl;
        cout << Form("  Pooled counts:  L↑=%.2f  R↑=%.2f  L↓=%.2f  R↓=%.2f",
                     wr.Lu, wr.Ru, wr.Ld, wr.Rd) << endl;
        cout << Form("  A_y = %.6f ± %.6f", wr.Ay, wr.dAy) << endl;
        if (theoAy != 0 && wr.dAy > 0)
            cout << Form("  Deviation from theory (%.3f): %.2f sigma",
                         theoAy, TMath::Abs(wr.Ay - theoAy) / wr.dAy) << endl;
    }

    SafeClose(hEU); SafeClose(hED);
    SafeClose(hIU); SafeClose(hID);

    cout << "\n=== Inspection ready! ===" << endl;
    cout << "  insp_draw_1d_raw()        - Raw 1D spectra" << endl;
    cout << "  insp_draw_1d_comparison() - Raw vs cut overlay with cut lines" << endl;
    cout << "  insp_draw_1d_spin()       - Spin-up vs spin-down after cuts" << endl;
    cout << "  insp_draw_2d()            - 2D correlation plots" << endl;
    cout << "  insp_draw_elas_inel()     - Elastic vs inelastic (scaled to R) RAW spectra with cut lines" << endl;
    cout << "  insp_draw_elas_inel_cuts()- Elastic vs inelastic (scaled to R) POST-CUT spectra" << endl;
    cout << "  insp_draw_all()           - All of the above" << endl;
}

// =========================================================
// ─── SECTION 10: INSPECTION DRAW WRAPPERS ────────────────
// =========================================================
void insp_draw_1d_raw() {
    if (!gInsp_hist) { cout << "Run inspect() first!" << endl; return; }
    TCanvas *c = new TCanvas(Form("c_raw_%s",  gInsp_label.Data()),
                             Form("Raw Spectra — %s", gInsp_label.Data()),
                             1400, 700);
    gInsp_hist->DrawRaw(c);
}

void insp_draw_1d_comparison() {
    if (!gInsp_hist || !gInsp_ana) { cout << "Run inspect() first!" << endl; return; }
    TCanvas *c = new TCanvas(Form("c_comp_%s", gInsp_label.Data()),
                             Form("Raw vs Cuts — %s", gInsp_label.Data()),
                             1400, 700);
    gInsp_hist->DrawComparison(c, gInsp_ana);
}

void insp_draw_1d_spin() {
    if (!gInsp_hist || !gInsp_ana) { cout << "Run inspect() first!" << endl; return; }
    TCanvas *c = new TCanvas(Form("c_spin_%s", gInsp_label.Data()),
                             Form("Spin#uparrow vs Spin#downarrow — %s", gInsp_label.Data()),
                             1400, 700);
    gInsp_hist->DrawSpinComparison(c, gInsp_ana);
}

void insp_draw_2d() {
    if (!gInsp_hist) { cout << "Run inspect() first!" << endl; return; }
    TCanvas *c = new TCanvas(Form("c_2d_%s",   gInsp_label.Data()),
                             Form("2D Correlations — %s", gInsp_label.Data()),
                             1400, 700);
    gInsp_hist->Draw2D(c);
}

void insp_draw_elas_inel() {
    if (!gInsp_hist_elas || !gInsp_hist_inel || !gInsp_ana) {
        cout << "Run inspect() first!" << endl; return;
    }
    TCanvas *c = new TCanvas(Form("c_elasinel_%s", gInsp_label.Data()),
                             Form("Elastic vs Inelastic (scaled to R) RAW — %s", gInsp_label.Data()),
                             1400, 700);
    gInsp_hist_elas->DrawElasInelComparison(c, gInsp_hist_inel, gInsp_ana, gInsp_R);
}

void insp_draw_elas_inel_cuts() {
    if (!gInsp_hist_elas || !gInsp_hist_inel || !gInsp_ana) {
        cout << "Run inspect() first!" << endl; return;
    }
    TCanvas *c = new TCanvas(Form("c_elasinelcut_%s", gInsp_label.Data()),
                             Form("Elastic vs Inelastic (scaled to R) AFTER CUTS — %s", gInsp_label.Data()),
                             1400, 700);
    gInsp_hist_elas->DrawElasInelCutComparison(c, gInsp_hist_inel, gInsp_ana, gInsp_R);
}

void insp_draw_all() {
    insp_draw_1d_raw();
    insp_draw_1d_comparison();
    insp_draw_1d_spin();
    insp_draw_2d();
    insp_draw_elas_inel();
    insp_draw_elas_inel_cuts();
}



// =========================================================
// ─── SECTION 11: COUNT-VS-PARAM PLOT (cache reader) ──────
//
// Plots scintillator counts vs scan parameter using the cache
// filled by analyze_blpmc_scan() — NO FILES ARE OPENED HERE.
// You must run analyze_blpmc_scan() at least once before calling
// this function; if the cache is empty it will tell you to do so.
//
// Two families of curves per scintillator pad, each independently
// selectable by reaction via `mode` ("elas", "inel", or "both"):
//
//  RAW counts   (branch > 0, single scintillator)
//    — per-PMT occupancy, independent of the analysis cuts.
//
//  FULL-CUT counts (showCuts=true, all scintillators in
//  coincidence within their energy windows — i.e. exactly the
//  same cut PolAnalysis::CountEvents applies during the scan)
//    — this is the number that actually feeds nL_up/nR_up/etc.
//      Compare raw vs full-cut for the SAME reaction to see how
//      much of that reaction's occupancy survives the coincidence
//      cut; compare elas vs inel (mode="both") to check whether
//      inelastic events leak through the elastic selection cuts.
//    — Inelastic counts (raw and full-cut) reflect the FULL
//      inelastic sample (no subsampling) — the same full-statistics
//      counts used to compute ε_inel before the cross-section-ratio
//      weighting is applied.
//
// Usage (after running analyze_blpmc_scan() at least once):
//   draw_count_vs_param()                  // elastic raw only (default)
//   draw_count_vs_param("elas", true)      // elastic raw + elastic full-cut
//   draw_count_vs_param("inel")            // inelastic raw only (full statistics)
//   draw_count_vs_param("inel", true)      // inelastic raw + inelastic full-cut
//   draw_count_vs_param("both", true)      // elastic + inelastic, raw + full-cut (16 curves/pad)
// =========================================================
void draw_count_vs_param(TString mode = "elas", Bool_t showCuts = false) {

    if (!gCountsValid_P16 && !gCountsValid_P12) {
        cout << "[Error] No cached count data. Run analyze_blpmc_scan() first!" << endl;
        return;
    }

    Bool_t wantElas     = (mode == "elas" || mode == "both");
    Bool_t wantInel     = (mode == "inel" || mode == "both");
    Bool_t wantElasCuts = showCuts && wantElas;
    Bool_t wantInelCuts = showCuts && wantInel;

    TString kScanParam, kScanUnit;
    ParseScanDir(gScanDir, kScanParam, kScanUnit);

    struct DetCfg {
        TString name, canvName;
        Int_t   nScint;
        Bool_t  valid;
        std::vector<TGraphErrors*> *cache;
    };

    std::vector<DetCfg> dets = {
        { "P12", "cCounts_P12", 3, gCountsValid_P12, &gCountGraphs_P12 },
        { "P16", "cCounts_P16", 4, gCountsValid_P16, &gCountGraphs_P16 }
    };

    for (auto &det : dets) {
        if (!det.valid) {
            cout << "[Warning] No cached counts for " << det.name
                 << " — was it included in the last analyze_blpmc_scan() run?" << endl;
            continue;
        }

        std::vector<TGraphErrors*> &cg = *det.cache;
        Int_t ns = det.nScint;

        // Style the cached graphs (styling is cheap; do it fresh each draw
        // call so colour/marker choices stay consistent regardless of how
        // the cache was populated). Layout: 0-3 elas raw, 4-7 elas cut,
        // 8-11 inel raw, 12-15 inel cut.
        for (Int_t is = 0; is < ns; is++) {
            // Elastic raw — blue (L) / red (R)
            cg[is*16+0]->SetMarkerStyle(20); cg[is*16+0]->SetMarkerColor(kBlue+1);
            cg[is*16+0]->SetLineColor(kBlue+1); cg[is*16+0]->SetLineStyle(1);
            cg[is*16+1]->SetMarkerStyle(24); cg[is*16+1]->SetMarkerColor(kBlue+1);
            cg[is*16+1]->SetLineColor(kBlue+1); cg[is*16+1]->SetLineStyle(2);
            cg[is*16+2]->SetMarkerStyle(21); cg[is*16+2]->SetMarkerColor(kRed+1);
            cg[is*16+2]->SetLineColor(kRed+1); cg[is*16+2]->SetLineStyle(1);
            cg[is*16+3]->SetMarkerStyle(25); cg[is*16+3]->SetMarkerColor(kRed+1);
            cg[is*16+3]->SetLineColor(kRed+1); cg[is*16+3]->SetLineStyle(2);

            // Elastic full-cut — green (L) / spring-green (R)
            cg[is*16+4]->SetMarkerStyle(22); cg[is*16+4]->SetMarkerColor(kGreen+2);
            cg[is*16+4]->SetLineColor(kGreen+2); cg[is*16+4]->SetLineStyle(1);
            cg[is*16+5]->SetMarkerStyle(26); cg[is*16+5]->SetMarkerColor(kGreen+2);
            cg[is*16+5]->SetLineColor(kGreen+2); cg[is*16+5]->SetLineStyle(2);
            cg[is*16+6]->SetMarkerStyle(23); cg[is*16+6]->SetMarkerColor(kSpring+10);
            cg[is*16+6]->SetLineColor(kSpring+10); cg[is*16+6]->SetLineStyle(1);
            cg[is*16+7]->SetMarkerStyle(32); cg[is*16+7]->SetMarkerColor(kSpring+10);
            cg[is*16+7]->SetLineColor(kSpring+10); cg[is*16+7]->SetLineStyle(2);

            // Inelastic raw (full statistics) — cyan (L) / azure (R)
            cg[is*16+8]->SetMarkerStyle(29);  cg[is*16+8]->SetMarkerColor(kCyan+2);
            cg[is*16+8]->SetLineColor(kCyan+2); cg[is*16+8]->SetLineStyle(1);
            cg[is*16+9]->SetMarkerStyle(30);  cg[is*16+9]->SetMarkerColor(kCyan+2);
            cg[is*16+9]->SetLineColor(kCyan+2); cg[is*16+9]->SetLineStyle(2);
            cg[is*16+10]->SetMarkerStyle(29); cg[is*16+10]->SetMarkerColor(kAzure+2);
            cg[is*16+10]->SetLineColor(kAzure+2); cg[is*16+10]->SetLineStyle(1);
            cg[is*16+11]->SetMarkerStyle(30); cg[is*16+11]->SetMarkerColor(kAzure+2);
            cg[is*16+11]->SetLineColor(kAzure+2); cg[is*16+11]->SetLineStyle(2);

            // Inelastic full-cut (full statistics) — orange (L) / magenta (R)
            cg[is*16+12]->SetMarkerStyle(33); cg[is*16+12]->SetMarkerColor(kOrange+7);
            cg[is*16+12]->SetLineColor(kOrange+7); cg[is*16+12]->SetLineStyle(1);
            cg[is*16+13]->SetMarkerStyle(27); cg[is*16+13]->SetMarkerColor(kOrange+7);
            cg[is*16+13]->SetLineColor(kOrange+7); cg[is*16+13]->SetLineStyle(2);
            cg[is*16+14]->SetMarkerStyle(34); cg[is*16+14]->SetMarkerColor(kMagenta+1);
            cg[is*16+14]->SetLineColor(kMagenta+1); cg[is*16+14]->SetLineStyle(1);
            cg[is*16+15]->SetMarkerStyle(28); cg[is*16+15]->SetMarkerColor(kMagenta+1);
            cg[is*16+15]->SetLineColor(kMagenta+1); cg[is*16+15]->SetLineStyle(2);
        }

        TString titleSuffix = Form(" (%s%s%s)",
                                   wantElas ? "elas" : "",
                                   (wantElas && wantInel) ? "+" : "",
                                   wantInel ? "inel" : "");
        if (showCuts) titleSuffix.Insert(titleSuffix.Length()-1, " raw+cut");
        else          titleSuffix.Insert(titleSuffix.Length()-1, " raw");

        TCanvas *c = new TCanvas(det.canvName.Data(),
                                 Form("Scintillator counts vs %s — %s%s",
                                      kScanParam.Data(), det.name.Data(), titleSuffix.Data()),
                                 400 * ns, 500);
        c->Divide(ns, 1);

        for (Int_t is = 0; is < ns; is++) {
            c->cd(is + 1);
            gPad->SetGrid();
            gPad->SetLeftMargin(0.14);
            gPad->SetBottomMargin(0.14);
            gPad->SetLogy();

            TMultiGraph *mg = new TMultiGraph();
            mg->SetTitle(Form("%s S%d;%s [%s];Counts",
                              det.name.Data(), is, kScanParam.Data(), kScanUnit.Data()));

            if (wantElas)     for (Int_t ig = 0;  ig < 4;  ig++) mg->Add(cg[is*16+ig], "PL");
            if (wantElasCuts) for (Int_t ig = 4;  ig < 8;  ig++) mg->Add(cg[is*16+ig], "PL");
            if (wantInel)     for (Int_t ig = 8;  ig < 12; ig++) mg->Add(cg[is*16+ig], "PL");
            if (wantInelCuts) for (Int_t ig = 12; ig < 16; ig++) mg->Add(cg[is*16+ig], "PL");

            mg->Draw("A");
            mg->GetXaxis()->SetTitle(Form("%s [%s]", kScanParam.Data(), kScanUnit.Data()));
            mg->GetYaxis()->SetTitle("Counts");
            mg->GetXaxis()->SetTitleSize(0.05);
            mg->GetYaxis()->SetTitleSize(0.05);
            mg->GetYaxis()->SetTitleOffset(1.2);

            if (is == 0) {
                Int_t nRows = (wantElas?1:0) + (wantElasCuts?1:0) + (wantInel?1:0) + (wantInelCuts?1:0);
                Double_t yTop = 0.88;
                Double_t yBot = yTop - 0.14 * nRows - 0.03;
                TLegend *leg = new TLegend(0.13, TMath::Max(yBot, 0.30), 0.62, yTop);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.032);
                if (wantElas) {
                    leg->AddEntry(cg[is*16+0], "Raw Elas L  #uparrow",   "lep");
                    leg->AddEntry(cg[is*16+1], "Raw Elas L  #downarrow", "lep");
                    leg->AddEntry(cg[is*16+2], "Raw Elas R  #uparrow",   "lep");
                    leg->AddEntry(cg[is*16+3], "Raw Elas R  #downarrow", "lep");
                }
                if (wantElasCuts) {
                    leg->AddEntry(cg[is*16+4], "Cut Elas L #uparrow",   "lep");
                    leg->AddEntry(cg[is*16+5], "Cut Elas L #downarrow", "lep");
                    leg->AddEntry(cg[is*16+6], "Cut Elas R #uparrow",   "lep");
                    leg->AddEntry(cg[is*16+7], "Cut Elas R #downarrow", "lep");
                }
                if (wantInel) {
                    leg->AddEntry(cg[is*16+8],  "Raw Inel(full) L  #uparrow",   "lep");
                    leg->AddEntry(cg[is*16+9],  "Raw Inel(full) L  #downarrow", "lep");
                    leg->AddEntry(cg[is*16+10], "Raw Inel(full) R  #uparrow",   "lep");
                    leg->AddEntry(cg[is*16+11], "Raw Inel(full) R  #downarrow", "lep");
                }
                if (wantInelCuts) {
                    leg->AddEntry(cg[is*16+12], "Cut Inel(full) L #uparrow",   "lep");
                    leg->AddEntry(cg[is*16+13], "Cut Inel(full) L #downarrow", "lep");
                    leg->AddEntry(cg[is*16+14], "Cut Inel(full) R #uparrow",   "lep");
                    leg->AddEntry(cg[is*16+15], "Cut Inel(full) R #downarrow", "lep");
                }
                leg->Draw();
            }
        }
        c->Update();
        cout << "Canvas " << det.canvName << " created (from cache, no files reopened)." << endl;
    }
}

// =========================================================
// ─── SECTION 12: SCAN CACHE PERSISTENCE (SAVE/LOAD) ──────
//
// Everything analyze_blpmc_scan() computes lives only in the runtime
// globals above (gGraph_P16/P12_*, gCountGraphs_P16/P12, gInelFactor_*,
// gBeamPol, gScanDir) — normally lost when the ROOT session ends, forcing
// a full re-scan (many files, many minutes) just to look at a plot again.
//
// SaveScanCache() / LoadScanCache() write/read exactly those globals to
// a plain TFile, so a later session can call LoadScanCache() and go
// straight to draw_scan_p16()/draw_scan_p12()/draw_scan_both()/
// draw_count_vs_param() — with ZERO Geant4 files reopened.
//
// This is saved ONLY on explicit request — analyze_blpmc_scan() does
// NOT auto-save. Call SaveScanCache() yourself when you're happy with
// a scan's results.
//
// inspect()'s per-point histograms (gInsp_hist, gInsp_ana, ...) are
// NOT included — a single inspect() call only opens 4 files and is
// already fast, so there's no re-scan cost to avoid there.
//
// File layout (plain TFile, inspectable with a TBrowser):
//   Top-level metadata (TNamed / TParameter<T>):
//     gScanDir, kScanParam, kScanUnit, gInelFactor_P16, gInelFactor_P12,
//     gBeamPol, gCountsValid_P16, gCountsValid_P12
//   A_y graphs, already named as they exist in memory:
//     gAy_P16_elas, gAy_P16_inel, gAy_P16_weighted,
//     gAy_P12_elas, gAy_P12_inel, gAy_P12_weighted
//   Count-cache graphs, named on the way out (slot layout matches the
//   comment above the gCountGraphs_P16/P12 declarations):
//     cnt_P16_S<is>_slot<00-15>, cnt_P12_S<is>_slot<00-15>
//
// Usage:
//   SaveScanCache()                 // auto: Cache/<gScanDir>.root
//   SaveScanCache("myscan.root")    // explicit filename
//   LoadScanCache()                 // needs gScanDir set; auto path as above
//   LoadScanCache("myscan.root")    // explicit filename
// =========================================================

void SaveScanCache(TString filename = "") {
    if (gScanDir == "") {
        cout << "[Error] gScanDir is not set — nothing to identify this cache by." << endl;
        return;
    }
    Bool_t haveAnyAy = gGraph_P16_weighted || gGraph_P12_weighted;
    if (!haveAnyAy && !gCountsValid_P16 && !gCountsValid_P12) {
        cout << "[Error] No scan data in memory. Run analyze_blpmc_scan() first!" << endl;
        return;
    }

    if (filename == "") {
        gSystem->mkdir("Cache", kTRUE);   // kTRUE = create parent dirs too
        filename = "Cache/" + gScanDir + ".root";
    }

    TFile *f = TFile::Open(filename.Data(), "RECREATE");
    if (!f || f->IsZombie()) {
        cout << "[Error] Could not open \"" << filename << "\" for writing." << endl;
        return;
    }

    // ── Metadata ──
    TString kScanParam, kScanUnit;
    ParseScanDir(gScanDir, kScanParam, kScanUnit);

    TNamed(  "gScanDir",         gScanDir.Data()  ).Write();
    TNamed(  "kScanParam",       kScanParam.Data()).Write();
    TNamed(  "kScanUnit",        kScanUnit.Data() ).Write();
    TParameter<Double_t>("gInelFactor_P16",  gInelFactor_P16      ).Write();
    TParameter<Double_t>("gInelFactor_P12",  gInelFactor_P12      ).Write();
    TParameter<Double_t>("gBeamPol",         gBeamPol             ).Write();
    TParameter<Int_t>(   "gCountsValid_P16", (Int_t)gCountsValid_P16).Write();
    TParameter<Int_t>(   "gCountsValid_P12", (Int_t)gCountsValid_P12).Write();

    // ── A_y graphs — already named, just write them ──
    auto writeGraph = [](TGraphErrors *g) { if (g) g->Write(); };
    writeGraph(gGraph_P16_elas); writeGraph(gGraph_P16_inel); writeGraph(gGraph_P16_weighted);
    writeGraph(gGraph_P12_elas); writeGraph(gGraph_P12_inel); writeGraph(gGraph_P12_weighted);

    // ── Count-cache graphs — name them on the way out ──
    auto writeCountGraphs = [](std::vector<TGraphErrors*> &cg, const TString &det,
                               Int_t nScint, Bool_t valid) {
        if (!valid) return;
        for (Int_t is = 0; is < nScint; is++) {
            for (Int_t slot = 0; slot < 16; slot++) {
                TGraphErrors *g = cg[is*16+slot];
                if (!g) continue;
                g->Write(Form("cnt_%s_S%d_slot%02d", det.Data(), is, slot));
            }
        }
    };
    writeCountGraphs(gCountGraphs_P16, "P16", 4, gCountsValid_P16);
    writeCountGraphs(gCountGraphs_P12, "P12", 3, gCountsValid_P12);

    f->Close();
    delete f;

    cout << "\n=== Scan cache saved: " << filename << " ===" << endl;
    cout << "  gScanDir       = " << gScanDir << endl;
    cout << "  A_y graphs     = " << (haveAnyAy ? "yes" : "no") << endl;
    cout << "  P16 count-cache = " << (gCountsValid_P16 ? "yes" : "no (not in this scan)") << endl;
    cout << "  P12 count-cache = " << (gCountsValid_P12 ? "yes" : "no (not in this scan)") << endl;
}

void LoadScanCache(TString filename = "") {
    if (filename == "") {
        if (gScanDir == "") {
            cout << "[Error] No filename given and gScanDir is not set — cannot guess the cache path." << endl;
            cout << "  Either set gScanDir = \"...\" first, or call LoadScanCache(\"path/to/file.root\")." << endl;
            return;
        }
        filename = "Cache/" + gScanDir + ".root";
    }

    TFile *f = TFile::Open(filename.Data(), "READ");
    if (!f || f->IsZombie()) {
        cout << "[Error] Could not open cache file: " << filename << endl;
        return;
    }

    // ── Metadata ──
    auto getNamed = [&](const char *key) -> TString {
        TNamed *n = (TNamed*)f->Get(key);
        return n ? TString(n->GetTitle()) : TString("");
    };
    auto getParamD = [&](const char *key, Double_t fallback) -> Double_t {
        TParameter<Double_t> *p = (TParameter<Double_t>*)f->Get(key);
        return p ? p->GetVal() : fallback;
    };
    auto getParamI = [&](const char *key, Int_t fallback) -> Int_t {
        TParameter<Int_t> *p = (TParameter<Int_t>*)f->Get(key);
        return p ? p->GetVal() : fallback;
    };

    TString savedScanDir = getNamed("gScanDir");
    if (savedScanDir == "") {
        cout << "[Warning] Cache file has no gScanDir metadata — is this really a scan cache file?" << endl;
    } else if (gScanDir != "" && savedScanDir != gScanDir) {
        cout << "[Warning] Cache's gScanDir (\"" << savedScanDir << "\") differs from the current"
             << " gScanDir (\"" << gScanDir << "\") — loading anyway, but check you meant to." << endl;
    }
    if (savedScanDir != "") gScanDir = savedScanDir;

    gInelFactor_P16  = getParamD("gInelFactor_P16", -1.0);
    gInelFactor_P12  = getParamD("gInelFactor_P12", -1.0);
    gBeamPol         = getParamD("gBeamPol",        -1.0);
    gCountsValid_P16 = (Bool_t)getParamI("gCountsValid_P16", 0);
    gCountsValid_P12 = (Bool_t)getParamI("gCountsValid_P12", 0);

    // ── Clear whatever is currently in memory before loading ──
    auto safeDelG = [](TGraphErrors* &g) { if (g) { delete g; g = nullptr; } };
    safeDelG(gGraph_P16_elas);     safeDelG(gGraph_P16_inel);     safeDelG(gGraph_P16_weighted);
    safeDelG(gGraph_P12_elas);     safeDelG(gGraph_P12_inel);     safeDelG(gGraph_P12_weighted);

    auto freeCountGraphs = [](std::vector<TGraphErrors*> &v) {
        for (auto &g : v) if (g) { delete g; g = nullptr; }
        v.clear();
    };
    freeCountGraphs(gCountGraphs_P16);
    freeCountGraphs(gCountGraphs_P12);

    // ── A_y graphs — Clone() so they survive f->Close() below ──
    auto getGraph = [&](const char *key) -> TGraphErrors* {
        TGraphErrors *g = (TGraphErrors*)f->Get(key);
        return g ? (TGraphErrors*)g->Clone() : nullptr;
    };
    gGraph_P16_elas     = getGraph("gAy_P16_elas");
    gGraph_P16_inel     = getGraph("gAy_P16_inel");
    gGraph_P16_weighted = getGraph("gAy_P16_weighted");
    gGraph_P12_elas     = getGraph("gAy_P12_elas");
    gGraph_P12_inel     = getGraph("gAy_P12_inel");
    gGraph_P12_weighted = getGraph("gAy_P12_weighted");

    if (!gGraph_P16_weighted && !gGraph_P12_weighted)
        cout << "[Warning] Neither P16 nor P12 weighted A_y graph found in cache — file may be incomplete." << endl;

    // ── Count-cache graphs ──
    auto loadCountGraphs = [&](std::vector<TGraphErrors*> &cg, const TString &det,
                               Int_t nScint, Bool_t valid) {
        if (!valid) return;
        cg.assign(nScint * 16, nullptr);
        Int_t nMissing = 0;
        for (Int_t is = 0; is < nScint; is++) {
            for (Int_t slot = 0; slot < 16; slot++) {
                TString key = Form("cnt_%s_S%d_slot%02d", det.Data(), is, slot);
                TGraphErrors *g = (TGraphErrors*)f->Get(key.Data());
                if (!g) { nMissing++; continue; }
                cg[is*16+slot] = (TGraphErrors*)g->Clone();
            }
        }
        if (nMissing > 0)
            cout << "[Warning] " << det << " count-cache: " << nMissing << " of " << (nScint*16)
                 << " expected graphs missing from file!" << endl;
    };
    loadCountGraphs(gCountGraphs_P16, "P16", 4, gCountsValid_P16);
    loadCountGraphs(gCountGraphs_P12, "P12", 3, gCountsValid_P12);

    f->Close();
    delete f;

    cout << "\n=== Scan cache loaded: " << filename << " ===" << endl;
    cout << "  gScanDir = " << gScanDir << endl;
    cout << Form("  R_P16 = %.6f   R_P12 = %.6f   P_beam = %.4f",
                 gInelFactor_P16, gInelFactor_P12, gBeamPol) << endl;
    cout << "  P16 count-cache: " << (gCountsValid_P16 ? "loaded" : "not present in cache") << endl;
    cout << "  P12 count-cache: " << (gCountsValid_P12 ? "loaded" : "not present in cache") << endl;
    cout << "No files reopened — call draw_scan_p16()/p12()/both() or draw_count_vs_param() now." << endl;
}