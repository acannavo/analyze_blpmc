#ifndef POLANALYSIS_H
#define POLANALYSIS_H

#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include <iostream>
#include <vector>

class PolAnalysis {
private:
    // Event counts
    Long64_t nL_up, nR_up;       // Spin-up: left and right arm counts
    Long64_t nL_down, nR_down;   // Spin-down: left and right arm counts

    // Configuration
    Double_t beamPolarization;
    Double_t beamPolarizationError;
    Double_t theoreticalAy;

    // Selection cuts (4 scintillators for P16, 3 for P12)
    std::vector<Double_t> cutMin;
    std::vector<Double_t> cutMax;
    Int_t nScintillators;
    TString detectorName;

public:
    PolAnalysis(TString name, Int_t nScint) {
        detectorName   = name;
        nScintillators = nScint;
        cutMin.resize(nScintillators);
        cutMax.resize(nScintillators);

        nL_up = 0; nR_up = 0;
        nL_down = 0; nR_down = 0;

        beamPolarization      = 1.0;
        beamPolarizationError = 0.0;
        theoreticalAy         = 0.0;
    }

    // -------------------------------------------------------
    // Configuration
    // -------------------------------------------------------
    void SetCuts(Int_t scint, Double_t min, Double_t max) {
        if (scint >= 0 && scint < nScintillators) {
            cutMin[scint] = min;
            cutMax[scint] = max;
        }
    }

    void SetBeamPolarization(Double_t pol)      { beamPolarization      = pol; }
    void SetBeamPolarizationError(Double_t err) { beamPolarizationError = err; }
    void SetTheoreticalAy(Double_t ay)          { theoreticalAy         = ay;  }

    // -------------------------------------------------------
    // Accessors (needed by PolHistograms)
    // -------------------------------------------------------
    Double_t GetCutMin(Int_t scint) const {
        return (scint >= 0 && scint < nScintillators) ? cutMin[scint] : 0.0;
    }
    Double_t GetCutMax(Int_t scint) const {
        return (scint >= 0 && scint < nScintillators) ? cutMax[scint] : 0.0;
    }
    Int_t GetNScintillators() const { return nScintillators; }

    // -------------------------------------------------------
    // Reset counters
    // -------------------------------------------------------
    void ResetCounts() {
        nL_up = 0; nR_up = 0;
        nL_down = 0; nR_down = 0;
    }

    // -------------------------------------------------------
    // CountEvents
    //
    // Uses TTree::Draw with cut strings for speed, but avoids
    // the ">>+" histogram accumulation syntax which triggers a
    // Cling JIT symbol lookup that crashes on ROOT 6.32.
    // Instead, a unique histogram name is generated each call
    // so ROOT always creates a fresh object with no lookup.
    // The temporary histograms are deleted immediately after
    // use to prevent memory buildup during long scans.
    //
    // maxEntries = -1 uses all entries;
    // maxEntries > 0  truncates (for inelastic normalization)
    // -------------------------------------------------------
    void CountEvents(TTree* tree, Bool_t isSpinUp, Long64_t maxEntries = -1) {
        if (!tree) {
            std::cout << "Error: Null tree pointer!" << std::endl;
            return;
        }

        // Build cut strings
        TString cutL = "", cutR = "";

        if (detectorName == "P16") {
            cutL = Form("eDepP16LS0>=%g && eDepP16LS0<=%g"
                        " && eDepP16LS1>=%g && eDepP16LS1<=%g"
                        " && eDepP16LS2>=%g && eDepP16LS2<=%g"
                        " && eDepP16LS3>=%g && eDepP16LS3<=%g",
                        cutMin[0], cutMax[0], cutMin[1], cutMax[1],
                        cutMin[2], cutMax[2], cutMin[3], cutMax[3]);
            cutR = Form("eDepP16RS0>=%g && eDepP16RS0<=%g"
                        " && eDepP16RS1>=%g && eDepP16RS1<=%g"
                        " && eDepP16RS2>=%g && eDepP16RS2<=%g"
                        " && eDepP16RS3>=%g && eDepP16RS3<=%g",
                        cutMin[0], cutMax[0], cutMin[1], cutMax[1],
                        cutMin[2], cutMax[2], cutMin[3], cutMax[3]);
        } else {
            cutL = Form("eDepP12LS0>=%g && eDepP12LS0<=%g"
                        " && eDepP12LS1>=%g && eDepP12LS1<=%g"
                        " && eDepP12LS2>=%g && eDepP12LS2<=%g",
                        cutMin[0], cutMax[0], cutMin[1], cutMax[1],
                        cutMin[2], cutMax[2]);
            cutR = Form("eDepP12RS0>=%g && eDepP12RS0<=%g"
                        " && eDepP12RS1>=%g && eDepP12RS1<=%g"
                        " && eDepP12RS2>=%g && eDepP12RS2<=%g",
                        cutMin[0], cutMax[0], cutMin[1], cutMax[1],
                        cutMin[2], cutMax[2]);
        }

        Long64_t nTotal = tree->GetEntries();
        Long64_t nUse   = (maxEntries > 0 && maxEntries < nTotal) ? maxEntries : nTotal;

        // Unique histogram names per call — avoids the ">>" accumulation
        // operator which requires Cling to look up an existing symbol and
        // crashes on ROOT 6.32 under heavy interpreter load.
        static Long64_t sDrawCount = 0;
        TString hnameL = Form("_hCntL_%lld", sDrawCount);
        TString hnameR = Form("_hCntR_%lld", sDrawCount);
        sDrawCount++;

        tree->Draw(Form("1>>%s(1,0,2)", hnameL.Data()), cutL, "goff", nUse, 0);
        Long64_t countL = (Long64_t)tree->GetSelectedRows();

        tree->Draw(Form("1>>%s(1,0,2)", hnameR.Data()), cutR, "goff", nUse, 0);
        Long64_t countR = (Long64_t)tree->GetSelectedRows();

        // Delete temporary histograms immediately to avoid memory buildup
        // across many CountEvents calls in a full parameter scan
        if (gDirectory) {
            TH1 *hL = (TH1*)gDirectory->FindObject(hnameL.Data());
            TH1 *hR = (TH1*)gDirectory->FindObject(hnameR.Data());
            if (hL) delete hL;
            if (hR) delete hR;
        }

        if (isSpinUp) { nL_up += countL; nR_up += countR; }
        else          { nL_down += countL; nR_down += countR; }
    }

    // -------------------------------------------------------
    // Asymmetry calculations
    // -------------------------------------------------------
    Double_t GetSimpleAsymmetry(Bool_t isSpinUp) {
        Long64_t nL = isSpinUp ? nL_up : nL_down;
        Long64_t nR = isSpinUp ? nR_up : nR_down;
        if (nL + nR == 0) return 0.0;
        return (1.0 * nL - 1.0 * nR) / (1.0 * nL + 1.0 * nR);
    }

    Double_t GetCrossRatioAsymmetry() {
        if (nL_up == 0 || nR_up == 0 || nL_down == 0 || nR_down == 0) {
            std::cout << "Warning: Zero counts detected, cannot compute cross-ratio!" << std::endl;
            return 0.0;
        }
        Double_t sqrt_LR = TMath::Sqrt(1.0 * nL_up * nR_down);
        Double_t sqrt_RL = TMath::Sqrt(1.0 * nR_up * nL_down);
        return (sqrt_LR - sqrt_RL) / (sqrt_LR + sqrt_RL);
    }

    // -------------------------------------------------------
    // Asymmetry error calculations
    // -------------------------------------------------------
    Double_t GetSimpleAsymmetryError(Bool_t isSpinUp) {
        Long64_t nL = isSpinUp ? nL_up : nL_down;
        Long64_t nR = isSpinUp ? nR_up : nR_down;
        Long64_t nS = nL + nR;
        if (nS == 0) return 0.0;
        // δε = 2√(nL·nR) / nS^(3/2)
        return 2.0 * TMath::Sqrt(1.0 * nL * nR) / (nS * TMath::Sqrt(1.0 * nS));
    }

    Double_t GetCrossRatioError() {
        if (nL_up == 0 || nR_up == 0 || nL_down == 0 || nR_down == 0) return 0.0;
        Double_t A = TMath::Sqrt(1.0 * nL_up  * nR_down);
        Double_t B = TMath::Sqrt(1.0 * nR_up  * nL_down);
        Double_t S = A + B;
        Double_t numerator = TMath::Sqrt(
            B*B * (nL_up + nR_down) + A*A * (nR_up + nL_down)
        );
        return numerator / (S * S);
    }

    // -------------------------------------------------------
    // Analyzing power calculations
    // -------------------------------------------------------
    Double_t GetAnalyzingPower_Simple(Bool_t isSpinUp) {
        if (beamPolarization == 0.0) return 0.0;
        return GetSimpleAsymmetry(isSpinUp) / beamPolarization;
    }

    Double_t GetAnalyzingPower_CrossRatio() {
        if (beamPolarization == 0.0) return 0.0;
        return GetCrossRatioAsymmetry() / beamPolarization;
    }

    Double_t GetAnalyzingPowerError_Simple(Bool_t isSpinUp) {
        if (beamPolarization == 0.0) return 0.0;
        Double_t eps_error = GetSimpleAsymmetryError(isSpinUp);
        Double_t eps       = GetSimpleAsymmetry(isSpinUp);
        if (beamPolarizationError == 0.0) return eps_error / beamPolarization;
        Double_t Ay = eps / beamPolarization;
        if (eps == 0.0) return 0.0;
        Double_t rel_eps_error = eps_error / TMath::Abs(eps);
        Double_t rel_pol_error = beamPolarizationError / beamPolarization;
        return TMath::Abs(Ay) * TMath::Sqrt(rel_eps_error * rel_eps_error +
                                            rel_pol_error * rel_pol_error);
    }

    Double_t GetAnalyzingPowerError_CrossRatio() {
        if (beamPolarization == 0.0) return 0.0;
        Double_t eps_error = GetCrossRatioError();
        Double_t eps       = GetCrossRatioAsymmetry();
        if (beamPolarizationError == 0.0) return eps_error / beamPolarization;
        Double_t Ay = eps / beamPolarization;
        if (eps == 0.0) return 0.0;
        Double_t rel_eps_error = eps_error / eps;
        Double_t rel_pol_error = beamPolarizationError / beamPolarization;
        return Ay * TMath::Sqrt(rel_eps_error * rel_eps_error +
                                rel_pol_error * rel_pol_error);
    }

    // -------------------------------------------------------
    // Significance of deviation from theory
    // -------------------------------------------------------
    Double_t GetSignificance() {
        if (theoreticalAy == 0.0) return 0.0;
        Double_t Ay_measured = GetAnalyzingPower_CrossRatio();
        Double_t Ay_error    = GetAnalyzingPowerError_CrossRatio();
        if (Ay_error == 0.0) return 0.0;
        return TMath::Abs(Ay_measured - theoreticalAy) / Ay_error;
    }

    // -------------------------------------------------------
    // Accessors for counts
    // -------------------------------------------------------
    Long64_t GetCountLeft(Bool_t isSpinUp)  { return isSpinUp ? nL_up : nL_down; }
    Long64_t GetCountRight(Bool_t isSpinUp) { return isSpinUp ? nR_up : nR_down; }
    Long64_t GetCountTotal(Bool_t isSpinUp) {
        return (isSpinUp ? nL_up : nL_down) + (isSpinUp ? nR_up : nR_down);
    }

    // -------------------------------------------------------
    // Print results
    // -------------------------------------------------------
    void PrintResults() {
        std::cout << "\n=== " << detectorName << " Analysis Results ===" << std::endl;

        std::cout << "\n--- Event Counts ---" << std::endl;
        std::cout << "Spin-Up:   Left = " << nL_up   << ", Right = " << nR_up
                  << ", Total = " << (nL_up + nR_up)   << std::endl;
        std::cout << "Spin-Down: Left = " << nL_down << ", Right = " << nR_down
                  << ", Total = " << (nL_down + nR_down) << std::endl;

        std::cout << "\n--- Simple Asymmetries ---" << std::endl;
        std::cout << "Spin-Up:   ε = " << GetSimpleAsymmetry(true)
                  << " ± " << GetSimpleAsymmetryError(true)  << std::endl;
        std::cout << "Spin-Down: ε = " << GetSimpleAsymmetry(false)
                  << " ± " << GetSimpleAsymmetryError(false) << std::endl;

        std::cout << "\n--- Cross-Ratio Method ---" << std::endl;
        std::cout << "Cross-Ratio: ε = " << GetCrossRatioAsymmetry()
                  << " ± " << GetCrossRatioError() << std::endl;

        if (beamPolarization > 0) {
            std::cout << "\n--- Analyzing Powers (P_beam = " << beamPolarization;
            if (beamPolarizationError > 0)
                std::cout << " ± " << beamPolarizationError;
            std::cout << ") ---" << std::endl;

            std::cout << "From Spin-Up:     A_y = " << GetAnalyzingPower_Simple(true)
                      << " ± " << GetAnalyzingPowerError_Simple(true)  << std::endl;
            std::cout << "From Spin-Down:   A_y = " << GetAnalyzingPower_Simple(false)
                      << " ± " << GetAnalyzingPowerError_Simple(false) << std::endl;
            std::cout << "From Cross-Ratio: A_y = " << GetAnalyzingPower_CrossRatio()
                      << " ± " << GetAnalyzingPowerError_CrossRatio()  << std::endl;

            if (theoreticalAy != 0) {
                std::cout << "\n--- Comparison to Theory ---" << std::endl;
                std::cout << "Theoretical A_y = " << theoreticalAy << std::endl;

                Double_t Ay_cross      = GetAnalyzingPower_CrossRatio();
                Double_t deviation     = Ay_cross - theoreticalAy;
                Double_t deviation_pct = (deviation / theoreticalAy) * 100.0;
                Double_t significance  = GetSignificance();

                std::cout << "Deviation: " << deviation
                          << " (" << deviation_pct << "%)" << std::endl;
                std::cout << "Statistical significance: " << significance << " σ" << std::endl;

                if      (significance < 1.0) std::cout << "  → Excellent agreement with theory!" << std::endl;
                else if (significance < 2.0) std::cout << "  → Good agreement with theory"       << std::endl;
                else if (significance < 3.0) std::cout << "  → Moderate disagreement (~2σ)"      << std::endl;
                else                         std::cout << "  → Significant disagreement (>3σ) - check systematics!" << std::endl;
            }
        }

        std::cout << "\n--- Selection Cuts ---" << std::endl;
        for (Int_t i = 0; i < nScintillators; i++)
            std::cout << "S" << i << ": " << cutMin[i] << " - " << cutMax[i] << " MeV" << std::endl;
        std::cout << "========================\n" << std::endl;
    }

    void PrintCutInfo() {
        std::cout << "\n=== " << detectorName << " Elastic Selection Cuts ===" << std::endl;
        for (Int_t i = 0; i < nScintillators; i++)
            std::cout << "S" << i << ": " << cutMin[i] << " - " << cutMax[i] << " MeV" << std::endl;
    }
};

#endif