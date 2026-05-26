#ifndef POLANALYSIS_H
#define POLANALYSIS_H

#include "TTree.h"
#include "TMath.h"
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
        detectorName = name;
        nScintillators = nScint;
        cutMin.resize(nScintillators);
        cutMax.resize(nScintillators);
        
        // Initialize counts
        nL_up = 0;
        nR_up = 0;
        nL_down = 0;
        nR_down = 0;
        
        // Default values
        beamPolarization = 1.0;  // 100% polarized (ideal)
        beamPolarizationError = 0.0;
        theoreticalAy = 0.0;     // To be set by user
    }
    
    // Configuration methods
    void SetCuts(Int_t scint, Double_t min, Double_t max) {
        if (scint >= 0 && scint < nScintillators) {
            cutMin[scint] = min;
            cutMax[scint] = max;
        }
    }
    
    void SetBeamPolarization(Double_t pol) {
        beamPolarization = pol;
    }
    
    void SetBeamPolarizationError(Double_t err) {
        beamPolarizationError = err;
    }
    
    void SetTheoreticalAy(Double_t ay) {
        theoreticalAy = ay;
    }
    
    // Accessor methods for cuts (needed by PolHistograms)
    Double_t GetCutMin(Int_t scint) const {
        if (scint >= 0 && scint < nScintillators) {
            return cutMin[scint];
        }
        return 0.0;
    }
    
    Double_t GetCutMax(Int_t scint) const {
        if (scint >= 0 && scint < nScintillators) {
            return cutMax[scint];
        }
        return 0.0;
    }
    
    Int_t GetNScintillators() const {
        return nScintillators;
    }
    
    // Reset all counters — call before reusing the same object across datasets
    void ResetCounts() {
        nL_up = 0; nR_up = 0;
        nL_down = 0; nR_down = 0;
    }

    // Event counting using TTree::GetEntries(cut) — fast ROOT engine, no manual loop
    // maxEntries = -1 uses all entries; > 0 truncates (for inelastic normalization)
    void CountEvents(TTree* tree, Bool_t isSpinUp, Long64_t maxEntries = -1) {
        if (!tree) {
            std::cout << "Error: Null tree pointer!" << std::endl;
            return;
        }

        // Build cut strings — ROOT's TTreeFormula engine evaluates these
        // much faster than a manual entry-by-entry loop
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
        } else if (detectorName == "P12") {
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

        // TTree::Draw with entry range — works in all ROOT versions
        // Draw into a temporary 1-bin histogram, count only entries in [0, nUse)
        Long64_t nTotal = tree->GetEntries();
        Long64_t nUse   = (maxEntries > 0 && maxEntries < nTotal) ? maxEntries : nTotal;

        // Use a dummy variable "1" just to trigger the cut evaluation;
        // GetSelectedRows() returns how many entries passed the cut
        tree->Draw("1>>+_hDummyL(1,0,2)", cutL, "goff", nUse, 0);
        Long64_t countL = (Long64_t)tree->GetSelectedRows();

        tree->Draw("1>>+_hDummyR(1,0,2)", cutR, "goff", nUse, 0);
        Long64_t countR = (Long64_t)tree->GetSelectedRows();

        // += accumulates across multiple CountEvents calls
        // (e.g. elastic then inelastic on the same object)
        if (isSpinUp) {
            nL_up += countL;
            nR_up += countR;
        } else {
            nL_down += countL;
            nR_down += countR;
        }
    }
    
    // Asymmetry calculations
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
    
    // Asymmetry error calculations
    Double_t GetSimpleAsymmetryError(Bool_t isSpinUp) {
        Long64_t nL = isSpinUp ? nL_up : nL_down;
        Long64_t nR = isSpinUp ? nR_up : nR_down;
        Long64_t nS = nL + nR;
        
        if (nS == 0) return 0.0;
        
        // Correct formula: δε = 2√(nL·nR) / nS^(3/2)
        return 2.0 * TMath::Sqrt(1.0 * nL * nR) / (nS * TMath::Sqrt(1.0 * nS));
    }
	
    
    Double_t GetCrossRatioError() {
        if (nL_up == 0 || nR_up == 0 || nL_down == 0 || nR_down == 0) {
            return 0.0;
        }
		Double_t A = TMath::Sqrt(1.0 * nL_up * nR_down);
		Double_t B = TMath::Sqrt(1.0 * nR_up * nL_down);
		Double_t S = A + B;

		Double_t numerator = TMath::Sqrt(
			B*B * (nL_up + nR_down) + A*A * (nR_up + nL_down)
		);

		return numerator / (S * S);
		
    }
    
    // Analyzing power calculations
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
		Double_t eps = GetSimpleAsymmetry(isSpinUp);
		
		// If no beam polarization error, simple propagation
		if (beamPolarizationError == 0.0) {
			return eps_error / beamPolarization;  // This is always positive
		}
		
		// Full error propagation
		Double_t Ay = eps / beamPolarization;
		
		if (eps == 0.0) return 0.0;
		
		Double_t rel_eps_error = eps_error / TMath::Abs(eps);  // Use absolute value!
		Double_t rel_pol_error = beamPolarizationError / beamPolarization;
		
		return TMath::Abs(Ay) * TMath::Sqrt(rel_eps_error * rel_eps_error + 
								rel_pol_error * rel_pol_error);
	}
    
    Double_t GetAnalyzingPowerError_CrossRatio() {
        if (beamPolarization == 0.0) return 0.0;
        
        Double_t eps_error = GetCrossRatioError();
        Double_t eps = GetCrossRatioAsymmetry();
        
        // If no beam polarization error, simple propagation
        if (beamPolarizationError == 0.0) {
            return eps_error / beamPolarization;
        }
        
        // Full error propagation
        Double_t Ay = eps / beamPolarization;
        
        if (eps == 0.0) return 0.0;
        
        Double_t rel_eps_error = eps_error / eps;
        Double_t rel_pol_error = beamPolarizationError / beamPolarization;
        
        return Ay * TMath::Sqrt(rel_eps_error * rel_eps_error + 
                                rel_pol_error * rel_pol_error);
    }
    
    // Calculate significance of deviation from theory
    Double_t GetSignificance() {
        if (theoreticalAy == 0.0) return 0.0;
        
        Double_t Ay_measured = GetAnalyzingPower_CrossRatio();
        Double_t Ay_error = GetAnalyzingPowerError_CrossRatio();
        
        if (Ay_error == 0.0) return 0.0;
        
        return TMath::Abs(Ay_measured - theoreticalAy) / Ay_error;
    }
    
    // Access to counts
    Long64_t GetCountLeft(Bool_t isSpinUp) { return isSpinUp ? nL_up : nL_down; }
    Long64_t GetCountRight(Bool_t isSpinUp) { return isSpinUp ? nR_up : nR_down; }
    Long64_t GetCountTotal(Bool_t isSpinUp) { 
        return (isSpinUp ? nL_up : nL_down) + (isSpinUp ? nR_up : nR_down); 
    }
    
    // Print results
    void PrintResults() {
        std::cout << "\n=== " << detectorName << " Analysis Results ===" << std::endl;
        std::cout << "\n--- Event Counts ---" << std::endl;
        std::cout << "Spin-Up:   Left = " << nL_up << ", Right = " << nR_up 
                  << ", Total = " << (nL_up + nR_up) << std::endl;
        std::cout << "Spin-Down: Left = " << nL_down << ", Right = " << nR_down 
                  << ", Total = " << (nL_down + nR_down) << std::endl;
        
        std::cout << "\n--- Simple Asymmetries ---" << std::endl;
        Double_t asym_up = GetSimpleAsymmetry(true);
        Double_t err_up = GetSimpleAsymmetryError(true);
        std::cout << "Spin-Up:   ε = " << asym_up << " ± " << err_up << std::endl;
        
        Double_t asym_down = GetSimpleAsymmetry(false);
        Double_t err_down = GetSimpleAsymmetryError(false);
        std::cout << "Spin-Down: ε = " << asym_down << " ± " << err_down << std::endl;
        
        std::cout << "\n--- Cross-Ratio Method ---" << std::endl;
        Double_t asym_cross = GetCrossRatioAsymmetry();
        Double_t err_cross = GetCrossRatioError();
        std::cout << "Cross-Ratio: ε = " << asym_cross << " ± " << err_cross << std::endl;
        
        if (beamPolarization > 0) {
            std::cout << "\n--- Analyzing Powers (P_beam = " << beamPolarization;
            if (beamPolarizationError > 0) {
                std::cout << " ± " << beamPolarizationError;
            }
            std::cout << ") ---" << std::endl;
            
            Double_t Ay_up = GetAnalyzingPower_Simple(true);
            Double_t dAy_up = GetAnalyzingPowerError_Simple(true);
            std::cout << "From Spin-Up:     A_y = " << Ay_up << " ± " << dAy_up << std::endl;
            
            Double_t Ay_down = GetAnalyzingPower_Simple(false);
            Double_t dAy_down = GetAnalyzingPowerError_Simple(false);
            std::cout << "From Spin-Down:   A_y = " << Ay_down << " ± " << dAy_down << std::endl;
            
            Double_t Ay_cross = GetAnalyzingPower_CrossRatio();
            Double_t dAy_cross = GetAnalyzingPowerError_CrossRatio();
            std::cout << "From Cross-Ratio: A_y = " << Ay_cross << " ± " << dAy_cross << std::endl;
            
            if (theoreticalAy != 0) {
                std::cout << "\n--- Comparison to Theory ---" << std::endl;
                std::cout << "Theoretical A_y = " << theoreticalAy << std::endl;
                
                Double_t deviation = Ay_cross - theoreticalAy;
                Double_t deviation_percent = (deviation / theoreticalAy) * 100.0;
                Double_t significance = GetSignificance();
                
                std::cout << "Deviation: " << deviation << " (" << deviation_percent << "%)" << std::endl;
                std::cout << "Statistical significance: " << significance << " σ" << std::endl;
                
                if (significance < 1.0) {
                    std::cout << "  → Excellent agreement with theory!" << std::endl;
                } else if (significance < 2.0) {
                    std::cout << "  → Good agreement with theory" << std::endl;
                } else if (significance < 3.0) {
                    std::cout << "  → Moderate disagreement (~2σ)" << std::endl;
                } else {
                    std::cout << "  → Significant disagreement (>3σ) - check systematics!" << std::endl;
                }
            }
        }
        
        std::cout << "\n--- Selection Cuts ---" << std::endl;
        for (Int_t i = 0; i < nScintillators; i++) {
            std::cout << "S" << i << ": " << cutMin[i] << " - " << cutMax[i] << " MeV" << std::endl;
        }
        std::cout << "========================\n" << std::endl;
    }
    
    void PrintCutInfo() {
        std::cout << "\n=== " << detectorName << " Elastic Selection Cuts ===" << std::endl;
        for (Int_t i = 0; i < nScintillators; i++) {
            std::cout << "S" << i << ": " << cutMin[i] << " - " << cutMax[i] << " MeV" << std::endl;
        }
    }
};

#endif