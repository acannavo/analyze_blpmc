// MeyerScattering.C
// Implementation of Meyer parameterization with energy interpolation

#include "MeyerScattering.h"
#include <iostream>

using namespace std;

// ====================================================================
// MeyerXS_Elastic: Energy-dependent elastic cross section
// ====================================================================
Double_t MeyerXS_Elastic(Double_t t, Double_t ekin)
{
    // CRITICAL: Convert t from GeV² to |t| in (MeV/c)²
    // Input t is NEGATIVE (e.g., -0.05 GeV²)
    // Splines expect POSITIVE |t| in (MeV/c)²
    Double_t t_MeV2 = TMath::Abs(t) * 1.0e6;  // |t| in (MeV/c)²
    
    // Energy-dependent interpolation
    if (ekin <= 160.0) {
        // Use 160 MeV spline directly
        return XSSpline_160MeV_elastic(t_MeV2);
    }
    else if (ekin >= 200.0) {
        // Use 200 MeV spline directly
        return XSSpline_200MeV_elastic(t_MeV2);
    }
    else {
        // Linear interpolation between 160 and 200 MeV
        Double_t xs_160 = XSSpline_160MeV_elastic(t_MeV2);
        Double_t xs_200 = XSSpline_200MeV_elastic(t_MeV2);
        
        // Weight: 0 at 160 MeV, 1 at 200 MeV
        Double_t weight = (ekin - 160.0) / 40.0;
        
        return xs_160 + weight * (xs_200 - xs_160);
    }
}

// ====================================================================
// MeyerXS_Inelastic: Energy-dependent inelastic cross section
// ====================================================================
Double_t MeyerXS_Inelastic(Double_t t, Double_t ekin)
{
    // Convert t from GeV² to |t| in (MeV/c)²
    Double_t t_MeV2 = TMath::Abs(t) * 1.0e6;
    
    // Energy-dependent interpolation
    if (ekin <= 160.0) {
        return XSSpline_160MeV_inelastic(t_MeV2);
    }
    else if (ekin >= 200.0) {
        return XSSpline_200MeV_inelastic(t_MeV2);
    }
    else {
        // Linear interpolation
        Double_t xs_160 = XSSpline_160MeV_inelastic(t_MeV2);
        Double_t xs_200 = XSSpline_200MeV_inelastic(t_MeV2);
        Double_t weight = (ekin - 160.0) / 40.0;
        
        return xs_160 + weight * (xs_200 - xs_160);
    }
}

// ====================================================================
// MeyerAP_Elastic: Energy-dependent elastic analyzing power
// ====================================================================
Double_t MeyerAP_Elastic(Double_t t, Double_t ekin)
{
    // Convert t from GeV² to |t| in (MeV/c)²
    Double_t t_MeV2 = TMath::Abs(t) * 1.0e6;
    
    // Energy-dependent interpolation
    if (ekin <= 160.0) {
        return APSpline_160MeV_elastic(t_MeV2);
    }
    else if (ekin >= 200.0) {
        return APSpline_200MeV_elastic(t_MeV2);
    }
    else {
        // Linear interpolation
        Double_t ap_160 = APSpline_160MeV_elastic(t_MeV2);
        Double_t ap_200 = APSpline_200MeV_elastic(t_MeV2);
        Double_t weight = (ekin - 160.0) / 40.0;
        
        return ap_160 + weight * (ap_200 - ap_160);
    }
}

// ====================================================================
// MeyerAP_Inelastic: Energy-dependent inelastic analyzing power
// ====================================================================
Double_t MeyerAP_Inelastic(Double_t t, Double_t ekin)
{
    // Convert t from GeV² to |t| in (MeV/c)²
    Double_t t_MeV2 = TMath::Abs(t) * 1.0e6;
    
    // Energy-dependent interpolation
    if (ekin <= 160.0) {
        return APSpline_160MeV_inelastic(t_MeV2);
    }
    else if (ekin >= 200.0) {
        return APSpline_200MeV_inelastic(t_MeV2);
    }
    else {
        // Linear interpolation
        Double_t ap_160 = APSpline_160MeV_inelastic(t_MeV2);
        Double_t ap_200 = APSpline_200MeV_inelastic(t_MeV2);
        Double_t weight = (ekin - 160.0) / 40.0;
        
        return ap_160 + weight * (ap_200 - ap_160);
    }
}

// ====================================================================
// CreateMeyerEnvelope_Elastic: Create envelope for elastic scattering
// ====================================================================
TH1D* CreateMeyerEnvelope_Elastic(Double_t t_min, Double_t t_max, Int_t nbins)
{
    // Create histogram spanning the t range
    TH1D* envelope = new TH1D("meyer_envelope_elastic", 
                              "Meyer Envelope (Elastic)", 
                              nbins, t_min, t_max);
    
    // Fill with MAX(XS_160MeV, XS_200MeV) for each t
    for (Int_t ibin = 1; ibin <= nbins; ibin++) {
        Double_t t_GeV2 = envelope->GetBinCenter(ibin);
        Double_t t_MeV2 = TMath::Abs(t_GeV2) * 1.0e6;
        
        // Get cross sections at both energies
        Double_t xs_160 = XSSpline_160MeV_elastic(t_MeV2);
        Double_t xs_200 = XSSpline_200MeV_elastic(t_MeV2);
        
        // Store maximum (envelope)
        Double_t xs_envelope = TMath::Max(xs_160, xs_200);
        envelope->SetBinContent(ibin, xs_envelope);
    }
    
    return envelope;
}

// ====================================================================
// CreateMeyerEnvelope_Inelastic: Create envelope for inelastic
// ====================================================================
TH1D* CreateMeyerEnvelope_Inelastic(Double_t t_min, Double_t t_max, Int_t nbins)
{
    TH1D* envelope = new TH1D("meyer_envelope_inelastic", 
                              "Meyer Envelope (Inelastic)", 
                              nbins, t_min, t_max);
    
    for (Int_t ibin = 1; ibin <= nbins; ibin++) {
        Double_t t_GeV2 = envelope->GetBinCenter(ibin);
        Double_t t_MeV2 = TMath::Abs(t_GeV2) * 1.0e6;
        
        Double_t xs_160 = XSSpline_160MeV_inelastic(t_MeV2);
        Double_t xs_200 = XSSpline_200MeV_inelastic(t_MeV2);
        
        Double_t xs_envelope = TMath::Max(xs_160, xs_200);
        envelope->SetBinContent(ibin, xs_envelope);
    }
    
    return envelope;
}

// ====================================================================
// SampleT_Meyer_Elastic: Accept/reject sampling for elastic
// ====================================================================
Double_t SampleT_Meyer_Elastic(Double_t ekin, TH1D* envelope)
{
    // Thread-safe random number generator
    // Note: For multithreading, each thread should have its own RNG
    // This is a simplified version - production code should use thread-local RNG
    static TRandom3 rng(0);  // 0 = time-based seed
    
    while (true) {
        // 1. Sample t from envelope histogram
        Double_t t_GeV2 = envelope->GetRandom();
        
        // 2. Calculate actual cross section at this energy
        Double_t xs_actual = MeyerXS_Elastic(t_GeV2, ekin);
        
        // 3. Get envelope value at this t
        Int_t bin = envelope->FindBin(t_GeV2);
        Double_t xs_envelope = envelope->GetBinContent(bin);
        
        // 4. Calculate acceptance ratio
        Double_t ratio = xs_actual / xs_envelope;
        
        // Sanity check: ratio should be ≤ 1
        // With 50k bins, ratio should stay very close to 1.0
        if (ratio > 1.0001) {  // Strict threshold to verify binning works
            cerr << "WARNING: MeyerScattering ratio > 1! " 
                 << "ratio = " << ratio 
                 << ", t = " << t_GeV2 
                 << ", E = " << ekin << endl;
        }
        
        // 5. Accept/reject
        if (rng.Rndm() < ratio) {
            return t_GeV2;  // ACCEPTED
        }
        // else: REJECTED, loop continues
    }
}

// ====================================================================
// SampleT_Meyer_Inelastic: Accept/reject sampling for inelastic
// ====================================================================
Double_t SampleT_Meyer_Inelastic(Double_t ekin, TH1D* envelope)
{
    static TRandom3 rng(0);
    
    while (true) {
        Double_t t_GeV2 = envelope->GetRandom();
        Double_t xs_actual = MeyerXS_Inelastic(t_GeV2, ekin);
        
        Int_t bin = envelope->FindBin(t_GeV2);
        Double_t xs_envelope = envelope->GetBinContent(bin);
        
        Double_t ratio = xs_actual / xs_envelope;
        
        if (ratio > 1.0001) {
            cerr << "WARNING: MeyerScattering (inelastic) ratio > 1! " 
                 << "ratio = " << ratio << endl;
        }
        
        if (rng.Rndm() < ratio) {
            return t_GeV2;
        }
    }
}

// ====================================================================
// MeyerScattering_Info - Test function
// ====================================================================
void MeyerScattering_Info()
{
    cout << "\n╔════════════════════════════════════════════════════════╗" << endl;
    cout << "║  Meyer Scattering Module                               ║" << endl;
    cout << "╚════════════════════════════════════════════════════════╝" << endl;
    cout << "Status: Module loaded and IMPLEMENTED" << endl;
    cout << "\nAvailable functions:" << endl;
    cout << "  - MeyerXS_Elastic(t, ekin)   : Elastic cross section" << endl;
    cout << "  - MeyerAP_Elastic(t, ekin)   : Elastic analyzing power" << endl;
    cout << "  - MeyerXS_Inelastic(t, ekin) : Inelastic cross section" << endl;
    cout << "  - MeyerAP_Inelastic(t, ekin) : Inelastic analyzing power" << endl;
    cout << "\nEnergy range: 160-200 MeV (linear interpolation)" << endl;
    cout << "Input: t in GeV² (negative), ekin in MeV" << endl;
    cout << "\nTest example:" << endl;
    cout << "  root [0] MeyerXS_Elastic(-0.05, 180.0)" << endl;
    
    // Quick self-test
    Double_t t_test = -0.05;  // GeV²
    Double_t ekin_test = 180.0;  // MeV
    Double_t xs = MeyerXS_Elastic(t_test, ekin_test);
    Double_t ap = MeyerAP_Elastic(t_test, ekin_test);
    
    cout << "\nSelf-test at t = " << t_test << " GeV², E = " << ekin_test << " MeV:" << endl;
    cout << "  XS = " << xs << " mb/GeV²" << endl;
    cout << "  AP = " << ap << endl;
    cout << endl;
}