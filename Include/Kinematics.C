// Kinematics.C
// Implementation of frame conversion and phase space mapping

#include "Kinematics.h"
#include <iostream>

using namespace std;

// ====================================================================
// CreatePhiSamplingFunction: Create TF1 for azimuthal angle sampling
// ====================================================================
TF1* CreatePhiSamplingFunction(Double_t polarization, Double_t analyzing_power, Int_t spin_state)
{
    // Create function: 1 + P·A_N·cos(φ)
    // Parameters: [0] = P·A_N (effective asymmetry)
    TF1 *fPhi = new TF1("fPhi", "1 + [0]*cos(x)", 0, 2*TMath::Pi());
    
    // Set the effective asymmetry parameter
    // For spin-up: P·A_N is positive → more at φ=0
    // For spin-down: -P·A_N is negative → more at φ=180
    Double_t effective_asymmetry = -spin_state * polarization * analyzing_power;
    fPhi->SetParameter(0, effective_asymmetry);
    
    return fPhi;
}

// ====================================================================
// ConvertPhiCMtoLab: Convert azimuthal angle from CM to LAB frame
// ====================================================================
Double_t ConvertPhiCMtoLab(Double_t ekin, Double_t theta_cm_deg, Double_t phi_cm_rad)
{
    // Calculate initial state
    Double_t ekinGeV = ekin/1000.;
    Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
    
    // Calculate CM momentum
    Double_t s = iState.Mag2();
    Double_t pcm = sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));
    
    // Create proton 4-momentum in CM frame
    Double_t theta_cm_rad = theta_cm_deg * TMath::DegToRad();
    TLorentzVector p_cm;
    p_cm.SetXYZM(
        pcm * TMath::Sin(theta_cm_rad) * TMath::Cos(phi_cm_rad),
        pcm * TMath::Sin(theta_cm_rad) * TMath::Sin(phi_cm_rad),
        pcm * TMath::Cos(theta_cm_rad),
        mp
    );
    
    // Boost to LAB frame
    p_cm.Boost(iState.BoostVector());
    
    // Get LAB phi
    Double_t phi_lab = p_cm.Phi();  // Returns value in [-π, +π]
    
    return phi_lab;
}

// ====================================================================
// ConvertThetaCMtoLab: Convert polar angle from CM to LAB frame
// ====================================================================
Double_t ConvertThetaCMtoLab(Double_t ekin, Double_t theta_cm_deg)
{
    // Calculate initial state
    Double_t ekinGeV = ekin/1000.;
    Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
    
    // Calculate CM momentum
    Double_t s = iState.Mag2();
    Double_t pcm = sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));
    
    // Create proton 4-momentum in CM frame
    // Use phi=0 since we only care about theta conversion
    Double_t theta_cm_rad = theta_cm_deg * TMath::DegToRad();
    TLorentzVector p_cm;
    p_cm.SetXYZM(
        pcm * TMath::Sin(theta_cm_rad),
        0.,
        pcm * TMath::Cos(theta_cm_rad),
        mp
    );
    
    // Boost to LAB frame
    p_cm.Boost(iState.BoostVector());
    
    // Get LAB theta
    Double_t theta_lab_rad = p_cm.Theta();
    Double_t theta_lab_deg = theta_lab_rad * TMath::RadToDeg();
    
    return theta_lab_deg;
}

// ====================================================================
// CalculateCMMomentum: Calculate CM momentum for elastic scattering
// ====================================================================
// For p+C elastic scattering, calculates the momentum of particles
// in the center-of-mass frame using Mandelstam s.
//
// Input:  ekin - beam kinetic energy [MeV]
// Output: p_CM - CM momentum [MeV/c]
// ====================================================================
Double_t CalculateCMMomentum(Double_t ekin)
{
    // Convert to GeV
    Double_t ekinGeV = ekin / 1000.;
    
    // Calculate beam momentum
    Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    
    // Initial state 4-momentum: beam proton + target carbon at rest
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
    
    // Mandelstam s (total energy squared in CM)
    Double_t s = iState.Mag2();
    
    // CM momentum formula for elastic scattering
    Double_t pcm = TMath::Sqrt((s-(mp+mC)*(mp+mC)) * 
                               (s-(mp-mC)*(mp-mC)) / (4.*s));
    
    // Convert to MeV/c
    return pcm * 1000.0;
}

// ====================================================================
// ConvertThetaCMtoT: Convert θ_CM to Mandelstam t
// ====================================================================
// Converts polar angle in CM frame to momentum transfer squared.
// 
// Formula: t = -2p²_CM (1 - cos(θ_CM))
//
// Inputs:
//   - theta_cm_deg: polar angle in CM frame [degrees]
//   - ekin: beam kinetic energy [MeV]
//
// Output:
//   - t: momentum transfer squared [GeV²]
//
// Note: t is negative for physical scattering (timelike transfer)
// ====================================================================
Double_t ConvertThetaCMtoT(Double_t theta_cm_deg, Double_t ekin)
{
    // Calculate CM momentum in GeV/c
    Double_t pcm_GeV = CalculateCMMomentum(ekin) / 1000.;
    
    // Convert angle to radians
    Double_t theta_rad = theta_cm_deg * TMath::DegToRad();
    
    // Calculate t = -2p² (1 - cos(θ))
    Double_t cos_theta = TMath::Cos(theta_rad);
    Double_t t = -2.0 * pcm_GeV * pcm_GeV * (1.0 - cos_theta);
    
    return t;  // GeV²
}

// ====================================================================
// ConvertTtoThetaCM: Convert Mandelstam t to θ_CM
// ====================================================================
// Converts momentum transfer squared to polar angle in CM frame.
//
// Formula: cos(θ_CM) = 1 + t/(2p²_CM)
//
// Inputs:
//   - t: momentum transfer squared [GeV²]
//   - ekin: beam kinetic energy [MeV]
//
// Output:
//   - theta_cm_deg: polar angle in CM frame [degrees]
//
// Note: Handles numerical errors by clamping cos(θ) to [-1, 1]
// ====================================================================
Double_t ConvertTtoThetaCM(Double_t t, Double_t ekin)
{
    // Calculate CM momentum in GeV/c
    Double_t pcm_GeV = CalculateCMMomentum(ekin) / 1000.;
    
    // Invert the formula: cos(θ) = 1 + t/(2p²)
    Double_t cos_theta = 1.0 + t / (2.0 * pcm_GeV * pcm_GeV);
    
    // Handle numerical errors (ensure cos_theta in [-1, 1])
    if (cos_theta > 1.0)  cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    
    // Convert to degrees
    Double_t theta_rad = TMath::ACos(cos_theta);
    Double_t theta_deg = theta_rad * TMath::RadToDeg();
    
    return theta_deg;
}

// ====================================================================
// ComputeUnifiedTRange: Calculate unified t range for energy interval
// ====================================================================
void ComputeUnifiedTRange(Double_t ekin_min, Double_t ekin_max,
                          Double_t theta_lab_center, Double_t theta_lab_window,
                          Double_t& t_min_out, Double_t& t_max_out)
{
    // Calculate CM angle range for minimum energy
    Double_t theta_cm_min_low, theta_cm_max_low;
    ComputeCMAngleRange(ekin_min, theta_lab_center, theta_lab_window,
                        theta_cm_min_low, theta_cm_max_low);
    
    // Convert to t range (note: larger theta → more negative t)
    Double_t t_min_low = ConvertThetaCMtoT(theta_cm_max_low, ekin_min);
    Double_t t_max_low = ConvertThetaCMtoT(theta_cm_min_low, ekin_min);
    
    // Calculate CM angle range for maximum energy
    Double_t theta_cm_min_high, theta_cm_max_high;
    ComputeCMAngleRange(ekin_max, theta_lab_center, theta_lab_window,
                        theta_cm_min_high, theta_cm_max_high);
    
    // Convert to t range
    Double_t t_min_high = ConvertThetaCMtoT(theta_cm_max_high, ekin_max);
    Double_t t_max_high = ConvertThetaCMtoT(theta_cm_min_high, ekin_max);
    
    // Take union of ranges (most negative to least negative)
    t_min_out = TMath::Min(t_min_low, t_min_high);
    t_max_out = TMath::Max(t_max_low, t_max_high);
    
    // Print summary
    cout << "\n=== Unified t-Range Calculation ===" << endl;
    cout << "Energy range: [" << ekin_min << ", " << ekin_max << "] MeV" << endl;
    cout << "Lab acceptance: " << theta_lab_center*TMath::RadToDeg() 
         << " ± " << theta_lab_window*TMath::RadToDeg() << " deg" << endl;
    cout << "Unified t range: [" << t_min_out << ", " << t_max_out << "] GeV²" << endl;
    cout << "Δt = " << (t_max_out - t_min_out) << " GeV²" << endl;
    cout << "====================================\n" << endl;
}

// ====================================================================
// ComputePhiCMRanges: Find CM phi ranges that map to LAB acceptance
// ====================================================================
void ComputePhiCMRanges(Double_t ekin, Double_t theta_cm_deg,
                        Double_t& phi_cm_min_0, Double_t& phi_cm_max_0,
                        Double_t& phi_cm_min_180, Double_t& phi_cm_max_180)
{
    // Lab acceptance windows
    const Double_t phi_lab_window = DETECTOR_PHI_WINDOW;  // 5 mrad
    
    // Window 1: around 0°
    const Double_t phi_lab_min_0 = -phi_lab_window;
    const Double_t phi_lab_max_0 = phi_lab_window;
    
    // Window 2: around 180° (near ±π)
    const Double_t phi_lab_min_180 = TMath::Pi() - phi_lab_window;
    const Double_t phi_lab_max_180 = -TMath::Pi() + phi_lab_window;
    
    // Create lookup table: phi_CM → phi_LAB
    const Int_t nPoints = 100000;  // Increased for better precision
    Double_t phi_cm_values[nPoints];
    Double_t phi_lab_values[nPoints];
    
    for (Int_t i = 0; i < nPoints; i++)
    {
        // Sample phi_CM from -π to +π
        Double_t phi_cm = -TMath::Pi() + i * 2.0 * TMath::Pi() / nPoints;
        phi_cm_values[i] = phi_cm;
        phi_lab_values[i] = ConvertPhiCMtoLab(ekin, theta_cm_deg, phi_cm);
    }
    
    // Find CM ranges for window 1 (around 0°)
    phi_cm_min_0 = -TMath::Pi();
    phi_cm_max_0 = TMath::Pi();
    Bool_t found_min_0 = false;
    Bool_t found_max_0 = false;
    
    for (Int_t i = 0; i < nPoints; i++)
    {
        if (phi_lab_values[i] >= phi_lab_min_0 && phi_lab_values[i] <= phi_lab_max_0)
        {
            if (!found_min_0) {
                phi_cm_min_0 = phi_cm_values[i];
                found_min_0 = true;
            }
            phi_cm_max_0 = phi_cm_values[i];
            found_max_0 = true;
        }
    }
    
    // Find CM ranges for window 2 (around 180°)
    // The acceptance is ±5 mrad around 180°, which appears as two regions:
    // near +π: [π - 5mrad, π]
    // near -π: [-π, -π + 5mrad]
    // Due to symmetry, phi_CM also splits into two regions
    
    // Since phi is nearly unchanged CM→LAB, we expect:
    // phi_CM_min ≈ π - 5.5mrad, phi_CM_max ≈ π
    // We'll report this as a symmetric range around ±π
    
    Double_t phi_cm_threshold_pos = TMath::Pi() - 2*phi_lab_window;  // Start looking from here
    
    // Find the smallest phi_CM that maps to the 180° window
    phi_cm_min_180 = TMath::Pi();
    phi_cm_max_180 = TMath::Pi();
    
    for (Int_t i = 0; i < nPoints; i++)
    {
        Bool_t near_plus_pi = (phi_lab_values[i] >= phi_lab_min_180);
        Bool_t near_minus_pi = (phi_lab_values[i] <= phi_lab_max_180);
        
        if ((near_plus_pi || near_minus_pi) && phi_cm_values[i] > 0) {
            // Only look at positive phi_CM values (near +π)
            if (phi_cm_values[i] < phi_cm_min_180) {
                phi_cm_min_180 = phi_cm_values[i];
            }
        }
    }
    
    // phi_cm_max_180 is just π
    phi_cm_max_180 = TMath::Pi();
    
    // By symmetry, the range near -π is the mirror:
    // [-π, -phi_cm_min_180]
    
    // Add small margin (10% of the acceptance window)
    Double_t margin = 0.1 * phi_lab_window;  // 10% of 5 mrad = 0.5 mrad
    phi_cm_min_0 -= margin;
    phi_cm_max_0 += margin;
    phi_cm_min_180 -= margin;
    // Don't add margin to phi_cm_max_180 since it's already at π
    // The margin is already included in phi_cm_min_180
    
    // Clamp detector 1 range only (detector 2 is already at boundaries)
    phi_cm_min_0 = TMath::Max(-TMath::Pi(), phi_cm_min_0);
    phi_cm_max_0 = TMath::Min(TMath::Pi(), phi_cm_max_0);
}

// ====================================================================
// ComputeCMAngleRange: Find CM theta range that maps to LAB acceptance
// ====================================================================
// ====================================================================
// IMPORTANT: CM ANGLE RANGE LIMITATION
// ====================================================================
// The lookup table samples θ_CM from 0° to 50°, which maps to 
// θ_Lab from 0° to ~40° for p+C at 200 MeV. This is sufficient for
// forward detector angles.
//
// If detector is positioned at larger Lab angles, extend this range:
//   - For θ_Lab up to ~74°:  Change 50.0 → 90.0
//   - For θ_Lab up to ~180°: Change 50.0 → 180.0
//
// Current range chosen for detector at θ_Lab ≈ 16.2° (default)
// ====================================================================
//for (Int_t i = 0; i < nPoints; i++)
//{
//   // Sample CM angles from 0° to 50° (sufficient for forward scattering)
//    Double_t theta_cm_deg = i * 50.0 / nPoints;  // ← MODIFY HERE IF NEEDED


void ComputeCMAngleRange(Double_t ekin, Double_t theta_lab_center, 
                         Double_t theta_lab_width, 
                         Double_t& theta_cm_min, Double_t& theta_cm_max)
{
    // Define lab acceptance window
    Double_t theta_lab_min_target = theta_lab_center - theta_lab_width;
    Double_t theta_lab_max_target = theta_lab_center + theta_lab_width;
    
    // Create lookup table: CM angle → Lab angle using our conversion function
    const Int_t nPoints = 10000;
    Double_t cm_angles[nPoints];
    Double_t lab_angles[nPoints];
    
    for (Int_t i = 0; i < nPoints; i++)
    {
        // Sample CM angles from 0° to 50° (sufficient for forward scattering)
        Double_t theta_cm_deg = i * 50.0 / nPoints;
        cm_angles[i] = theta_cm_deg;
        
        // Use our conversion function
        lab_angles[i] = ConvertThetaCMtoLab(ekin, theta_cm_deg);
    }
    
    // Find CM angles corresponding to lab acceptance edges
    theta_cm_min = 0.;
    theta_cm_max = 50.;
    Bool_t found_min = false;
    Bool_t found_max = false;
    
    for (Int_t i = 0; i < nPoints; i++)
    {
        // Find first CM angle that maps to minimum lab angle
        if (!found_min && lab_angles[i] >= theta_lab_min_target * TMath::RadToDeg()) {
            theta_cm_min = cm_angles[i];
            found_min = true;
        }
        // Find last CM angle that maps to maximum lab angle
        if (lab_angles[i] <= theta_lab_max_target * TMath::RadToDeg()) {
            theta_cm_max = cm_angles[i];
            found_max = true;
        }
    }
    
    // Add 10% margin to ensure we cover the full acceptance
    Double_t margin = 0.1 * (theta_cm_max - theta_cm_min);
    theta_cm_min = TMath::Max(0., theta_cm_min - margin);
    theta_cm_max = theta_cm_max + margin;
    
    // Print summary
    cout << "\n=== CM Angle Range Calculation ===" << endl;
    cout << "Beam energy: " << ekin << " MeV" << endl;
    cout << "Lab acceptance: " << theta_lab_center*TMath::RadToDeg() 
         << " ± " << theta_lab_width*TMath::RadToDeg() << " deg" << endl;
    cout << "               " << theta_lab_min_target*TMath::RadToDeg()
         << " to " << theta_lab_max_target*TMath::RadToDeg() << " deg" << endl;
    cout << "CM angle range: [" << theta_cm_min << ", " 
         << theta_cm_max << "] degrees" << endl;
    cout << "===================================\n" << endl;
}