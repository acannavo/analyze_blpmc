// Kinematics.h
// Frame conversion and phase space mapping for p+C scattering

#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "TF1.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "PhysicalConstants.h"
#include "DetectorConfig.h"

// ====================================================================
// CreatePhiSamplingFunction: Create TF1 for azimuthal angle sampling
// ====================================================================
// The azimuthal distribution for polarized beam scattering is:
//   dN/dφ ∝ 1 + P·A_N·cos(φ)
// 
// Where:
//   - P = beam polarization (magnitude, 0 to 1)
//   - A_N = analyzing power (depends on energy and theta)
//   - φ = azimuthal angle in CM frame
//   - cos(0°) = +1 (right), cos(180°) = -1 (left)
//
// For spin-up (+1):   more scattering to RIGHT (φ=0°)
// For spin-down (-1): more scattering to LEFT (φ=180°)
//
// Inputs:
//   - polarization: beam polarization magnitude (0 to 1)
//   - analyzing_power: A_N value for this scattering angle
//   - spin_state: +1 for up, -1 for down
//
// Returns:
//   - TF1* that can be used with GetRandom() to sample φ
// ====================================================================
TF1* CreatePhiSamplingFunction(Double_t polarization, Double_t analyzing_power, Int_t spin_state);

// ====================================================================
// ConvertPhiCMtoLab: Convert azimuthal angle from CM to LAB frame
// ====================================================================
// Unlike theta, phi DOES change under Lorentz boost because the 
// transverse momentum components (px, py) are affected by the boost.
//
// Inputs:
//   - ekin: beam kinetic energy [MeV]
//   - theta_cm_deg: polar angle in CM frame [degrees]
//   - phi_cm_rad: azimuthal angle in CM frame [radians]
//
// Output:
//   - phi_lab: azimuthal angle in LAB frame [radians]
// ====================================================================
Double_t ConvertPhiCMtoLab(Double_t ekin, Double_t theta_cm_deg, Double_t phi_cm_rad);

// ====================================================================
// ConvertThetaCMtoLab: Convert polar angle from CM to LAB frame
// ====================================================================
// Inputs:
//   - ekin: beam kinetic energy [MeV]
//   - theta_cm_deg: polar angle in CM frame [degrees]
//
// Output:
//   - theta_lab_deg: polar angle in LAB frame [degrees]
// ====================================================================
Double_t ConvertThetaCMtoLab(Double_t ekin, Double_t theta_cm_deg);

// ====================================================================
// CalculateCMMomentum: Calculate CM momentum for elastic scattering
// ====================================================================
// Calculate the momentum of particles in the center-of-mass frame
// for p+C elastic scattering at given beam energy.
//
// Input:  ekin - beam kinetic energy [MeV]
// Output: p_CM - CM momentum [MeV/c]
// ====================================================================
Double_t CalculateCMMomentum(Double_t ekin);

// ====================================================================
// ConvertThetaCMtoT: Convert θ_CM to Mandelstam t
// ====================================================================
// Convert polar angle in CM frame to momentum transfer squared.
// 
// Formula: t = -2p²_CM (1 - cos(θ_CM))
//
// Inputs:
//   - theta_cm_deg: polar angle in CM frame [degrees]
//   - ekin: beam kinetic energy [MeV]
//
// Output:
//   - t: momentum transfer squared [GeV²] (negative for physical scattering)
// ====================================================================
Double_t ConvertThetaCMtoT(Double_t theta_cm_deg, Double_t ekin);

// ====================================================================
// ConvertTtoThetaCM: Convert Mandelstam t to θ_CM
// ====================================================================
// Convert momentum transfer squared to polar angle in CM frame.
//
// Formula: cos(θ_CM) = 1 + t/(2p²_CM)
//
// Inputs:
//   - t: momentum transfer squared [GeV²]
//   - ekin: beam kinetic energy [MeV]
//
// Output:
//   - theta_cm_deg: polar angle in CM frame [degrees]
// ====================================================================
Double_t ConvertTtoThetaCM(Double_t t, Double_t ekin);


// ====================================================================
// ComputeUnifiedTRange: Calculate unified t range for energy interval
// ====================================================================
// Calculate the unified momentum transfer range that covers the detector
// acceptance for all energies in [ekin_min, ekin_max].
//
// For single energy studies, use ekin_min = ekin_max = your energy.
//
// Inputs:
//   - ekin_min, ekin_max: beam energy range [MeV]
//   - theta_lab_center: central lab angle [radians]
//   - theta_lab_window: half-width of acceptance [radians]
//
// Outputs:
//   - t_min_out, t_max_out: unified t range [GeV²]
// ====================================================================
void ComputeUnifiedTRange(Double_t ekin_min, Double_t ekin_max,
                          Double_t theta_lab_center, Double_t theta_lab_window,
                          Double_t& t_min_out, Double_t& t_max_out);


// ====================================================================
// ComputePhiCMRanges: Find CM phi ranges that map to LAB acceptance
// ====================================================================
// The detector accepts events at:
//   - phi_lab ≈ 0° ± 5 mrad
//   - phi_lab ≈ 180° ± 5 mrad
//
// This function finds which CM phi values will boost into these windows.
//
// Inputs:
//   - ekin: beam kinetic energy [MeV]
//   - theta_cm_deg: polar angle in CM (approximately constant in our case)
//
// Outputs:
//   - phi_cm_min_0, phi_cm_max_0: CM range for 0° detector [radians]
//   - phi_cm_min_180, phi_cm_max_180: CM range for 180° detector [radians]
// ====================================================================
void ComputePhiCMRanges(Double_t ekin, Double_t theta_cm_deg,
                        Double_t& phi_cm_min_0, Double_t& phi_cm_max_0,
                        Double_t& phi_cm_min_180, Double_t& phi_cm_max_180);

// ====================================================================
// ComputeCMAngleRange: Find CM theta range that maps to LAB acceptance
// ====================================================================
// Given a lab acceptance window, find the corresponding CM angle range.
//
// Inputs:
//   - ekin: beam kinetic energy [MeV]
//   - theta_lab_center: central lab angle [radians]
//   - theta_lab_width: half-width of acceptance [radians]
//
// Outputs:
//   - theta_cm_min, theta_cm_max: CM angle range [degrees]
// ====================================================================
void ComputeCMAngleRange(Double_t ekin, Double_t theta_lab_center, 
                         Double_t theta_lab_width, 
                         Double_t& theta_cm_min, Double_t& theta_cm_max);

#endif // KINEMATICS_H