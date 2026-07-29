// MeyerScattering.h
// Energy-dependent t-based cross sections and analyzing powers
// Uses Meyer parameterization with splines at specific energies

#ifndef MEYER_SCATTERING_H
#define MEYER_SCATTERING_H

#include "TMath.h"
#include "TH1D.h"
#include "TRandom3.h"

// ====================================================================
// Meyer Spline Files
// ====================================================================
// These files contain TSpline3 definitions for cross sections (XS)
// and analyzing powers (AP) at 160 MeV and 200 MeV.
//
// All splines use |t| in (MeV/c)² as input (positive values)
// ====================================================================

// Cross section splines
#include "XSSpline_160MeV_elastic.C"
#include "XSSpline_200MeV_elastic.C"
#include "XSSpline_160MeV_inelastic.C"
#include "XSSpline_200MeV_inelastic.C"

// Analyzing power splines
#include "APSpline_160MeV_elastic.C"
#include "APSpline_200MeV_elastic.C"
#include "APSpline_160MeV_inelastic.C"
#include "APSpline_200MeV_inelastic.C"

// ====================================================================
// Meyer Cross Section Functions
// ====================================================================

// ====================================================================
// MeyerXS_Elastic: Energy-dependent elastic cross section
// ====================================================================
// Returns differential cross section for elastic scattering using
// Meyer parameterization with linear interpolation between 160-200 MeV.
//
// Inputs:
//   - t: momentum transfer squared [GeV²] (negative for physical scattering)
//   - ekin: beam kinetic energy [MeV]
//
// Output:
//   - dσ/dt [mb/GeV²]
//
// Energy handling:
//   - ekin ≤ 160 MeV: use 160 MeV spline
//   - 160 < ekin < 200 MeV: linear interpolation
//   - ekin ≥ 200 MeV: use 200 MeV spline
// ====================================================================
Double_t MeyerXS_Elastic(Double_t t, Double_t ekin);

// ====================================================================
// MeyerXS_Inelastic: Energy-dependent inelastic cross section
// ====================================================================
// Returns differential cross section for inelastic scattering (4.43 MeV)
// using Meyer parameterization with linear interpolation.
//
// Inputs:
//   - t: momentum transfer squared [GeV²] (negative)
//   - ekin: beam kinetic energy [MeV]
//
// Output:
//   - dσ/dt [mb/GeV²]
// ====================================================================
Double_t MeyerXS_Inelastic(Double_t t, Double_t ekin);

// ====================================================================
// Meyer Analyzing Power Functions
// ====================================================================

// ====================================================================
// MeyerAP_Elastic: Energy-dependent elastic analyzing power
// ====================================================================
// Returns analyzing power A_N for elastic scattering using Meyer
// parameterization with linear interpolation between 160-200 MeV.
//
// Inputs:
//   - t: momentum transfer squared [GeV²] (negative)
//   - ekin: beam kinetic energy [MeV]
//
// Output:
//   - A_N (dimensionless, typically -0.5 to +0.5)
// ====================================================================
Double_t MeyerAP_Elastic(Double_t t, Double_t ekin);

// ====================================================================
// MeyerAP_Inelastic: Energy-dependent inelastic analyzing power
// ====================================================================
// Returns analyzing power A_N for inelastic scattering (4.43 MeV)
// using Meyer parameterization with linear interpolation.
//
// Inputs:
//   - t: momentum transfer squared [GeV²] (negative)
//   - ekin: beam kinetic energy [MeV]
//
// Output:
//   - A_N (dimensionless)
// ====================================================================
Double_t MeyerAP_Inelastic(Double_t t, Double_t ekin);

// ====================================================================
// Utility Functions
// ====================================================================

// ====================================================================
// CreateMeyerEnvelope_Elastic: Create envelope histogram for elastic
// ====================================================================
// Creates a histogram containing MAX(XS_160MeV, XS_200MeV) for each t.
// This envelope is used for accept/reject sampling at any energy.
//
// Inputs:
//   - t_min, t_max: momentum transfer range [GeV²]
//   - nbins: number of histogram bins (default 10000)
//
// Output:
//   - TH1D* containing envelope (caller must delete!)
//
// Usage: Call once per thread at initialization
// ====================================================================
TH1D* CreateMeyerEnvelope_Elastic(Double_t t_min, Double_t t_max, Int_t nbins = 50000);

// ====================================================================
// CreateMeyerEnvelope_Inelastic: Create envelope for inelastic
// ====================================================================
TH1D* CreateMeyerEnvelope_Inelastic(Double_t t_min, Double_t t_max, Int_t nbins = 50000);

// ====================================================================
// SampleT_Meyer_Elastic: Sample t using envelope accept/reject
// ====================================================================
// Samples momentum transfer t from the interpolated cross section
// at the specified energy using accept/reject method.
//
// Inputs:
//   - ekin: beam kinetic energy [MeV]
//   - envelope: envelope histogram (from CreateMeyerEnvelope_Elastic)
//
// Output:
//   - t: sampled momentum transfer [GeV²]
//
// Method:
//   1. Sample t from envelope histogram
//   2. Calculate xs_actual = MeyerXS_Elastic(t, ekin)
//   3. Calculate xs_envelope from histogram
//   4. Accept with probability = xs_actual / xs_envelope
//   5. Repeat until accepted
//
// Usage: Call once per event
// ====================================================================
Double_t SampleT_Meyer_Elastic(Double_t ekin, TH1D* envelope);

// ====================================================================
// SampleT_Meyer_Inelastic: Sample t for inelastic scattering
// ====================================================================
Double_t SampleT_Meyer_Inelastic(Double_t ekin, TH1D* envelope);

// Test function for compilation verification
void MeyerScattering_Info();

#endif // MEYER_SCATTERING_H