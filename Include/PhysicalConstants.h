// PhysicalConstants.h
// Physical constants for p+C scattering

#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

// ====================================================================
// Physical constants
// ====================================================================

// Proton mass in GeV/c²
// Pluto database value (PParticle("p").M()): 0.93827231
// Our value (higher precision):
const Double_t mp = 0.9382720813;  

// Carbon-12 mass in GeV/c² (not available in Pluto database)
const Double_t mC = 11.174862;

// Excitation energy for C-12 2+ state
const Double_t E_EXCITATION = 0.00443;  // 4.43 MeV in GeV

#endif // PHYSICAL_CONSTANTS_H