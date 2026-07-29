// ====================================================================
// DetectorConfig.h
// ====================================================================
// Detector geometry and acceptance parameters
// 
// This header declares global variables for detector configuration.
// The actual values are defined in DetectorConfig.C and can be 
// modified at runtime using SetDetectorConfig().
//
// Usage:
//   SetDetectorConfig(16.2);  // Set detector to 16.4 degrees
//   PrintDetectorConfig();     // Display current configuration
//
// Dependencies: TMath.h
// Used by: All modules requiring detector acceptance cuts
// ====================================================================

#ifndef DETECTOR_CONFIG_H
#define DETECTOR_CONFIG_H

#include <TMath.h>

// ====================================================================
// PRIMARY DETECTOR PARAMETERS
// ====================================================================
// These can be modified using SetDetectorConfig() in DetectorConfig.C
// ====================================================================

extern Double_t DETECTOR_THETA_CENTER;    // Center angle [degrees]
extern Double_t DETECTOR_THETA_WINDOW;    // Angular acceptance [radians]
extern Double_t DETECTOR_PHI_WINDOW;      // Azimuthal acceptance [radians]
extern Double_t POLARIMETER_DISTANCE;     // Distance from target [meters]

// ====================================================================
// DERIVED CONSTANTS
// ====================================================================
// These are automatically recalculated when SetDetectorConfig() is called
// ====================================================================

extern Double_t DETECTOR_THETA_CENTER_RAD;  // Center angle [radians]
extern Double_t DETECTOR_THETA_MIN;         // Min acceptance angle [radians]
extern Double_t DETECTOR_THETA_MAX;         // Max acceptance angle [radians]
extern Double_t DETECTOR_PHI_MIN_0;         // Min phi at 0° [radians]
extern Double_t DETECTOR_PHI_MAX_0;         // Max phi at 0° [radians]
extern Double_t DETECTOR_PHI_MIN_180;       // Min phi at 180° [radians]
extern Double_t DETECTOR_PHI_MAX_180;       // Max phi at 180° [radians]

// ====================================================================
// FUNCTION DECLARATIONS
// ====================================================================

// Set detector configuration and recalculate derived constants
void SetDetectorConfig(Double_t theta_center,           // degrees
                       Double_t theta_window = 0.005,   // radians (default ±4°)
                       Double_t phi_window = 0.005,     // radians (default ±1°)
                       Double_t distance = 2.18);       // meters

// Display current detector configuration
void PrintDetectorConfig();

// Reset to default configuration (theta = 16.2°)
void ResetDetectorConfig();

#endif // DETECTOR_CONFIG_H