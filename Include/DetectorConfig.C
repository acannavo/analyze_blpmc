// ====================================================================
// DetectorConfig.C
// ====================================================================
// Implementation of detector configuration management
//
// This file defines the global detector configuration variables and
// provides functions to modify them at runtime.
//
// Key Features:
// - Dynamic detector positioning (for systematic studies)
// - Automatic recalculation of derived constants
// - Configuration validation and reporting
//
// Dependencies: DetectorConfig.h, TMath.h, iostream
//
// ====================================================================
// IMPORTANT NOTE FOR KINEMATICS MODULE
// ====================================================================
// The function ComputeCMAngleRange() in Kinematics.C uses a lookup
// table that samples θ_CM from 0° to 50° (line 126). This range is
// sufficient for the default detector position (θ = 16.2°) but must
// be extended if detector is repositioned to larger angles.
//
// Required modifications for different detector angles:
//   DETECTOR_THETA_CENTER < 35° → No change needed (default 50° OK)
//   DETECTOR_THETA_CENTER > 40° → Change to 90° range in Kinematics.C
//   DETECTOR_THETA_CENTER > 70° → Change to 180° range in Kinematics.C
//
// To modify: Edit line 126 in Kinematics.C:
//   Double_t theta_cm_deg = i * 50.0 / nPoints;  // Change 50.0
// ====================================================================



#include "DetectorConfig.h"
#include <TMath.h>
#include <iostream>

// ====================================================================
// GLOBAL DETECTOR CONFIGURATION VARIABLES
// ====================================================================
// Default values correspond to nominal experimental setup
// Can be modified using SetDetectorConfig()
// ====================================================================

Double_t DETECTOR_THETA_CENTER = 16.2;      // degrees
Double_t DETECTOR_THETA_WINDOW = 0.005;     // radians (±1.0°)
Double_t DETECTOR_PHI_WINDOW = 0.005;       // radians (±1.0°)
Double_t POLARIMETER_DISTANCE = 2.18;       // meters


// ====================================================================
// DERIVED CONSTANTS
// ====================================================================
// Automatically calculated from primary parameters
// ====================================================================

Double_t DETECTOR_THETA_CENTER_RAD = DETECTOR_THETA_CENTER * TMath::DegToRad();
Double_t DETECTOR_THETA_MIN = DETECTOR_THETA_CENTER_RAD - DETECTOR_THETA_WINDOW;
Double_t DETECTOR_THETA_MAX = DETECTOR_THETA_CENTER_RAD + DETECTOR_THETA_WINDOW;
Double_t DETECTOR_PHI_MIN_0 = -DETECTOR_PHI_WINDOW;
Double_t DETECTOR_PHI_MAX_0 = +DETECTOR_PHI_WINDOW;
Double_t DETECTOR_PHI_MIN_180 = TMath::Pi() - DETECTOR_PHI_WINDOW;
Double_t DETECTOR_PHI_MAX_180 = TMath::Pi() + DETECTOR_PHI_WINDOW;

// ====================================================================
// FUNCTION: SetDetectorConfig
// ====================================================================
// Purpose: Update detector configuration and recalculate derived constants
//
// Parameters:
//   theta_center  - Detector center angle [degrees]
//   theta_window  - Angular acceptance window [radians] (default: 0.070 rad ≈ 4°)
//   phi_window    - Azimuthal acceptance window [radians] (default: 0.017 rad ≈ 1°)
//   distance      - Polarimeter distance from target [meters] (default: 2.18 m)
//
// Usage:
//   SetDetectorConfig(14.0);              // Set to 14.0°, use defaults for others
//   SetDetectorConfig(16.2, 0.05);        // Set to 16.2°, custom theta window
//   SetDetectorConfig(18.0, 0.070, 0.02); // Custom theta and phi windows
//
// Notes:
//   - All derived constants are automatically recalculated
//   - Configuration is printed to screen for verification
//   - No validation checks (assumes user provides sensible values)
// ====================================================================

void SetDetectorConfig(Double_t theta_center,
                       Double_t theta_window,
                       Double_t phi_window,
                       Double_t distance) {
    
    // Update primary parameters
    DETECTOR_THETA_CENTER = theta_center;
    DETECTOR_THETA_WINDOW = theta_window;
    DETECTOR_PHI_WINDOW = phi_window;
    POLARIMETER_DISTANCE = distance;
    
    // Recalculate derived constants
    DETECTOR_THETA_CENTER_RAD = DETECTOR_THETA_CENTER * TMath::DegToRad();
    DETECTOR_THETA_MIN = DETECTOR_THETA_CENTER_RAD - DETECTOR_THETA_WINDOW;
    DETECTOR_THETA_MAX = DETECTOR_THETA_CENTER_RAD + DETECTOR_THETA_WINDOW;
    DETECTOR_PHI_MIN_0 = -DETECTOR_PHI_WINDOW;
    DETECTOR_PHI_MAX_0 = +DETECTOR_PHI_WINDOW;
    DETECTOR_PHI_MIN_180 = TMath::Pi() - DETECTOR_PHI_WINDOW;
    DETECTOR_PHI_MAX_180 = TMath::Pi() + DETECTOR_PHI_WINDOW;
    
    // Print confirmation
    std::cout << "\n==========================================" << std::endl;
    std::cout << "   Detector Configuration Updated" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "  Theta center:  " << DETECTOR_THETA_CENTER << " deg" << std::endl;
    std::cout << "  Theta window:  ±" << DETECTOR_THETA_WINDOW << " rad";
    std::cout << " (±" << DETECTOR_THETA_WINDOW * TMath::RadToDeg() << " deg)" << std::endl;
    std::cout << "  Phi window:    ±" << DETECTOR_PHI_WINDOW << " rad";
    std::cout << " (±" << DETECTOR_PHI_WINDOW * TMath::RadToDeg() << " deg)" << std::endl;
    std::cout << "  Distance:      " << POLARIMETER_DISTANCE << " m" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "  Theta range:   [" << DETECTOR_THETA_MIN * TMath::RadToDeg();
    std::cout << ", " << DETECTOR_THETA_MAX * TMath::RadToDeg() << "] deg" << std::endl;
    std::cout << "==========================================" << std::endl << std::endl;
}

// ====================================================================
// FUNCTION: PrintDetectorConfig
// ====================================================================
// Purpose: Display current detector configuration (detailed output)
//
// Usage:
//   PrintDetectorConfig();
//
// Output: Shows all primary and derived configuration parameters
// ====================================================================

void PrintDetectorConfig() {
    std::cout << "\n==========================================" << std::endl;
    std::cout << "   Current Detector Configuration" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "\nPrimary Parameters:" << std::endl;
    std::cout << "  DETECTOR_THETA_CENTER = " << DETECTOR_THETA_CENTER << " deg" << std::endl;
    std::cout << "  DETECTOR_THETA_WINDOW = " << DETECTOR_THETA_WINDOW << " rad";
    std::cout << " (±" << DETECTOR_THETA_WINDOW * TMath::RadToDeg() << " deg)" << std::endl;
    std::cout << "  DETECTOR_PHI_WINDOW   = " << DETECTOR_PHI_WINDOW << " rad";
    std::cout << " (±" << DETECTOR_PHI_WINDOW * TMath::RadToDeg() << " deg)" << std::endl;
    std::cout << "  POLARIMETER_DISTANCE  = " << POLARIMETER_DISTANCE << " m" << std::endl;
    
    std::cout << "\nDerived Constants:" << std::endl;
    std::cout << "  DETECTOR_THETA_CENTER_RAD = " << DETECTOR_THETA_CENTER_RAD << " rad" << std::endl;
    std::cout << "  DETECTOR_THETA_MIN        = " << DETECTOR_THETA_MIN << " rad";
    std::cout << " (" << DETECTOR_THETA_MIN * TMath::RadToDeg() << " deg)" << std::endl;
    std::cout << "  DETECTOR_THETA_MAX        = " << DETECTOR_THETA_MAX << " rad";
    std::cout << " (" << DETECTOR_THETA_MAX * TMath::RadToDeg() << " deg)" << std::endl;
    std::cout << "  DETECTOR_PHI_MIN_0        = " << DETECTOR_PHI_MIN_0 << " rad" << std::endl;
    std::cout << "  DETECTOR_PHI_MAX_0        = " << DETECTOR_PHI_MAX_0 << " rad" << std::endl;
    std::cout << "  DETECTOR_PHI_MIN_180      = " << DETECTOR_PHI_MIN_180 << " rad" << std::endl;
    std::cout << "  DETECTOR_PHI_MAX_180      = " << DETECTOR_PHI_MAX_180 << " rad" << std::endl;
    
    std::cout << "\nAcceptance Summary:" << std::endl;
    std::cout << "  Theta range: [" << DETECTOR_THETA_MIN * TMath::RadToDeg();
    std::cout << ", " << DETECTOR_THETA_MAX * TMath::RadToDeg() << "] deg" << std::endl;
    std::cout << "  Phi @ 0°:    [" << DETECTOR_PHI_MIN_0 * TMath::RadToDeg();
    std::cout << ", " << DETECTOR_PHI_MAX_0 * TMath::RadToDeg() << "] deg" << std::endl;
    std::cout << "  Phi @ 180°:  [" << DETECTOR_PHI_MIN_180 * TMath::RadToDeg();
    std::cout << ", " << DETECTOR_PHI_MAX_180 * TMath::RadToDeg() << "] deg" << std::endl;
    std::cout << "==========================================" << std::endl << std::endl;
}

// ====================================================================
// FUNCTION: ResetDetectorConfig
// ====================================================================
// Purpose: Reset detector configuration to default values
//
// Usage:
//   ResetDetectorConfig();
//
// Notes: Default configuration is 16.2° center with standard windows
// ====================================================================

void ResetDetectorConfig() {
    std::cout << "\nResetting detector configuration to defaults..." << std::endl;
    SetDetectorConfig(16.2, 0.005, 0.005, 2.18);
}