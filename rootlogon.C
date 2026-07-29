{
    cout << "=== Custom ROOT Environment Loading ===" << endl;

    // ── Multithreading ──────────────────────────────────────────────
    // Enables ROOT's implicit multithreading globally. RDataFrame
    // (used by PolAnalysis::CountEvents for scintillator-cut counting)
    // automatically parallelizes its Filter/Count actions across
    // worker threads once this is on — this is what makes scanning
    // many large ROOT files fast.
    //
    // Thread count defaults to the machine's hardware concurrency.
    // Override by setting BLPMC_NTHREADS before starting ROOT, e.g.:
    //   BLPMC_NTHREADS=8 root
    {
        UInt_t nThreads = 0;   // 0 = let ROOT pick (hardware_concurrency)
        const char *envN = gSystem->Getenv("BLPMC_NTHREADS");
        if (envN && envN[0] != '\0') nThreads = (UInt_t)atoi(envN);

        ROOT::EnableImplicitMT(nThreads);
        cout << "Enabled ROOT implicit MT with "
             << ROOT::GetThreadPoolSize() << " threads" << endl;
    }

    // ── Include path so all #include "..." inside the files resolve ───
    gInterpreter->AddIncludePath("Include/");
    gROOT->ProcessLine(".I Include/");

    // ── PhysicalConstants: header-only, compiled once ─────────────────
    // Defines mp, mC, E_EXCITATION etc. used by Kinematics and Meyer.
    cout << "Compiling PhysicalConstants..." << endl;
    gROOT->ProcessLine(".L Include/PhysicalConstants.h+");

    // ── DetectorConfig: has both .h and .C — load the .C compiled ─────
    // Defines DETECTOR_THETA_CENTER, DETECTOR_THETA_WINDOW,
    // DETECTOR_PHI_WINDOW etc. used by Kinematics.C and ComputeNormalization.
    // Must be compiled (not just the .h) so the symbols are exported
    // into the shared library that Kinematics_C.so links against.
    cout << "Compiling DetectorConfig..." << endl;
    gROOT->ProcessLine(".L Include/DetectorConfig.h+");
    gROOT->ProcessLine(".L Include/DetectorConfig.C+");

    // ── Kinematics: loaded INTERPRETED (no +) ─────────────────────────
    // Reason: Kinematics.C uses DETECTOR_PHI_WINDOW which is defined in
    // DetectorConfig.h.  When ACLiC compiles Kinematics.C into a .so it
    // cannot resolve that symbol at link time (header constants are not
    // exported as shared library symbols).  Loading interpreted avoids
    // the linker step entirely — the interpreter resolves everything at
    // runtime from the already-loaded DetectorConfig.
    cout << "Loading Kinematics (interpreted)..." << endl;
    gROOT->ProcessLine(".L Include/Kinematics.C");

    // ── MeyerScattering: loaded INTERPRETED (no +) ────────────────────
    // Reason: MeyerScattering.h directly #includes the spline data files
    //   APSpline_160MeV_elastic.C, XSSpline_200MeV_inelastic.C, etc.
    // ACLiC would try to compile those as part of the translation unit
    // and cannot find them even with -I because rootcling does not
    // honour the runtime include path the same way.  Loading interpreted
    // works because the interpreter processes each #include inline.
    cout << "Loading MeyerScattering (interpreted)..." << endl;
    gROOT->ProcessLine(".L Include/MeyerScattering.C");

    // ── Analysis classes: compiled for speed ──────────────────────────
    cout << "Compiling PolHistograms..." << endl;
    gROOT->ProcessLine(".L Include/PolHistograms.h+");

    cout << "Compiling PolAnalysis..." << endl;
    gROOT->ProcessLine(".L Include/PolAnalysis.h+");

    // ── Main script: interpreted so functions are callable interactively
    cout << "Loading analyze_blpmc..." << endl;
    gROOT->ProcessLine(".L analyze_blpmc.C");

    cout << "\n=== Environment ready ===" << endl;
    cout << "  analyze_blpmc_scan(\"ROTY_m0p2_p0p1_0p01\", 199.5)" << endl;
    cout << "  inspect(\"16p2\", 0.0, 0, 199.5)" << endl;
    cout << "  draw_scan_p16()   draw_scan_p16(\"elas\")   draw_scan_p16(\"all\")" << endl;
    cout << "  draw_scan_p12()   draw_scan_both()   draw_scan_all_modes()" << endl;
    cout << "==========================================" << endl;
}