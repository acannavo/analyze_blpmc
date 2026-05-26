// =========================================================
// analyze_blpmc.C — TRANSX Parameter Scan
//
// Automatically scans ../Data/TRANSX/ folder structure:
//   TRANSX/
//     pC_{Elas|Inel443}_200MeV_MT_P80_{SpinUp|SpinDown}_{12p0|16p2}_MEYER/
//       {0p0 .. 20p0}/
//         blpmc_*-*-*_*-*-*.root
//
// For each (angle, TRANSX) point: accumulates elastic + inelastic
// (with normalization factor) into one PolAnalysis object.
// Results: TGraphErrors of Ay vs TRANSX for P12 and P16.
// =========================================================

#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>

// =========================================================
// Configuration
// =========================================================
const Double_t kInelFactor   = 0.1077;   // fraction of inelastic events to use

// P16 (16.2 deg, 4 scintillators) cuts
const Double_t kP16_S0_min = 1.5,  kP16_S0_max = 2.5;
const Double_t kP16_S1_min = 10.0, kP16_S1_max = 14.0;
const Double_t kP16_S2_min = 11.0, kP16_S2_max = 21.0;
const Double_t kP16_S3_min = 15.0, kP16_S3_max = 35.0;

// P12 (12.0 deg, 3 scintillators) cuts
const Double_t kP12_S0_min = 1.5,  kP12_S0_max = 2.5;
const Double_t kP12_S1_min = 3.5,  kP12_S1_max = 5.7;
const Double_t kP12_S2_min = 6.5,  kP12_S2_max = 9.0;

// Beam polarization and theoretical Ay
const Double_t kBeamPol      = 0.80;
const Double_t kP16_TheoAy  = 0.993;
const Double_t kP12_TheoAy  = 0.786;

// =========================================================
// Scan result containers (filled by analyze_blpmc_scan)
// =========================================================
TGraphErrors *gGraph_P16 = nullptr;
TGraphErrors *gGraph_P12 = nullptr;

// =========================================================
// Helper: find blpmc_*.root inside a folder
// =========================================================
TString FindOutputRoot(const TString &runDir) {
    TSystemDirectory dir(runDir.Data(), runDir.Data());
    TList *files = dir.GetListOfFiles();
    if (!files) return "";
    TSystemFile *f;
    TIter next(files);
    while ((f = (TSystemFile*)next())) {
        TString fname = f->GetName();
        if (fname.EndsWith(".root") && fname.BeginsWith("blpmc_"))
            return runDir + "/" + fname;
    }
    return "";
}

// =========================================================
// Helper: open TFile and retrieve TTree "tree"
// =========================================================
TTree* OpenTree(const TString &rootFile, TFile* &fileHandle) {
    fileHandle = TFile::Open(rootFile);
    if (!fileHandle || fileHandle->IsZombie()) {
        cout << "[Error] Cannot open: " << rootFile << endl;
        fileHandle = nullptr;
        return nullptr;
    }
    TTree *t = (TTree*)fileHandle->Get("tree");
    if (!t) {
        cout << "[Error] TTree 'tree' not found in: " << rootFile << endl;
        fileHandle->Close();
        fileHandle = nullptr;
        return nullptr;
    }
    cout << "  Opened: " << rootFile
         << "  (" << t->GetEntries() << " entries)" << endl;
    return t;
}

// =========================================================
// Helper: parse "Xp0" or "Xp2" style strings → Double_t
//   e.g. "0p0" → 0.0,  "16p2" → 16.2,  "1p0" → 1.0
// =========================================================
Double_t ParsePValue(const TString &s) {
    TString tmp = s;
    tmp.ReplaceAll("p", ".");
    return tmp.Atof();
}

// =========================================================
// Helper: close a TFile safely
// =========================================================
void SafeClose(TFile* &f) {
    if (f) { f->Close(); f = nullptr; }
}

// =========================================================
// MAIN SCAN FUNCTION
// =========================================================
void analyze_blpmc_scan() {

    // ----------------------------------------------------------
    // Relative path from Analysis/ to the TRANSX data folder
    // ----------------------------------------------------------
    TString transxDir = "../Data/TRANSX";

    // ----------------------------------------------------------
    // 1. Scan TRANSX/ for subdirectories matching our naming
    //    convention and group by (angle_str, transx_str)
    //
    //    Key format: angle_str → e.g. "12p0" or "16p2"
    //    Sub-key:    transx_str → e.g. "0p0", "1p0", ..., "20p0"
    //    Value:      map of role → full folder path
    //      roles: "elas_up", "elas_down", "inel_up", "inel_down"
    // ----------------------------------------------------------

    // data[angle][transx][role] = folder path
    std::map<TString, std::map<TString, std::map<TString,TString>>> data;

    TSystemDirectory topDir(transxDir.Data(), transxDir.Data());
    TList *plutoFolders = topDir.GetListOfFiles();
    if (!plutoFolders) {
        cout << "[Error] Cannot list: " << transxDir << endl;
        return;
    }

    TSystemFile *plutoEntry;
    TIter nextPluto(plutoFolders);
    while ((plutoEntry = (TSystemFile*)nextPluto())) {
        TString pname = plutoEntry->GetName();
        if (!plutoEntry->IsDirectory()) continue;
        if (pname == "." || pname == "..") continue;

        // Parse folder name: pC_{reac}_200MeV_MT_P80_{spin}_{angle}_MEYER
        // Must contain angle token (12p0 or 16p2) and spin token
        Bool_t isUp    = pname.Contains("SpinUp");
        Bool_t isDown  = pname.Contains("SpinDown");
        Bool_t isElas  = pname.Contains("Elas") && !pname.Contains("Inel");
        Bool_t isInel  = pname.Contains("Inel");
        Bool_t is12    = pname.Contains("_12p0_");
        Bool_t is16    = pname.Contains("_16p2_");

        if ((!isUp && !isDown) || (!isElas && !isInel) || (!is12 && !is16)) {
            cout << "[Skip] Unrecognised folder: " << pname << endl;
            continue;
        }

        TString angle = is12 ? "12p0" : "16p2";
        TString role;
        if      (isElas && isUp)   role = "elas_up";
        else if (isElas && isDown) role = "elas_down";
        else if (isInel && isUp)   role = "inel_up";
        else if (isInel && isDown) role = "inel_down";

        // Now scan subfolders for TRANSX values (e.g. 0p0, 1p0, ...)
        TString plutoPath = transxDir + "/" + pname;
        TSystemDirectory plutoDir(plutoPath.Data(), plutoPath.Data());
        TList *txFolders = plutoDir.GetListOfFiles();
        if (!txFolders) continue;

        TSystemFile *txEntry;
        TIter nextTx(txFolders);
        while ((txEntry = (TSystemFile*)nextTx())) {
            TString txname = txEntry->GetName();
            if (!txEntry->IsDirectory()) continue;
            if (txname == "." || txname == "..") continue;
            // Expect names like "0p0", "1p0", ..., "20p0"
            if (!txname.Contains("p")) continue;

            TString fullPath = plutoPath + "/" + txname;
            data[angle][txname][role] = fullPath;
        }
    }

    if (data.empty()) {
        cout << "[Error] No valid folders found under " << transxDir << endl;
        return;
    }

    // ----------------------------------------------------------
    // 2. Collect sorted TRANSX values and angles
    // ----------------------------------------------------------
    std::vector<TString> angles;
    for (auto &a : data) angles.push_back(a.first);
    sort(angles.begin(), angles.end());

    // Get all transx keys from the first angle
    std::vector<TString> txKeys;
    for (auto &t : data[angles[0]]) txKeys.push_back(t.first);
    // Sort by numeric value
    sort(txKeys.begin(), txKeys.end(), [](const TString &a, const TString &b){
        return ParsePValue(a) < ParsePValue(b);
    });

    cout << "\n=== TRANSX Scan ===" << endl;
    cout << "Angles found:  ";
    for (auto &a : angles) cout << ParsePValue(a) << "  ";
    cout << endl;
    cout << "TRANSX values: ";
    for (auto &t : txKeys) cout << ParsePValue(t) << "  ";
    cout << endl;
    cout << "Inelastic normalization factor: " << kInelFactor << endl;
    cout << "===================" << endl;

    // ----------------------------------------------------------
    // 3. Build one PolAnalysis per angle (reused across TRANSX)
    //    and two TGraphErrors for the output
    // ----------------------------------------------------------
    Int_t nPoints = (Int_t)txKeys.size();

    gGraph_P16 = new TGraphErrors(nPoints);
    gGraph_P16->SetName("gAy_P16");
    gGraph_P16->SetTitle("Analyzing Power vs TRANSX (P16, 16.2#circ);TRANSX [mm];A_{y}");
    gGraph_P16->SetMarkerStyle(21);
    gGraph_P16->SetMarkerColor(kBlue+1);
    gGraph_P16->SetLineColor(kBlue+1);

    gGraph_P12 = new TGraphErrors(nPoints);
    gGraph_P12->SetName("gAy_P12");
    gGraph_P12->SetTitle("Analyzing Power vs TRANSX (P12, 12.0#circ);TRANSX [mm];A_{y}");
    gGraph_P12->SetMarkerStyle(20);
    gGraph_P12->SetMarkerColor(kRed+1);
    gGraph_P12->SetLineColor(kRed+1);

    // ----------------------------------------------------------
    // 4. Loop over TRANSX values
    // ----------------------------------------------------------
    for (Int_t ipt = 0; ipt < nPoints; ipt++) {
        TString txKey = txKeys[ipt];
        Double_t txVal = ParsePValue(txKey);

        cout << "\n=============================" << endl;
        cout << " TRANSX = " << txVal << endl;
        cout << "=============================" << endl;

        // We'll compute Ay for each angle and store in the graphs
        // Process P16 (16p2) and P12 (12p0)
        struct AngleResult { Double_t ay, ayErr; };
        std::map<TString, AngleResult> results;

        for (auto &angle : angles) {
            cout << "\n  --- Angle: " << ParsePValue(angle) << " deg ---" << endl;

            auto &roleMap = data[angle][txKey];

            // Check all 4 roles exist
            Bool_t ok = (roleMap.count("elas_up")   &&
                         roleMap.count("elas_down") &&
                         roleMap.count("inel_up")   &&
                         roleMap.count("inel_down"));
            if (!ok) {
                cout << "[Warning] Missing folders for angle=" << angle
                     << " TRANSX=" << txKey << " — skipping." << endl;
                results[angle] = {0.0, 0.0};
                continue;
            }

            // Find ROOT files
            TString fEU = FindOutputRoot(roleMap["elas_up"]);
            TString fED = FindOutputRoot(roleMap["elas_down"]);
            TString fIU = FindOutputRoot(roleMap["inel_up"]);
            TString fID = FindOutputRoot(roleMap["inel_down"]);

            if (fEU.IsNull()||fED.IsNull()||fIU.IsNull()||fID.IsNull()) {
                cout << "[Warning] Missing ROOT file for angle=" << angle
                     << " TRANSX=" << txKey << " — skipping." << endl;
                results[angle] = {0.0, 0.0};
                continue;
            }

            // Open trees
            TFile *hEU=nullptr, *hED=nullptr, *hIU=nullptr, *hID=nullptr;
            TTree *tEU = OpenTree(fEU, hEU);
            TTree *tED = OpenTree(fED, hED);
            TTree *tIU = OpenTree(fIU, hIU);
            TTree *tID = OpenTree(fID, hID);

            if (!tEU||!tED||!tIU||!tID) {
                cout << "[Warning] Failed to open trees for angle=" << angle
                     << " TRANSX=" << txKey << " — skipping." << endl;
                SafeClose(hEU); SafeClose(hED);
                SafeClose(hIU); SafeClose(hID);
                results[angle] = {0.0, 0.0};
                continue;
            }

            // Inelastic normalization
            Long64_t nIU_total = tIU->GetEntries();
            Long64_t nID_total = tID->GetEntries();
            Long64_t maxIU = (Long64_t)TMath::Nint(kInelFactor * nIU_total);
            Long64_t maxID = (Long64_t)TMath::Nint(kInelFactor * nID_total);

            cout << "  Inel SpinUp:   " << nIU_total << " total -> using " << maxIU << endl;
            cout << "  Inel SpinDown: " << nID_total << " total -> using " << maxID << endl;

            // Create PolAnalysis for this angle
            Bool_t isP16 = angle.Contains("16");
            TString detName = isP16 ? "P16" : "P12";
            Int_t   nScint  = isP16 ? 4 : 3;

            PolAnalysis *ana = new PolAnalysis(detName, nScint);
            if (isP16) {
                ana->SetCuts(0, kP16_S0_min, kP16_S0_max);
                ana->SetCuts(1, kP16_S1_min, kP16_S1_max);
                ana->SetCuts(2, kP16_S2_min, kP16_S2_max);
                ana->SetCuts(3, kP16_S3_min, kP16_S3_max);
                ana->SetTheoreticalAy(kP16_TheoAy);
            } else {
                ana->SetCuts(0, kP12_S0_min, kP12_S0_max);
                ana->SetCuts(1, kP12_S1_min, kP12_S1_max);
                ana->SetCuts(2, kP12_S2_min, kP12_S2_max);
                ana->SetTheoreticalAy(kP12_TheoAy);
            }
            ana->SetBeamPolarization(kBeamPol);
            ana->SetBeamPolarizationError(0.0);

            // Count events: elastic (all) + inelastic (normalized)
            cout << "  Counting SpinUp elastic..." << endl;
            ana->CountEvents(tEU, true);
            cout << "  Counting SpinUp inelastic (normalized)..." << endl;
            ana->CountEvents(tIU, true,  maxIU);
            cout << "  Counting SpinDown elastic..." << endl;
            ana->CountEvents(tED, false);
            cout << "  Counting SpinDown inelastic (normalized)..." << endl;
            ana->CountEvents(tID, false, maxID);

            ana->PrintResults();

            results[angle] = { ana->GetAnalyzingPower_CrossRatio(),
                               ana->GetAnalyzingPowerError_CrossRatio() };

            delete ana;
            SafeClose(hEU); SafeClose(hED);
            SafeClose(hIU); SafeClose(hID);
        }

        // Store in graphs
        for (auto &angle : angles) {
            Bool_t isP16 = angle.Contains("16");
            TGraphErrors *g = isP16 ? gGraph_P16 : gGraph_P12;
            g->SetPoint(ipt, txVal, results[angle].ay);
            g->SetPointError(ipt, 0.0, results[angle].ayErr);
        }
    }

    cout << "\n=== Scan complete ===" << endl;
    cout << "Call draw_scan_p16() or draw_scan_p12() to view results." << endl;
}

// =========================================================
// DRAW FUNCTIONS
// =========================================================
void draw_scan_p16() {
    if (!gGraph_P16) {
        cout << "Please run analyze_blpmc_scan() first!" << endl;
        return;
    }
    TCanvas *c = new TCanvas("cScan_P16", "Ay vs TRANSX — P16 (16.2 deg)", 900, 600);
    c->SetGrid();
    gGraph_P16->Draw("AP");
    gGraph_P16->GetXaxis()->SetTitle("TRANSX [mm]");
    gGraph_P16->GetYaxis()->SetTitle("A_{y}");

    // Draw horizontal line at theoretical Ay
    TLine *lTheo = new TLine(gGraph_P16->GetXaxis()->GetXmin(), kP16_TheoAy,
                             gGraph_P16->GetXaxis()->GetXmax(), kP16_TheoAy);
    lTheo->SetLineStyle(2);
    lTheo->SetLineColor(kGray+2);
    lTheo->SetLineWidth(2);

    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(gGraph_P16, "Measured A_{y}", "lep");
    leg->AddEntry(lTheo, Form("Theory A_{y} = %.3f", kP16_TheoAy), "l");
    leg->Draw();

    c->Update();
}

void draw_scan_p12() {
    if (!gGraph_P12) {
        cout << "Please run analyze_blpmc_scan() first!" << endl;
        return;
    }
    TCanvas *c = new TCanvas("cScan_P12", "Ay vs TRANSX — P12 (12.0 deg)", 900, 600);
    c->SetGrid();
    gGraph_P12->Draw("AP");
    gGraph_P12->GetXaxis()->SetTitle("TRANSX [mm]");
    gGraph_P12->GetYaxis()->SetTitle("A_{y}");

    TLine *lTheo = new TLine(gGraph_P12->GetXaxis()->GetXmin(), kP12_TheoAy,
                             gGraph_P12->GetXaxis()->GetXmax(), kP12_TheoAy);
    lTheo->SetLineStyle(2);
    lTheo->SetLineColor(kGray+2);
    lTheo->SetLineWidth(2);
    lTheo->Draw();

    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(gGraph_P12, "Measured A_{y}", "lep");
    leg->AddEntry(lTheo, Form("Theory A_{y} = %.3f", kP12_TheoAy), "l");
    leg->Draw();

    c->Update();
}

void draw_scan_both() {
    draw_scan_p16();
    draw_scan_p12();
}

// =========================================================
// INSPECTION MODE
// Opens one (angle, TRANSX) point, fills PolHistograms,
// and exposes all 1D/2D draw functions for that point.
//
// Usage:
//   inspect("12p0", 0.0)   → P12, TRANSX=0
//   inspect("16p2", 10.0)  → P16, TRANSX=10
//   inspect("16p2", 0.0)   → P16, TRANSX=0
//
// After calling inspect(), use:
//   insp_draw_1d_raw()
//   insp_draw_1d_comparison()
//   insp_draw_2d()
//   insp_draw_all()
// =========================================================

// Global inspection objects (reused across calls)
PolHistograms *gInsp_hist    = nullptr;
PolAnalysis   *gInsp_ana     = nullptr;
TString        gInsp_label   = "";

// Helper: fill histograms from one tree (elastic or inelastic)
void InspFillTree(TTree *tree, Long64_t maxEntries,
                  PolHistograms *hist, PolAnalysis *ana,
                  Int_t nScint) {

    Double_t eDepP16LS0, eDepP16LS1, eDepP16LS2, eDepP16LS3;
    Double_t eDepP16RS0, eDepP16RS1, eDepP16RS2, eDepP16RS3;
    Double_t eDepP12LS0, eDepP12LS1, eDepP12LS2;
    Double_t eDepP12RS0, eDepP12RS1, eDepP12RS2;

    tree->SetBranchAddress("eDepP16LS0", &eDepP16LS0);
    tree->SetBranchAddress("eDepP16LS1", &eDepP16LS1);
    tree->SetBranchAddress("eDepP16LS2", &eDepP16LS2);
    tree->SetBranchAddress("eDepP16LS3", &eDepP16LS3);
    tree->SetBranchAddress("eDepP16RS0", &eDepP16RS0);
    tree->SetBranchAddress("eDepP16RS1", &eDepP16RS1);
    tree->SetBranchAddress("eDepP16RS2", &eDepP16RS2);
    tree->SetBranchAddress("eDepP16RS3", &eDepP16RS3);
    tree->SetBranchAddress("eDepP12LS0", &eDepP12LS0);
    tree->SetBranchAddress("eDepP12LS1", &eDepP12LS1);
    tree->SetBranchAddress("eDepP12LS2", &eDepP12LS2);
    tree->SetBranchAddress("eDepP12RS0", &eDepP12RS0);
    tree->SetBranchAddress("eDepP12RS1", &eDepP12RS1);
    tree->SetBranchAddress("eDepP12RS2", &eDepP12RS2);

    Long64_t nentries = tree->GetEntries();
    if (maxEntries > 0 && maxEntries < nentries) nentries = maxEntries;

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        if (nScint == 4) {
            Double_t e[8] = {eDepP16LS0, eDepP16LS1, eDepP16LS2, eDepP16LS3,
                             eDepP16RS0, eDepP16RS1, eDepP16RS2, eDepP16RS3};
            for (Int_t j = 0; j < 8; j++) hist->Fill(j, e);
            hist->Fill2D(e);
            for (Int_t j = 0; j < 8; j++) hist->FillWithCuts(j, e, ana);
        } else {
            Double_t e[6] = {eDepP12LS0, eDepP12LS1, eDepP12LS2,
                             eDepP12RS0, eDepP12RS1, eDepP12RS2};
            for (Int_t j = 0; j < 6; j++) hist->Fill(j, e);
            hist->Fill2D(e);
            for (Int_t j = 0; j < 6; j++) hist->FillWithCuts(j, e, ana);
        }
    }
}

void inspect(TString angleKey, Double_t txVal, Int_t spin = 0) {
    // spin:  0 = all combined (default), +1 = SpinUp only, -1 = SpinDown only

    // ----------------------------------------------------------
    // Build folder tokens from arguments
    // ----------------------------------------------------------
    // angleKey: "12p0" or "16p2"
    Bool_t isP16   = angleKey.Contains("16");
    TString detName = isP16 ? "P16" : "P12";
    Int_t   nScint  = isP16 ? 4 : 3;
    TString angleFull = isP16 ? "16p2" : "12p0";

    // Convert txVal to folder name: 0.0→"0p0", 10.0→"10p0"
    TString txKey = Form("%gp0", txVal);
    txKey.ReplaceAll(".", "p");  // in case of decimal input

    TString dataPath = "../Data/TRANSX";
    TString label    = Form("%s_TRANSX%g", detName.Data(), txVal);
    TString spinTag = (spin == 1) ? "_SpinUp" : (spin == -1) ? "_SpinDown" : "_AllSpin";
    gInsp_label = label + spinTag;

    cout << "\n=== Inspection: " << detName << "  TRANSX=" << txVal << " ===" << endl;

    // ----------------------------------------------------------
    // Build the 4 folder paths
    // ----------------------------------------------------------
    auto MakePath = [&](TString reac, TString spin) -> TString {
        return Form("%s/pC_%s_200MeV_MT_P80_%s_%s_MEYER/%s",
                    dataPath.Data(), reac.Data(), spin.Data(),
                    angleFull.Data(), txKey.Data());
    };

    TString dirEU = MakePath("Elas",    "SpinUp");
    TString dirED = MakePath("Elas",    "SpinDown");
    TString dirIU = MakePath("Inel443", "SpinUp");
    TString dirID = MakePath("Inel443", "SpinDown");

    TString fEU = FindOutputRoot(dirEU);
    TString fED = FindOutputRoot(dirED);
    TString fIU = FindOutputRoot(dirIU);
    TString fID = FindOutputRoot(dirID);

    if (fEU.IsNull()||fED.IsNull()||fIU.IsNull()||fID.IsNull()) {
        cout << "[Error] Could not find ROOT files. Check angle/TRANSX values." << endl;
        cout << "  Looked in:" << endl;
        cout << "    " << dirEU << endl;
        cout << "    " << dirED << endl;
        cout << "    " << dirIU << endl;
        cout << "    " << dirID << endl;
        return;
    }

    // ----------------------------------------------------------
    // Open trees
    // ----------------------------------------------------------
    TFile *hEU=nullptr, *hED=nullptr, *hIU=nullptr, *hID=nullptr;
    TTree *tEU = OpenTree(fEU, hEU);
    TTree *tED = OpenTree(fED, hED);
    TTree *tIU = OpenTree(fIU, hIU);
    TTree *tID = OpenTree(fID, hID);

    if (!tEU||!tED||!tIU||!tID) {
        cout << "[Error] Failed to open one or more trees." << endl;
        SafeClose(hEU); SafeClose(hED);
        SafeClose(hIU); SafeClose(hID);
        return;
    }

    Long64_t maxIU = (Long64_t)TMath::Nint(kInelFactor * tIU->GetEntries());
    Long64_t maxID = (Long64_t)TMath::Nint(kInelFactor * tID->GetEntries());
    cout << "  Inel SpinUp:   " << tIU->GetEntries() << " total -> using " << maxIU << endl;
    cout << "  Inel SpinDown: " << tID->GetEntries() << " total -> using " << maxID << endl;

    // ----------------------------------------------------------
    // (Re)create global inspection objects
    // ----------------------------------------------------------
    if (gInsp_hist) { delete gInsp_hist; gInsp_hist = nullptr; }
    if (gInsp_ana)  { delete gInsp_ana;  gInsp_ana  = nullptr; }

    Double_t emin = 0, emax = 50;
    Int_t    nbins = 500;

    TH1::AddDirectory(kFALSE);  // histograms won't be owned by any TFile
    gInsp_hist = new PolHistograms(detName, detName, nScint, nbins, emin, emax);
    gInsp_hist->Create2DHistograms();
    TH1::AddDirectory(kTRUE);   // restore default for everything else

    gInsp_ana = new PolAnalysis(detName, nScint);
    if (isP16) {
        gInsp_ana->SetCuts(0, kP16_S0_min, kP16_S0_max);
        gInsp_ana->SetCuts(1, kP16_S1_min, kP16_S1_max);
        gInsp_ana->SetCuts(2, kP16_S2_min, kP16_S2_max);
        gInsp_ana->SetCuts(3, kP16_S3_min, kP16_S3_max);
        gInsp_ana->SetTheoreticalAy(kP16_TheoAy);
    } else {
        gInsp_ana->SetCuts(0, kP12_S0_min, kP12_S0_max);
        gInsp_ana->SetCuts(1, kP12_S1_min, kP12_S1_max);
        gInsp_ana->SetCuts(2, kP12_S2_min, kP12_S2_max);
        gInsp_ana->SetTheoreticalAy(kP12_TheoAy);
    }
    gInsp_ana->SetBeamPolarization(kBeamPol);
    gInsp_ana->SetBeamPolarizationError(0.0);

    // ----------------------------------------------------------
    // Fill histograms — filtered by spin argument
    //   spin =  0: all four trees (default)
    //   spin = +1: SpinUp only
    //   spin = -1: SpinDown only
    // ----------------------------------------------------------
    if (spin >= 0) {
        cout << "  Filling SpinUp elastic..." << endl;
        InspFillTree(tEU, -1,    gInsp_hist, gInsp_ana, nScint);
        cout << "  Filling SpinUp inelastic (normalized)..." << endl;
        InspFillTree(tIU, maxIU, gInsp_hist, gInsp_ana, nScint);
    }
    if (spin <= 0) {
        cout << "  Filling SpinDown elastic..." << endl;
        InspFillTree(tED, -1,    gInsp_hist, gInsp_ana, nScint);
        cout << "  Filling SpinDown inelastic (normalized)..." << endl;
        InspFillTree(tID, maxID, gInsp_hist, gInsp_ana, nScint);
    }

    // ----------------------------------------------------------
    // Count events for analysis (always use both spins for Ay)
    // ----------------------------------------------------------
    gInsp_ana->ResetCounts();
    gInsp_ana->CountEvents(tEU, true);
    gInsp_ana->CountEvents(tIU, true,  maxIU);
    gInsp_ana->CountEvents(tED, false);
    gInsp_ana->CountEvents(tID, false, maxID);
    gInsp_ana->PrintResults();

    SafeClose(hEU); SafeClose(hED);
    SafeClose(hIU); SafeClose(hID);

    cout << "\n=== Inspection ready! ===" << endl;
    cout << "  insp_draw_1d_raw()        - Raw 1D spectra" << endl;
    cout << "  insp_draw_1d_comparison() - Raw vs cut overlay with cut lines" << endl;
    cout << "  insp_draw_2d()            - 2D correlation plots" << endl;
    cout << "  insp_draw_all()           - All of the above" << endl;
}

// =========================================================
// Inspection draw functions
// =========================================================
void insp_draw_1d_raw() {
    if (!gInsp_hist) { cout << "Run inspect() first!" << endl; return; }
    TString title = Form("Raw Spectra — %s", gInsp_label.Data());
    TCanvas *c = new TCanvas(Form("c_raw_%s", gInsp_label.Data()), title, 1400, 700);
    gInsp_hist->DrawRaw(c);
}

void insp_draw_1d_comparison() {
    if (!gInsp_hist || !gInsp_ana) { cout << "Run inspect() first!" << endl; return; }
    TString title = Form("Raw vs Cuts — %s", gInsp_label.Data());
    TCanvas *c = new TCanvas(Form("c_comp_%s", gInsp_label.Data()), title, 1400, 700);
    gInsp_hist->DrawComparison(c, gInsp_ana);
}

void insp_draw_2d() {
    if (!gInsp_hist) { cout << "Run inspect() first!" << endl; return; }
    TString title = Form("2D Correlations — %s", gInsp_label.Data());
    TCanvas *c = new TCanvas(Form("c_2d_%s", gInsp_label.Data()), title, 1400, 700);
    gInsp_hist->Draw2D(c);
}

void insp_draw_all() {
    insp_draw_1d_raw();
    insp_draw_1d_comparison();
    insp_draw_2d();
}