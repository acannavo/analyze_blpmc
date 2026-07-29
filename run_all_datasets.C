// Grabs a canvas by name, saves it, then closes it so the next
// same-named canvas (e.g. cCounts_P16 on the next mode) doesn't collide.
void SaveCanvas(const char *canvName, const TString &outPath) {
    TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvName);
    if (!c) { cout << "[Warning] Canvas not found: " << canvName << endl; return; }
    c->SaveAs(outPath);
    c->Close();
}

void run_all_datasets() {
    std::vector<TString> datasets = {
        "ROTY_m1p5_p1p5_0p1_Pin_200MeV_P80",
        "ROTY_m1p5_p1p5_0p1_200MeV_MT_P80_dTh20_dPh20_MEYER",
        "TRANSX_m1_p1_0p1_200MeV_MT_P80_dTh20_dPh20_MEYER",
        "TRANSY_m1_p1_0p1_200MeV_MT_P80_dTh20_dPh20_MEYER",
        "TRANSY_m1_p1_0p1_Pin_200MeV_P80",
        "TRANSX_m1_p1_0p1_Pin_200MeV_P80",
        "ROTY_m0p2_p0p2_0p01_Pin_200MeV_P80",
        "ROTY_m0p2_p0p2_0p01_200MeV_MT_P80_dTh20_dPh20_MEYER"
    };

    Double_t beamEnergy = 199.5;   // confirm this is the true beam KE, see note below

    for (auto &dataset : datasets) {
        TString outDir = "/home/acannavo/Analysis/" + dataset;
        gSystem->mkdir(outDir, kTRUE);   // kTRUE = create recursively, ok if it already exists

        cout << "\n=== " << dataset << " ===" << endl;
        analyze_blpmc_scan(dataset, beamEnergy);

        draw_scan_p16();       SaveCanvas("cScan_P16_weighted", outDir + "/P16_weighted.png");
        draw_scan_p16("elas"); SaveCanvas("cScan_P16_elas",     outDir + "/P16_elas.png");
        draw_scan_p16("inel"); SaveCanvas("cScan_P16_inel",     outDir + "/P16_inel.png");
        draw_scan_p16("all");  SaveCanvas("cScan_P16_all",      outDir + "/P16_all.png");

        draw_scan_p12();       SaveCanvas("cScan_P12_weighted", outDir + "/P12_weighted.png");
        draw_scan_p12("elas"); SaveCanvas("cScan_P12_elas",     outDir + "/P12_elas.png");
        draw_scan_p12("inel"); SaveCanvas("cScan_P12_inel",     outDir + "/P12_inel.png");
        draw_scan_p12("all");  SaveCanvas("cScan_P12_all",      outDir + "/P12_all.png");

        draw_count_vs_param("elas", true);
        SaveCanvas("cCounts_P12", outDir + "/Counts_P12_elas.png");
        SaveCanvas("cCounts_P16", outDir + "/Counts_P16_elas.png");

        draw_count_vs_param("inel", true);
        SaveCanvas("cCounts_P12", outDir + "/Counts_P12_inel.png");
        SaveCanvas("cCounts_P16", outDir + "/Counts_P16_inel.png");

        draw_count_vs_param("both", true);
        SaveCanvas("cCounts_P12", outDir + "/Counts_P12_both.png");
        SaveCanvas("cCounts_P16", outDir + "/Counts_P16_both.png");
    }
    cout << "\nAll datasets done." << endl;
}