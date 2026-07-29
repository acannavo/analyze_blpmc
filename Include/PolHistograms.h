#ifndef POLHISTOGRAMS_H
#define POLHISTOGRAMS_H

#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"
#include <vector>
#include <iostream>

// Forward declaration
class PolAnalysis;

class PolHistograms {
private:
    TString name;
    TString title;
    Int_t nbins;
    Double_t emin, emax;

    std::vector<TH1F*> hRaw;
    std::vector<TH1F*> hCut;
    std::vector<TH1F*> hCutUp;    // after cuts, spin-up
    std::vector<TH1F*> hCutDown;  // after cuts, spin-down
    std::vector<TH2D*> h2D;

    Int_t nScintillators;

public:
    PolHistograms(TString _name, TString _title, Int_t _nScint,
                  Int_t _nbins, Double_t _emin, Double_t _emax) {
        name  = _name;
        title = _title;
        nScintillators = _nScint;
        nbins = _nbins;
        emin  = _emin;
        emax  = _emax;

        for (Int_t i = 0; i < nScintillators * 2; i++) {
            TString side  = (i < nScintillators) ? "L" : "R";
            Int_t   scint = i % nScintillators;

            // Raw
            TH1F *h = new TH1F(Form("h%s%sS%d", name.Data(), side.Data(), scint),
                               Form("%s %s S%d;Energy [MeV];Counts",
                                    title.Data(), side.Data(), scint),
                               nbins, emin, emax);
            hRaw.push_back(h);

            // Cut (both spins combined — kept for DrawComparison)
            TH1F *hc = new TH1F(Form("h%s%sS%d_cut", name.Data(), side.Data(), scint),
                                Form("%s %s S%d (with cuts);Energy [MeV];Counts",
                                     title.Data(), side.Data(), scint),
                                nbins, emin, emax);
            hc->SetLineColor(kRed);
            hCut.push_back(hc);

            // Cut spin-up (blue)
            TH1F *hu = new TH1F(Form("h%s%sS%d_up", name.Data(), side.Data(), scint),
                                Form("%s %s S%d Spin#uparrow;Energy [MeV];Counts",
                                     title.Data(), side.Data(), scint),
                                nbins, emin, emax);
            hu->SetLineColor(kBlue+1);
            hu->SetLineWidth(2);
            hCutUp.push_back(hu);

            // Cut spin-down (red)
            TH1F *hd = new TH1F(Form("h%s%sS%d_down", name.Data(), side.Data(), scint),
                                Form("%s %s S%d Spin#downarrow;Energy [MeV];Counts",
                                     title.Data(), side.Data(), scint),
                                nbins, emin, emax);
            hd->SetLineColor(kRed+1);
            hd->SetLineWidth(2);
            hCutDown.push_back(hd);
        }
    }

    void Fill(Int_t idx, Double_t *energies, Bool_t applyNonZero = true) {
        if (applyNonZero && energies[idx] <= 0) return;
        hRaw[idx]->Fill(energies[idx]);
    }

    // Original FillWithCuts — fills the combined hCut (both spins)
    void FillWithCuts(Int_t idx, Double_t *energies, PolAnalysis* analysis);

    // New: fills hCutUp or hCutDown depending on isSpinUp
    void FillWithCutsSpin(Int_t idx, Double_t *energies,
                          PolAnalysis* analysis, Bool_t isSpinUp);

    void Create2DHistograms() {
        if (nScintillators == 4) {
            h2D.push_back(new TH2D(Form("h%sLS1S0", name.Data()),
                                   Form("%s Left S1 vs S0;S0 Energy [MeV];S1 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sLS2S1", name.Data()),
                                   Form("%s Left S2 vs S1;S1 Energy [MeV];S2 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sLS3S2", name.Data()),
                                   Form("%s Left S3 vs S2;S2 Energy [MeV];S3 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sLS3S1", name.Data()),
                                   Form("%s Left S3 vs S1;S1 Energy [MeV];S3 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS1S0", name.Data()),
                                   Form("%s Right S1 vs S0;S0 Energy [MeV];S1 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS2S1", name.Data()),
                                   Form("%s Right S2 vs S1;S1 Energy [MeV];S2 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS3S2", name.Data()),
                                   Form("%s Right S3 vs S2;S2 Energy [MeV];S3 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS3S1", name.Data()),
                                   Form("%s Right S3 vs S1;S1 Energy [MeV];S3 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
        } else if (nScintillators == 3) {
            h2D.push_back(new TH2D(Form("h%sLS1S0", name.Data()),
                                   Form("%s Left S1 vs S0;S0 Energy [MeV];S1 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sLS2S1", name.Data()),
                                   Form("%s Left S2 vs S1;S1 Energy [MeV];S2 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sLS2S0", name.Data()),
                                   Form("%s Left S2 vs S0;S0 Energy [MeV];S2 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sDummy1", name.Data()),
                                   "Dummy",
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS1S0", name.Data()),
                                   Form("%s Right S1 vs S0;S0 Energy [MeV];S1 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS2S1", name.Data()),
                                   Form("%s Right S2 vs S1;S1 Energy [MeV];S2 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS2S0", name.Data()),
                                   Form("%s Right S2 vs S0;S0 Energy [MeV];S2 Energy [MeV]",
                                        title.Data()),
                                   nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sDummy2", name.Data()),
                                   "Dummy",
                                   nbins, emin, emax, nbins, emin, emax));
        }
    }

    void Fill2D(Double_t *energies) {
        Double_t threshold = 0.1;

        if (nScintillators == 4) {
            if (energies[0] > threshold && energies[1] > threshold) h2D[0]->Fill(energies[0], energies[1]);
            if (energies[1] > threshold && energies[2] > threshold) h2D[1]->Fill(energies[1], energies[2]);
            if (energies[2] > threshold && energies[3] > threshold) h2D[2]->Fill(energies[2], energies[3]);
            if (energies[1] > threshold && energies[3] > threshold) h2D[3]->Fill(energies[1], energies[3]);
            if (energies[4] > threshold && energies[5] > threshold) h2D[4]->Fill(energies[4], energies[5]);
            if (energies[5] > threshold && energies[6] > threshold) h2D[5]->Fill(energies[5], energies[6]);
            if (energies[6] > threshold && energies[7] > threshold) h2D[6]->Fill(energies[6], energies[7]);
            if (energies[5] > threshold && energies[7] > threshold) h2D[7]->Fill(energies[5], energies[7]);
        } else if (nScintillators == 3) {
            if (energies[0] > threshold && energies[1] > threshold) h2D[0]->Fill(energies[0], energies[1]);
            if (energies[1] > threshold && energies[2] > threshold) h2D[1]->Fill(energies[1], energies[2]);
            if (energies[0] > threshold && energies[2] > threshold) h2D[2]->Fill(energies[0], energies[2]);
            if (energies[3] > threshold && energies[4] > threshold) h2D[4]->Fill(energies[3], energies[4]);
            if (energies[4] > threshold && energies[5] > threshold) h2D[5]->Fill(energies[4], energies[5]);
            if (energies[3] > threshold && energies[5] > threshold) h2D[6]->Fill(energies[3], energies[5]);
        }
    }

    void DrawRaw(TCanvas *c) {
        Int_t nx = (nScintillators == 4) ? 4 : 3;
        c->Divide(nx, 2);
        for (Int_t i = 0; i < nScintillators * 2; i++) {
            c->cd(i + 1);
            hRaw[i]->Draw();
        }
        c->Update();
    }

    void DrawCuts(TCanvas *c, PolAnalysis* analysis);

    void DrawComparison(TCanvas *c, PolAnalysis* analysis);

    // New: spin-up (blue) vs spin-down (red) after cuts
    void DrawSpinComparison(TCanvas *c, PolAnalysis* analysis);

    // New: overlay THIS object's raw spectrum (treated as elastic) against
    // another PolHistograms' raw spectrum (treated as inelastic), with cut
    // window lines from `analysis`. Used to visually separate which
    // population (elastic vs inelastic) contributes inside the cut window —
    // something a single combined histogram cannot show.
    // R > 0: the inelastic histogram is rescaled (display only, per
    // scintillator) so its integral relative to the elastic one equals
    // the physical cross-section ratio R. R <= 0: drawn unscaled.
    void DrawElasInelComparison(TCanvas *c, PolHistograms* inelHist, PolAnalysis* analysis,
                                 Double_t R = -1.0);

    // Same as DrawElasInelComparison but shows POST-CUT spectra instead
    // of raw — lets you see directly which events from each population
    // actually survive the full coincidence cut.
    void DrawElasInelCutComparison(TCanvas *c, PolHistograms* inelHist, PolAnalysis* analysis,
                                    Double_t R = -1.0);

    // Public accessor needed by DrawElasInelComparison to read another
    // instance's raw histograms (hRaw is private).
    TH1F* GetRawHist(Int_t idx) const {
        return (idx >= 0 && idx < (Int_t)hRaw.size()) ? hRaw[idx] : nullptr;
    }

    // Public accessor for post-cut histograms (hCut is private).
    // Used by DrawElasInelCutComparison.
    TH1F* GetCutHist(Int_t idx) const {
        return (idx >= 0 && idx < (Int_t)hCut.size()) ? hCut[idx] : nullptr;
    }

    void Draw2D(TCanvas *c) {
        c->Divide(4, 2);
        for (Int_t i = 0; i < 8; i++) {
            c->cd(i + 1);
            gPad->SetLogz();
            h2D[i]->Draw("colz");
        }
        c->Update();
    }
};

// Include PolAnalysis after PolHistograms is defined
#include "PolAnalysis.h"

// =========================================================
// FillWithCuts — combined (both spins), used by DrawComparison
// =========================================================
void PolHistograms::FillWithCuts(Int_t idx, Double_t *energies,
                                  PolAnalysis* analysis) {
    if (!analysis) return;
    if (energies[idx] <= 0) return;

    Int_t  scint  = idx % nScintillators;
    Bool_t isLeft = (idx < nScintillators);
    Int_t  offset = isLeft ? 0 : nScintillators;

    for (Int_t i = 0; i < nScintillators; i++) {
        if (i == scint) continue;
        Double_t e = energies[offset + i];
        if (e < analysis->GetCutMin(i) || e > analysis->GetCutMax(i)) return;
    }
    hCut[idx]->Fill(energies[idx]);
}

// =========================================================
// FillWithCutsSpin — fills hCutUp or hCutDown
// =========================================================
void PolHistograms::FillWithCutsSpin(Int_t idx, Double_t *energies,
                                      PolAnalysis* analysis, Bool_t isSpinUp) {
    if (!analysis) return;
    if (energies[idx] <= 0) return;

    Int_t  scint  = idx % nScintillators;
    Bool_t isLeft = (idx < nScintillators);
    Int_t  offset = isLeft ? 0 : nScintillators;

    for (Int_t i = 0; i < nScintillators; i++) {
        if (i == scint) continue;
        Double_t e = energies[offset + i];
        if (e < analysis->GetCutMin(i) || e > analysis->GetCutMax(i)) return;
    }

    if (isSpinUp) hCutUp[idx]->Fill(energies[idx]);
    else          hCutDown[idx]->Fill(energies[idx]);
}

// =========================================================
// DrawCuts
// =========================================================
void PolHistograms::DrawCuts(TCanvas *c, PolAnalysis* analysis) {
    if (!analysis) {
        std::cout << "Error: PolAnalysis pointer is null!" << std::endl;
        return;
    }

    Int_t nx = (nScintillators == 4) ? 4 : 3;
    c->Divide(nx, 2);
    TLatex latex;
    latex.SetTextSize(0.04);

    for (Int_t i = 0; i < nScintillators * 2; i++) {
        c->cd(i + 1);
        hCut[i]->Draw();

        Int_t    scint = i % nScintillators;
        Double_t ypos  = 0.85;
        for (Int_t j = 0; j < nScintillators; j++) {
            if (j == scint) continue;
            latex.DrawLatexNDC(0.15, ypos, Form("S%d: %.1f-%.1f", j,
                               analysis->GetCutMin(j), analysis->GetCutMax(j)));
            ypos -= 0.05;
        }
    }
    c->Update();
}

// =========================================================
// DrawComparison — raw (black) vs cut (red) with cut lines
// =========================================================
void PolHistograms::DrawComparison(TCanvas *c, PolAnalysis* analysis) {
    if (!analysis) {
        std::cout << "Error: PolAnalysis pointer is null!" << std::endl;
        return;
    }

    Int_t nx = (nScintillators == 4) ? 4 : 3;
    c->Divide(nx, 2);

    for (Int_t i = 0; i < nScintillators * 2; i++) {
        c->cd(i + 1);
        hRaw[i]->Draw();
        hCut[i]->Draw("same");

        Int_t scint = i % nScintillators;
        Double_t ymax = hRaw[i]->GetMaximum() * 1.1;

        TLine *lineMin = new TLine(analysis->GetCutMin(scint), 0,
                                   analysis->GetCutMin(scint), ymax);
        lineMin->SetLineStyle(2);
        lineMin->SetLineColor(kBlue);
        lineMin->SetLineWidth(2);
        lineMin->Draw();

        TLine *lineMax = new TLine(analysis->GetCutMax(scint), 0,
                                   analysis->GetCutMax(scint), ymax);
        lineMax->SetLineStyle(2);
        lineMax->SetLineColor(kBlue);
        lineMax->SetLineWidth(2);
        lineMax->Draw();

        if (i == 0) {
            TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
            leg->AddEntry(hRaw[i], "Raw",        "l");
            leg->AddEntry(hCut[i], "With cuts",  "l");
            leg->Draw();
        }
    }
    c->Update();
}

// =========================================================
// DrawSpinComparison — spin-up (blue) vs spin-down (red)
//                      after energy cuts
// =========================================================
void PolHistograms::DrawSpinComparison(TCanvas *c, PolAnalysis* analysis) {
    if (!analysis) {
        std::cout << "Error: PolAnalysis pointer is null!" << std::endl;
        return;
    }

    Int_t nx = (nScintillators == 4) ? 4 : 3;
    c->Divide(nx, 2);

    for (Int_t i = 0; i < nScintillators * 2; i++) {
        c->cd(i + 1);

        // Scale to the taller of the two so both fit in the frame
        Double_t ymax = TMath::Max(hCutUp[i]->GetMaximum(),
                                   hCutDown[i]->GetMaximum()) * 1.15;
        if (ymax <= 0) ymax = 1.0;

        hCutUp[i]->SetMaximum(ymax);
        hCutUp[i]->Draw("hist");
        hCutDown[i]->Draw("hist same");

        // Draw cut boundaries for this scintillator
        Int_t scint = i % nScintillators;

        TLine *lineMin = new TLine(analysis->GetCutMin(scint), 0,
                                   analysis->GetCutMin(scint), ymax);
        lineMin->SetLineStyle(2);
        lineMin->SetLineColor(kGray+1);
        lineMin->SetLineWidth(1);
        lineMin->Draw();

        TLine *lineMax = new TLine(analysis->GetCutMax(scint), 0,
                                   analysis->GetCutMax(scint), ymax);
        lineMax->SetLineStyle(2);
        lineMax->SetLineColor(kGray+1);
        lineMax->SetLineWidth(1);
        lineMax->Draw();

        // Legend on the first pad only
        if (i == 0) {
            TLegend *leg = new TLegend(0.55, 0.75, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(hCutUp[i],   "Spin#uparrow",   "l");
            leg->AddEntry(hCutDown[i], "Spin#downarrow", "l");
            leg->Draw();
        }
    }
    c->Update();
}

// =========================================================
// DrawElasInelComparison — elastic raw spectrum (black) vs
// inelastic raw spectrum (orange) with cut window lines, per
// scintillator. THIS object supplies the elastic spectrum;
// `inelHist` supplies the inelastic one. Both are filled from
// FULL statistics (no subsampling) — the raw MC sample sizes do
// NOT reflect the true physical elastic/inelastic cross-section
// ratio, so for display purposes only, each inelastic histogram
// is independently rescaled so that its integral relative to the
// elastic one equals the physical ratio R:
//
//   R > 0   →  hInel drawn scaled by R * Integral(hElas)/Integral(hInel)
//   R <= 0  →  hInel drawn unscaled (raw MC statistics)
//
// The underlying stored histograms are never modified — a scaled
// clone is drawn instead.
// =========================================================
void PolHistograms::DrawElasInelComparison(TCanvas *c, PolHistograms* inelHist,
                                            PolAnalysis* analysis, Double_t R) {
    if (!analysis) {
        std::cout << "Error: PolAnalysis pointer is null!" << std::endl;
        return;
    }
    if (!inelHist) {
        std::cout << "Error: inelastic PolHistograms pointer is null!" << std::endl;
        return;
    }

    Int_t nx = (nScintillators == 4) ? 4 : 3;
    c->Divide(nx, 2);

    for (Int_t i = 0; i < nScintillators * 2; i++) {
        c->cd(i + 1);
        gPad->SetLogy();   // inelastic tail is often orders of magnitude below elastic peak

        TH1F *hElas = hRaw[i];
        TH1F *hInel = inelHist->GetRawHist(i);
        if (!hElas || !hInel) continue;

        Bool_t scaled = false;
        if (R > 0.0 && hInel->Integral() > 0.0 && hElas->Integral() > 0.0) {
            Double_t scale = R * hElas->Integral() / hInel->Integral();
            hInel = (TH1F*)hInel->Clone(Form("%s_scaled", hInel->GetName()));
            hInel->SetDirectory(nullptr);
            hInel->Scale(scale);
            scaled = true;
        }

        hInel->SetLineColor(kOrange+7);
        hInel->SetLineWidth(2);
        hElas->SetLineColor(kBlack);
        hElas->SetLineWidth(1);

        Double_t ymax = TMath::Max(hElas->GetMaximum(), hInel->GetMaximum()) * 1.3;
        if (ymax <= 0) ymax = 1.0;
        hElas->SetMaximum(ymax);
        hElas->SetMinimum(0.5);   // log scale floor
        hElas->Draw("hist");
        hInel->Draw("hist same");

        Int_t scint = i % nScintillators;

        TLine *lineMin = new TLine(analysis->GetCutMin(scint), 0.5,
                                   analysis->GetCutMin(scint), ymax);
        lineMin->SetLineStyle(2);
        lineMin->SetLineColor(kBlue);
        lineMin->SetLineWidth(2);
        lineMin->Draw();

        TLine *lineMax = new TLine(analysis->GetCutMax(scint), 0.5,
                                   analysis->GetCutMax(scint), ymax);
        lineMax->SetLineStyle(2);
        lineMax->SetLineColor(kBlue);
        lineMax->SetLineWidth(2);
        lineMax->Draw();

        if (i == 0) {
            TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(hElas, "Elastic",     "l");
            leg->AddEntry(hInel, scaled ? "Inelastic (scaled to R)" : "Inelastic (raw stats)", "l");
            leg->AddEntry(lineMin, "Cut window", "l");
            leg->Draw();
        }
    }
    c->Update();
}


// =========================================================
// DrawElasInelCutComparison — same layout as DrawElasInelComparison
// but uses hCut (post-full-coincidence-cut) instead of hRaw.
// THIS object = elastic (black), inelHist = inelastic (orange).
// Cut lines are still drawn for reference — they confirm the
// window boundaries visually even though all events already passed.
// See DrawElasInelComparison for the meaning of R.
// =========================================================
void PolHistograms::DrawElasInelCutComparison(TCanvas *c, PolHistograms* inelHist,
                                               PolAnalysis* analysis, Double_t R) {
    if (!analysis || !inelHist) {
        std::cout << "Error: null pointer passed to DrawElasInelCutComparison" << std::endl;
        return;
    }

    Int_t nx = (nScintillators == 4) ? 4 : 3;
    c->Divide(nx, 2);

    for (Int_t i = 0; i < nScintillators * 2; i++) {
        c->cd(i + 1);
        gPad->SetLogy();

        TH1F *hElas = hCut[i];
        TH1F *hInel = inelHist->GetCutHist(i);
        if (!hElas || !hInel) continue;

        Bool_t scaled = false;
        if (R > 0.0 && hInel->Integral() > 0.0 && hElas->Integral() > 0.0) {
            Double_t scale = R * hElas->Integral() / hInel->Integral();
            hInel = (TH1F*)hInel->Clone(Form("%s_scaled", hInel->GetName()));
            hInel->SetDirectory(nullptr);
            hInel->Scale(scale);
            scaled = true;
        }

        hInel->SetLineColor(kOrange+7);
        hInel->SetLineWidth(2);
        hElas->SetLineColor(kBlack);
        hElas->SetLineWidth(1);

        Double_t ymax = TMath::Max(hElas->GetMaximum(), hInel->GetMaximum()) * 1.3;
        if (ymax <= 0) ymax = 1.0;
        hElas->SetMaximum(ymax);
        hElas->SetMinimum(0.5);
        hElas->Draw("hist");
        hInel->Draw("hist same");

        // Draw cut lines for reference
        Int_t scint = i % nScintillators;
        TLine *lineMin = new TLine(analysis->GetCutMin(scint), 0.5,
                                   analysis->GetCutMin(scint), ymax);
        lineMin->SetLineStyle(2); lineMin->SetLineColor(kBlue); lineMin->SetLineWidth(2);
        lineMin->Draw();
        TLine *lineMax = new TLine(analysis->GetCutMax(scint), 0.5,
                                   analysis->GetCutMax(scint), ymax);
        lineMax->SetLineStyle(2); lineMax->SetLineColor(kBlue); lineMax->SetLineWidth(2);
        lineMax->Draw();

        if (i == 0) {
            TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
            leg->SetBorderSize(0); leg->SetFillStyle(0);
            leg->AddEntry(hElas,   "Elastic (after cuts)",    "l");
            leg->AddEntry(hInel,   scaled ? "Inelastic (scaled to R, after cuts)"
                                           : "Inelastic (raw stats, after cuts)", "l");
            leg->AddEntry(lineMin, "Cut window",               "l");
            leg->Draw();
        }
    }
    c->Update();
}
#endif