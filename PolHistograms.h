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
    std::vector<TH2D*> h2D;
    
    Int_t nScintillators;
    
public:
    PolHistograms(TString _name, TString _title, Int_t _nScint, Int_t _nbins, Double_t _emin, Double_t _emax) {
        name = _name;
        title = _title;
        nScintillators = _nScint;
        nbins = _nbins;
        emin = _emin;
        emax = _emax;
        
        // Create raw histograms
        for (Int_t i = 0; i < nScintillators * 2; i++) {
            TString side = (i < nScintillators) ? "L" : "R";
            Int_t scint = i % nScintillators;
            
            TH1F *h = new TH1F(Form("h%s%sS%d", name.Data(), side.Data(), scint),
                              Form("%s %s S%d;Energy [MeV];Counts", title.Data(), side.Data(), scint),
                              nbins, emin, emax);
            hRaw.push_back(h);
            
            TH1F *hc = new TH1F(Form("h%s%sS%d_cut", name.Data(), side.Data(), scint),
                               Form("%s %s S%d (with cuts);Energy [MeV];Counts", title.Data(), side.Data(), scint),
                               nbins, emin, emax);
            hc->SetLineColor(kRed);
            hCut.push_back(hc);
        }
    }
    
    void Fill(Int_t idx, Double_t *energies, Bool_t applyNonZero = true) {
        if (applyNonZero && energies[idx] <= 0) return;
        hRaw[idx]->Fill(energies[idx]);
    }
    
    void FillWithCuts(Int_t idx, Double_t *energies, PolAnalysis* analysis);
    
    void Create2DHistograms() {
        if (nScintillators == 4) {
            // Left arm 2D correlations
            h2D.push_back(new TH2D(Form("h%sLS1S0", name.Data()), 
                                  Form("%s Left S1 vs S0;S0 Energy [MeV];S1 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sLS2S1", name.Data()), 
                                  Form("%s Left S2 vs S1;S1 Energy [MeV];S2 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sLS3S2", name.Data()), 
                                  Form("%s Left S3 vs S2;S2 Energy [MeV];S3 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sLS3S1", name.Data()), 
                                  Form("%s Left S3 vs S1;S1 Energy [MeV];S3 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            
            // Right arm 2D correlations
            h2D.push_back(new TH2D(Form("h%sRS1S0", name.Data()), 
                                  Form("%s Right S1 vs S0;S0 Energy [MeV];S1 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS2S1", name.Data()), 
                                  Form("%s Right S2 vs S1;S1 Energy [MeV];S2 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS3S2", name.Data()), 
                                  Form("%s Right S3 vs S2;S2 Energy [MeV];S3 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS3S1", name.Data()), 
                                  Form("%s Right S3 vs S1;S1 Energy [MeV];S3 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
        } else if (nScintillators == 3) {
            // Left arm 2D correlations for 3 scintillators
            h2D.push_back(new TH2D(Form("h%sLS1S0", name.Data()), 
                                  Form("%s Left S1 vs S0;S0 Energy [MeV];S1 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sLS2S1", name.Data()), 
                                  Form("%s Left S2 vs S1;S1 Energy [MeV];S2 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sLS2S0", name.Data()), 
                                  Form("%s Left S2 vs S0;S0 Energy [MeV];S2 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sDummy1", name.Data()), 
                                  "Dummy", nbins, emin, emax, nbins, emin, emax));
            
            // Right arm 2D correlations for 3 scintillators
            h2D.push_back(new TH2D(Form("h%sRS1S0", name.Data()), 
                                  Form("%s Right S1 vs S0;S0 Energy [MeV];S1 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS2S1", name.Data()), 
                                  Form("%s Right S2 vs S1;S1 Energy [MeV];S2 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sRS2S0", name.Data()), 
                                  Form("%s Right S2 vs S0;S0 Energy [MeV];S2 Energy [MeV]", title.Data()), 
                                  nbins, emin, emax, nbins, emin, emax));
            h2D.push_back(new TH2D(Form("h%sDummy2", name.Data()), 
                                  "Dummy", nbins, emin, emax, nbins, emin, emax));
        }
    }
    
    void Fill2D(Double_t *energies) {
        Double_t threshold = 0.1;
        
        if (nScintillators == 4) {
            // Left arm correlations
            if (energies[0] > threshold && energies[1] > threshold) h2D[0]->Fill(energies[0], energies[1]);
            if (energies[1] > threshold && energies[2] > threshold) h2D[1]->Fill(energies[1], energies[2]);
            if (energies[2] > threshold && energies[3] > threshold) h2D[2]->Fill(energies[2], energies[3]);
            if (energies[1] > threshold && energies[3] > threshold) h2D[3]->Fill(energies[1], energies[3]);
            
            // Right arm correlations
            if (energies[4] > threshold && energies[5] > threshold) h2D[4]->Fill(energies[4], energies[5]);
            if (energies[5] > threshold && energies[6] > threshold) h2D[5]->Fill(energies[5], energies[6]);
            if (energies[6] > threshold && energies[7] > threshold) h2D[6]->Fill(energies[6], energies[7]);
            if (energies[5] > threshold && energies[7] > threshold) h2D[7]->Fill(energies[5], energies[7]);
        } else if (nScintillators == 3) {
            // Left arm correlations
            if (energies[0] > threshold && energies[1] > threshold) h2D[0]->Fill(energies[0], energies[1]);
            if (energies[1] > threshold && energies[2] > threshold) h2D[1]->Fill(energies[1], energies[2]);
            if (energies[0] > threshold && energies[2] > threshold) h2D[2]->Fill(energies[0], energies[2]);
            
            // Right arm correlations
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

// Now implement methods that depend on PolAnalysis
void PolHistograms::FillWithCuts(Int_t idx, Double_t *energies, PolAnalysis* analysis) {
    if (!analysis) return;
    if (energies[idx] <= 0) return;
    
    Int_t scint = idx % nScintillators;
    Bool_t passCuts = true;
    
    Bool_t isLeft = (idx < nScintillators);
    Int_t offset = isLeft ? 0 : nScintillators;
    
    for (Int_t i = 0; i < nScintillators; i++) {
        if (i == scint) continue;
        Double_t e = energies[offset + i];
        if (e < analysis->GetCutMin(i) || e > analysis->GetCutMax(i)) {
            passCuts = false;
            break;
        }
    }
    
    if (passCuts) {
        hCut[idx]->Fill(energies[idx]);
    }
}

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
        
        Int_t scint = i % nScintillators;
        Double_t ypos = 0.85;
        for (Int_t j = 0; j < nScintillators; j++) {
            if (j == scint) continue;
            latex.DrawLatexNDC(0.15, ypos, Form("S%d: %.1f-%.1f", j, 
                              analysis->GetCutMin(j), analysis->GetCutMax(j)));
            ypos -= 0.05;
        }
    }
    c->Update();
}

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
        
        // Draw vertical lines ONLY for the current scintillator's cuts
        Double_t ymax = hRaw[i]->GetMaximum() * 1.1;
        
        TLine *lineMin = new TLine(analysis->GetCutMin(scint), 0, analysis->GetCutMin(scint), ymax);
        lineMin->SetLineStyle(2);
        lineMin->SetLineColor(kBlue);
        lineMin->SetLineWidth(2);
        lineMin->Draw();
        
        TLine *lineMax = new TLine(analysis->GetCutMax(scint), 0, analysis->GetCutMax(scint), ymax);
        lineMax->SetLineStyle(2);
        lineMax->SetLineColor(kBlue);
        lineMax->SetLineWidth(2);
        lineMax->Draw();
        
        if (i == 0) {
            TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
            leg->AddEntry(hRaw[i], "Raw", "l");
            leg->AddEntry(hCut[i], "With Cuts", "l");
            leg->Draw();
        }
    }
    c->Update();
}

#endif