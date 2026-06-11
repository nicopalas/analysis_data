#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <iostream>

void flux_analysis() {

    // ── open files ────────────────────────────────────────────────────────────
    TFile *f_ear1 = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/flux_data/Energy_EAR1_FLUKA_500bdp.root",
        "READ");
    if (!f_ear1 || f_ear1->IsZombie()) {
        std::cerr << "[ERROR] Cannot open EAR1 file\n"; return;
    }

    TFile *f_eval = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/flux_data/evalFlux_prelim.root",
        "READ");
    if (!f_eval || f_eval->IsZombie()) {
        std::cerr << "[ERROR] Cannot open evalFlux file\n"; return;
    }

    // ── get histograms ────────────────────────────────────────────────────────
    TH1D *h_ear1 = (TH1D*)f_ear1->Get("h_Flux");
    if (!h_ear1) { std::cerr << "[ERROR] h_Flux not found\n"; return; }

    TH1D *h_eval = (TH1D*)f_eval->Get("hEval_Abs");
    if (!h_eval) { std::cerr << "[ERROR] hEval_Abs not found\n"; return; }

    // clone to avoid ownership issues
    h_ear1 = (TH1D*)h_ear1->Clone("h_ear1"); h_ear1->SetDirectory(nullptr);
    h_eval = (TH1D*)h_eval->Clone("h_eval"); h_eval->SetDirectory(nullptr);
    f_ear1->Close();
    f_eval->Close();

    // ── check and align binning ───────────────────────────────────────────────
    std::cout << "[INFO] EAR1  bins=" << h_ear1->GetNbinsX()
              << "  xmin=" << h_ear1->GetXaxis()->GetXmin()
              << "  xmax=" << h_ear1->GetXaxis()->GetXmax() << "\n";
    std::cout << "[INFO] eval  bins=" << h_eval->GetNbinsX()
              << "  xmin=" << h_eval->GetXaxis()->GetXmin()
              << "  xmax=" << h_eval->GetXaxis()->GetXmax() << "\n";

    // rebin eval onto EAR1 binning if they differ
    TH1D *h_eval_rebinned = nullptr;
    bool same_binning = (h_ear1->GetNbinsX()            == h_eval->GetNbinsX()  &&
                         h_ear1->GetXaxis()->GetXmin()   == h_eval->GetXaxis()->GetXmin() &&
                         h_ear1->GetXaxis()->GetXmax()   == h_eval->GetXaxis()->GetXmax());

    if (same_binning) {
        std::cout << "[INFO] Binning matches — dividing directly\n";
        h_eval_rebinned = h_eval;
    } else {
        std::cout << "[INFO] Binning differs — rebinning eval onto EAR1 grid\n";

        // build new histogram with EAR1 bin edges
        int    n   = h_ear1->GetNbinsX();
        double lo  = h_ear1->GetXaxis()->GetXmin();
        double hi  = h_ear1->GetXaxis()->GetXmax();

        // check if EAR1 uses variable bin widths
        if (h_ear1->GetXaxis()->IsVariableBinSize()) {
            const Double_t *edges = h_ear1->GetXaxis()->GetXbins()->GetArray();
            h_eval_rebinned = new TH1D("h_eval_rebinned", "", n, edges);
        } else {
            h_eval_rebinned = new TH1D("h_eval_rebinned", "", n, lo, hi);
        }
        h_eval_rebinned->SetDirectory(nullptr);

        // fill by integrating eval over each EAR1 bin
        for (int i = 1; i <= n; ++i) {
            double xlo = h_ear1->GetBinLowEdge(i);
            double xhi = h_ear1->GetBinLowEdge(i) + h_ear1->GetBinWidth(i);

            int b1 = h_eval->FindBin(xlo + 1e-9);
            int b2 = h_eval->FindBin(xhi - 1e-9);
            b1 = std::max(b1, 1);
            b2 = std::min(b2, h_eval->GetNbinsX());

            double val  = 0.0, err2 = 0.0;
            for (int b = b1; b <= b2; ++b) {
                // scale by overlap fraction if partial bin
                double blo     = h_eval->GetBinLowEdge(b);
                double bhi     = blo + h_eval->GetBinWidth(b);
                double overlap = (std::min(xhi, bhi) - std::max(xlo, blo))
                               / h_eval->GetBinWidth(b);
                overlap = std::max(overlap, 0.0);
                val  += h_eval->GetBinContent(b) * overlap;
                err2 += std::pow(h_eval->GetBinError(b) * overlap, 2);
            }
            h_eval_rebinned->SetBinContent(i, val);
            h_eval_rebinned->SetBinError  (i, std::sqrt(err2));
        }
    }

    // ── divide: eval / EAR1 ───────────────────────────────────────────────────
    TH1D *h_ratio = (TH1D*)h_eval_rebinned->Clone("h_ratio");
    h_ratio->SetDirectory(nullptr);
    h_ratio->SetTitle("evalFlux / EAR1;E (eV);ratio");

    for (int i = 1; i <= h_ratio->GetNbinsX(); ++i) {
        double num    = h_eval_rebinned->GetBinContent(i);
        double u_num  = h_eval_rebinned->GetBinError(i);
        double den    = h_ear1->GetBinContent(i);
        double u_den  = h_ear1->GetBinError(i);

        if (den <= 0.0 || num <= 0.0) {
            h_ratio->SetBinContent(i, 0.0);
            h_ratio->SetBinError  (i, 0.0);
            continue;
        }

        double ratio = num / den;
        double u_rel = std::sqrt((u_num/num)*(u_num/num) +
                                 (u_den/den)*(u_den/den));
        h_ratio->SetBinContent(i, ratio);
        h_ratio->SetBinError  (i, ratio * u_rel);
    }

    // ── save ──────────────────────────────────────────────────────────────────
    TFile *fout = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/flux_data/flux_ratio.root",
        "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[ERROR] Cannot create output file\n"; return;
    }

    h_ear1->Write("h_ear1");
    h_eval_rebinned->Write("h_eval_rebinned");
    h_ratio->Write("h_ratio");

    TCanvas *c = new TCanvas("c_ratio", "Flux ratio", 900, 600);
    c->SetLogx();
    h_ratio->SetLineColor(kBlue+1);
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerColor(kBlue+1);
    h_ratio->Draw("E");
    c->Write();

    fout->Close();
    std::cout << "[DONE] flux_ratio.root\n";
}