#pragma once
#include "types.h" 
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TF1.h"
#include <vector>
#include <string>

static void plotBackgroundFits(
    std::vector<TH1D*>& hists_tof,
    std::vector<TH1D*>& hists_sub,
    std::vector<TF1*>& fits,
    int   nbins,
    const std::string& outname)
{
    TCanvas* c = new TCanvas("c_bkg", "Background fits", 1200, 800);
    c->Divide(2, nbins);
    for(int i = 0; i < nbins; ++i){
        c->cd(2*i + 1); gPad->SetLogy();
        hists_tof[i]->Draw();
        fits[i]->SetLineColor(kRed);
        fits[i]->Draw("same");
        c->cd(2*i + 2); gPad->SetLogy();
        hists_sub[i]->Draw();
    }
    c->SaveAs(outname.c_str());
}

static void plotEfficiency(
    const std::vector<EfficiencyResult>& eff,
    int nbins,
    int nbins_det,
    const std::vector<double>& energy_bins,
    const std::string& outname)
{
    double centers[10];
    for(int i = 0; i < nbins_det; ++i) centers[i] = (i + 0.5) * (1.0/nbins_det);

    int colors[] = {kOrange, kSpring+5, kAzure+3, kMagenta+2};

    TCanvas* c = new TCanvas("c_eff", "Efficiency", 900, 700);
    gStyle->SetOptStat(0);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.15);

    TLegend* leg = new TLegend(0.15, 0.70, 0.45, 0.90);
    leg->SetBorderSize(0);

    for(int e = 0; e < nbins; ++e){
        TGraphErrors* gr = new TGraphErrors(
            nbins_det, centers, eff[e].eps.data(), nullptr, eff[e].u_eps.data());
        gr->SetMarkerStyle(21);
        gr->SetMarkerColor(colors[e]);
        gr->SetLineColor(colors[e]);
        gr->SetLineWidth(2);

        if(e == 0){
            gr->SetTitle("");
            gr->GetXaxis()->SetTitle("cos(#theta')");
            gr->GetYaxis()->SetTitle("#varepsilon(cos#theta')");
            gr->GetYaxis()->SetRangeUser(0.0, 1.2);
            gr->Draw("ALP");
        } else {
            gr->Draw("LP SAME");
        }
        leg->AddEntry(gr,
            Form("%.0f-%.0f MeV", energy_bins[e], energy_bins[e+1]), "lp");
    }
    leg->Draw();
    c->SaveAs(outname.c_str());
}