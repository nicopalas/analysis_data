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

static void plotAnisotropy(
    const std::vector<AnisotropyResult>& aniso,
    int nbins,
    int nbins_beam,
    const std::vector<double>& energy_bins,
    const std::string& outname)
{
    TCanvas* c = new TCanvas("c_aniso", "Anisotropy", 900, 700);
    gStyle->SetOptStat(0);
    c->SetLeftMargin(0.15);

    int colors[] = {kOrange, kSpring+5, kAzure+3, kMagenta+2};
    TLegend* leg = new TLegend(0.15, 0.70, 0.55, 0.90);
    leg->SetBorderSize(0);

    for(int e = 0; e < nbins; ++e){
        std::vector<double> x(nbins_beam), y(nbins_beam), ex(nbins_beam, 0.0);
        for(int i = 0; i < nbins_beam; ++i){
            x[i] = (i + 0.5) * dcos_beam;
            y[i] = aniso[e].w[i];
        }
        TGraphErrors* g = new TGraphErrors(
            nbins_beam,
            x.data(), y.data(),
            ex.data(), aniso[e].u_w.data());

        g->SetMarkerStyle(20 + e);
        g->SetMarkerColor(colors[e]);
        g->SetLineColor(colors[e]);
        g->SetLineWidth(2);

        if(e == 0){
            g->SetTitle(";cos(#theta_{beam});W(#theta)/W(90^{#circ})");
            g->GetYaxis()->SetRangeUser(0.5, 1.5);
            g->Draw("ALP");
        } else {
            g->Draw("LP SAME");
        }
        leg->AddEntry(g, Form("%.0f-%.0f MeV",
                              energy_bins[e], energy_bins[e+1]), "lp");
    }
    leg->Draw();
    c->SaveAs(outname.c_str());
}

static void plotAnisotropyRatio(
    const std::vector<AnisotropyResult>& aniso,
    int nbins,
    const std::vector<double>& energy_bins,
    const std::string& outname)
{
    // W0/W90 is just aniso[e].w[0] since everything is normalized to W(90)=1
    // W0 is the last beam bin (cos_theta ~ 1, i.e. 0 degrees)
    int bin_0  = 0;              // cos_theta ~ 0 → 90 degrees from beam
    int bin_90 = nbins_beam - 1; // cos_theta ~ 1 → 0 degrees (beam direction)

    std::vector<double> E_centers(nbins), ratio(nbins), u_ratio(nbins), ex(nbins, 0.0);
    for(int e = 0; e < nbins; ++e){
        E_centers[e] = std::sqrt(energy_bins[e] * energy_bins[e+1]);
        ratio[e]     = aniso[e].w[bin_90];    // W(0 deg) normalized to W(90)=1
        u_ratio[e]   = aniso[e].u_w[bin_90];
    }

    TGraphErrors* g = new TGraphErrors(
        nbins,
        E_centers.data(), ratio.data(),
        ex.data(),        u_ratio.data());

    g->SetTitle(";E_{n} (MeV);W_{0^{#circ}}/W_{90^{#circ}}");
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.2);
    g->SetMarkerColor(kAzure+3);
    g->SetLineColor(kAzure+3);
    g->SetLineWidth(2);

    TCanvas* c = new TCanvas("c_ratio", "W0/W90 vs Energy", 900, 700);
    gStyle->SetOptStat(0);
    c->SetLogx();
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    gPad->SetGrid(1, 1);

    // reference line at 1
    g->GetYaxis()->SetRangeUser(0.5, 2.0);
    g->Draw("APL");

    TLine* line = new TLine(energy_bins.front(), 1.0, energy_bins.back(), 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray+1);
    line->SetLineWidth(1);
    line->Draw();

    TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(g,    "^{238}U  W_{0}/W_{90}", "lp");
    leg->AddEntry(line, "Isotropic (=1)",         "l");
    leg->Draw();

    c->SaveAs(outname.c_str());
}