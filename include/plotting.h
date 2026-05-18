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
        if (e%2 == 0){
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
        g->SetMarkerColor(colors[int(e/10)]);
        g->SetLineColor(colors[int(e/10)]);
        g->SetLineWidth(2);

        if(e == 0){
            g->SetTitle(";cos(#theta_{beam});W(#theta)/W(90^{#circ})");
            g->SetMinimum(0.8);
            g->SetMaximum(2.0);
            g->Draw("ALP");
        } else {
            g->Draw("LP SAME");
        }
        leg->AddEntry(g, Form("%.0f-%.0f MeV",
                              energy_bins[e], energy_bins[e+1]), "lp");
        }
    }
    leg->Draw();
    c->SaveAs(outname.c_str());
}

static void plotAnisotropyRatio(
    const std::vector<AnisotropyResult>& aniso,
    int nbins,
    int nbins_beam,
    const std::vector<double>& energy_bins,
    const std::string& outname)
{
    int bin_0 = nbins_beam - 1;

    std::vector<double> E_centers(nbins), ratio(nbins), u_ratio(nbins), ex(nbins, 0.0);
    for(int e = 0; e < nbins; ++e){
        E_centers[e] = std::sqrt(energy_bins[e] * energy_bins[e+1]);
        ratio[e]     = aniso[e].w[bin_0];
        u_ratio[e]   = aniso[e].u_w[bin_0];
    }

    // ── style ──────────────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(1);

    // fonts: 43 = helvetica, size in pixels
    gStyle->SetTextFont(43);
    gStyle->SetLabelFont(43, "XY");
    gStyle->SetTitleFont(43, "XY");
    gStyle->SetLabelSize(22, "XY");
    gStyle->SetTitleSize(26, "XY");

    gStyle->SetTickLength(0.025, "XY");
    gStyle->SetNdivisions(510, "X");
    gStyle->SetNdivisions(505, "Y");

    // ── canvas ─────────────────────────────────────────────────────────
    TCanvas* c = new TCanvas("c_ratio", "", 800, 650);
    c->SetLogx();
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.04);
    c->SetBottomMargin(0.13);

    // ── graph ──────────────────────────────────────────────────────────
    TGraphErrors* g = new TGraphErrors(
        nbins,
        E_centers.data(), ratio.data(),
        ex.data(),        u_ratio.data());

    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.9);
    g->SetMarkerColor(kAzure+2);
    g->SetLineColor(kAzure+2);
    g->SetLineWidth(1);

    g->GetXaxis()->SetTitle("Neutron energy  E_{n}  (MeV)");
    g->GetYaxis()->SetTitle("W(0^{#circ}) / W(90^{#circ})");
    g->GetXaxis()->SetTitleOffset(1.2);
    g->GetYaxis()->SetTitleOffset(1.3);
    g->GetYaxis()->SetRangeUser(0.5, 2.2);
    g->GetXaxis()->SetMoreLogLabels();
    g->GetXaxis()->SetNoExponent();

    g->Draw("AP");

    // ── isotropic reference ────────────────────────────────────────────
    TLine* line = new TLine(energy_bins.front(), 1.0, energy_bins.back(), 1.0);
    line->SetLineStyle(7);       // dashed
    line->SetLineColor(kGray+1);
    line->SetLineWidth(2);
    line->Draw();

    // ── shaded band ±5% around isotropy ───────────────────────────────
    int nband = 200;
    std::vector<double> xb(nband), yhi(nband), ylo(nband);
    double logA = std::log10(energy_bins.front());
    double logB = std::log10(energy_bins.back());
    for(int i = 0; i < nband; ++i){
        xb[i]  = std::pow(10.0, logA + (logB - logA) * i / (nband - 1));
        yhi[i] = 1.05;
        ylo[i] = 0.95;
    }
    TGraph* band = new TGraph(2 * nband);
    for(int i = 0;       i < nband; ++i) band->SetPoint(i,          xb[i],        yhi[i]);
    for(int i = 0;       i < nband; ++i) band->SetPoint(nband + i,  xb[nband-1-i], ylo[nband-1-i]);
    band->SetFillColorAlpha(kGray, 0.25);
    band->SetLineWidth(0);
    band->Draw("F SAME");

    // redraw graph on top of band
    g->Draw("P SAME");

    // ── legend ─────────────────────────────────────────────────────────
    TLegend* leg = new TLegend(0.55, 0.74, 0.93, 0.93);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(43);
    leg->SetTextSize(20);
    leg->AddEntry(g,    "^{238}U(n,f)  W(0^{#circ})/W(90^{#circ})", "lp");
    leg->AddEntry(line, "Isotropic",                                  "l");
    leg->AddEntry(band, "#pm5% band",                                 "f");
    leg->Draw();

    c->RedrawAxis();
    c->SaveAs(outname.c_str());
}