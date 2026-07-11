#pragma once
#include "types.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TColor.h"
#include "TLegend.h"
#include "TFile.h"
#include "TVirtualPad.h"
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

// ========================================================================
// House style for all publication figures in this analysis.
//
// One style function (setPubStyle) and one color source (okabeIto, the
// colorblind-safe 8-color qualitative palette standard in modern
// particle-physics papers) are shared by every plotting function below.
// Each figure gets its own accent color from that same palette so that,
// placed side by side in a thesis/note, they read as one consistent set
// rather than four different plotting styles glued together.
// ========================================================================

static int okabeIto(double r, double g, double b)
{
    return TColor::GetColor((Float_t)(r/255.), (Float_t)(g/255.), (Float_t)(b/255.));
}

// full qualitative palette (Okabe & Ito, 2008), CUD-safe, with sky blue and
// vermillion dropped per house choice — navy blue and burnt orange lead instead
static std::vector<int> okabeItoPalette()
{
    return {
        okabeIto(  0, 114, 178),  // navy blue      (primary accent 1)
        okabeIto(230, 159,   0),  // burnt orange   (primary accent 2)
        okabeIto(  0, 158, 115),  // bluish green
        okabeIto(204, 121, 167),  // reddish purple
        okabeIto(  0,   0,   0),  // black
        okabeIto(240, 228,  66)   // yellow (use sparingly — low contrast on white)
    };
}

// one accent color per figure, all drawn from the palette above
static const int kAnisoColor  = okabeIto(  0, 114, 178);  // navy blue      -> W(theta) panels
static const int kRatioColor  = okabeIto(128,   0,  32);  // burgundy       -> W(0)/W(90) vs E
static const int kBkgColor    = okabeIto(  0,   0,   0);  // black          -> raw TOF spectra
static const int kBkgFitColor = okabeIto(230, 159,   0);  // burnt orange   -> background fit
static const int kSubColor    = okabeIto(  0, 158, 115);  // bluish green   -> subtracted spectra

static void setPubStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetTickLength(0.025, "XY");
    gStyle->SetNdivisions(510, "X");
    gStyle->SetNdivisions(505, "Y");
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetTextFont(42);
    gStyle->SetLabelSize(0.045, "xyz");
    gStyle->SetTitleSize(0.050, "xyz");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendTextSize(0.032);
}

// grid layout that stays roughly square for a given number of pads
static void gridLayout(int n, int& ncol, int& nrow)
{
    ncol = (int)std::ceil(std::sqrt((double)n));
    nrow = (int)std::ceil((double)n / ncol);
}

static void stylePad(TVirtualPad* pad)
{
    pad->SetLeftMargin(0.16);
    pad->SetBottomMargin(0.14);
    pad->SetTopMargin(0.06);
    pad->SetRightMargin(0.04);
}

// ------------------------------------------------------------------------
// Background fits: one row per energy bin, raw spectrum + fit | subtracted
// spectrum. energy_bins is optional so existing call sites keep compiling;
// pass it to get per-row energy labels instead of a bare bin index.
// ------------------------------------------------------------------------
static void plotBackgroundFits(
    std::vector<TH1D*>& hists_tof,
    std::vector<TH1D*>& hists_sub,
    std::vector<TF1*>&  fits,
    int   nbins,
    const std::string& outname,
    const std::vector<double>& energy_bins = {})
{
    setPubStyle();

    TCanvas* c = new TCanvas("c_bkg", "Background fits", 1100, 380*nbins);
    c->Divide(2, nbins, 0.0002, 0.0002);

    const bool haveLabels = (int)energy_bins.size() == nbins + 1;

    for(int i = 0; i < nbins; ++i){

        std::string label = haveLabels
            ? Form("%.0f-%.0f MeV", energy_bins[i], energy_bins[i+1])
            : Form("bin %d", i);

        // --- left: raw spectrum + fit ---
        TVirtualPad* pL = c->cd(2*i + 1);
        stylePad(pL);
        pL->SetLogy();

        hists_tof[i]->SetTitle(";#Delta t (ns);counts");
        hists_tof[i]->SetLineColor(kBkgColor);
        hists_tof[i]->SetMarkerColor(kBkgColor);
        hists_tof[i]->SetMarkerStyle(20);
        hists_tof[i]->SetMarkerSize(0.2);
        hists_tof[i]->GetXaxis()->SetTitleOffset(1.1);
        hists_tof[i]->GetYaxis()->SetTitleOffset(1.4);
        hists_tof[i]->Draw("PE");

        fits[i]->SetLineColor(kBkgFitColor);
        fits[i]->SetLineWidth(2);
        fits[i]->Draw("SAME");

        TLatex latL; latL.SetNDC(); latL.SetTextFont(42);
        latL.SetTextSize(0.06); latL.SetTextAlign(33);
        latL.DrawLatex(0.94, 0.90, label.c_str());

        // --- right: background-subtracted spectrum ---
        TVirtualPad* pR = c->cd(2*i + 2);
        stylePad(pR);

        hists_sub[i]->SetTitle(";#Delta t (ns);counts (bkg. subtracted)");
        hists_sub[i]->SetLineColor(kSubColor);
        hists_sub[i]->SetMarkerColor(kSubColor);
        hists_sub[i]->SetMarkerStyle(20);
        hists_sub[i]->SetMarkerSize(0.2);
        hists_sub[i]->GetXaxis()->SetTitleOffset(1.1);
        hists_sub[i]->GetYaxis()->SetTitleOffset(1.4);
        hists_sub[i]->Draw("PE");

        TLine* zero = new TLine(hists_sub[i]->GetXaxis()->GetXmin(), 0.0,
                                hists_sub[i]->GetXaxis()->GetXmax(), 0.0);
        zero->SetLineStyle(2);
        zero->SetLineColor(kGray+1);
        zero->Draw("SAME");
    }
    c->SaveAs(outname.c_str());
}

// ------------------------------------------------------------------------
// Detector efficiency vs cos(theta'), one overlaid series per energy bin.
// Kept as an overlay (not per-panel) since these curves are meant to be
// compared directly; styled to match the rest of the house style.
// ------------------------------------------------------------------------
static void plotEfficiency(
    const std::vector<EfficiencyResult>& eff,
    int nbins,
    int nbins_det,
    const std::vector<double>& energy_bins,
    const std::string& outname)
{
    setPubStyle();

    std::vector<double> centers(nbins_det);
    for(int i = 0; i < nbins_det; ++i) centers[i] = (i + 0.5) * (1.0/nbins_det);

    std::vector<int> palette = okabeItoPalette();
    const int kSquareMarker = 21;  // filled square, per user request

    TCanvas* c = new TCanvas("c_eff", "Efficiency", 850, 680);
    stylePad(gPad);

    TLegend* leg = new TLegend(0.18, 0.68, 0.50, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    for(int e = 0; e < nbins; ++e){
        int color = palette[e % (int)palette.size()];

        TGraphErrors* gr = new TGraphErrors(
            nbins_det, centers.data(), eff[e].eps.data(),
            nullptr, eff[e].u_eps.data());
        gr->SetMarkerStyle(kSquareMarker);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color);
        gr->SetMarkerSize(1.);
        gr->SetLineWidth(2);

        if(e == 0){
            gr->SetTitle(";cos(#theta');#varepsilon(cos#theta')");
            gr->GetXaxis()->SetTitleOffset(1.1);
            gr->GetYaxis()->SetTitleOffset(1.4);
            gr->SetMinimum(0.0);
            gr->SetMaximum(1.25);
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

// ------------------------------------------------------------------------
// Relative resolution of the efficiency, sigma(eps)/eps, vs cos(theta'),
// one overlaid series per energy bin. Same layout and palette as
// plotEfficiency so the two figures read as a matched pair (value +
// uncertainty) rather than two unrelated plots. Log-y by default since
// the relative uncertainty typically spans more than a decade across the
// angular range (it blows up wherever eps -> 0, e.g. near cos(theta')=0).
// ------------------------------------------------------------------------
static void plotEfficiencyResolution(
    const std::vector<EfficiencyResult>& eff,
    int nbins,
    int nbins_det,
    const std::vector<double>& energy_bins,
    const std::string& outname,
    bool logy = true)
{
    setPubStyle();

    std::vector<double> centers(nbins_det);
    for(int i = 0; i < nbins_det; ++i) centers[i] = (i + 0.5) * (1.0/nbins_det);

    std::vector<int> palette = okabeItoPalette();
    // circles and squares only, filled then open, cycling — no triangles
    int markers[] = {20, 21, 24, 25};

    TCanvas* c = new TCanvas("c_eff_res", "Efficiency resolution", 850, 680);
    stylePad(gPad);
    if(logy) gPad->SetLogy();

    TLegend* leg = new TLegend(0.18, 0.68, 0.50, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    double ymax_seen = 0.0;

    std::vector<TGraph*> graphs;
    graphs.reserve(nbins);

    for(int e = 0; e < nbins; ++e){
        int color = palette[e % (int)palette.size()];

        // relative resolution: sigma(eps)/eps, only where eps > 0
        std::vector<double> x, y;
        x.reserve(nbins_det); y.reserve(nbins_det);
        for(int i = 0; i < nbins_det; ++i){
            double eps = eff[e].eps[i];
            if(eps <= 0.0) continue;
            double rel = eff[e].u_eps[i] / eps;
            x.push_back(centers[i]);
            y.push_back(rel);
            if(rel > ymax_seen) ymax_seen = rel;
        }

        TGraph* gr = new TGraph((int)x.size(), x.data(), y.data());
        graphs.push_back(gr);
        gr->SetMarkerStyle(markers[e % 4]);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color);
        gr->SetMarkerSize(1.0);
        gr->SetLineWidth(2);

        if(e == 0){
            gr->SetTitle(";cos(#theta');#sigma_{#varepsilon} / #varepsilon");
            gr->GetXaxis()->SetLimits(0.0, 1.0);
            gr->GetXaxis()->SetTitleOffset(1.1);
            gr->GetYaxis()->SetTitleOffset(1.4);
            gr->Draw("ALP");
        } else {
            gr->Draw("LP SAME");
        }
        leg->AddEntry(gr,
            Form("%.0f-%.0f MeV", energy_bins[e], energy_bins[e+1]), "lp");
    }

    // fix the y-range once all series are known, with headroom for the legend
    if(!graphs.empty()){
        double ymin = logy ? 1e-3 : 0.0;
        double ymax = logy ? ymax_seen * 3.0 : ymax_seen * 1.3;
        graphs[0]->SetMinimum(ymin);
        graphs[0]->SetMaximum(ymax);
        gPad->Modified();
        gPad->Update();
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
    setPubStyle();

    int ncol, nrow;
    gridLayout(nbins, ncol, nrow);

    TCanvas* c = new TCanvas("c_aniso", "Anisotropy", 480*ncol, 420*nrow);
    c->Divide(ncol, nrow, 0.0002, 0.0002);

    for(int e = 0; e < nbins; ++e){

        std::vector<double> x(nbins_beam), y(nbins_beam), ex(nbins_beam, 0.0);
        for(int i = 0; i < nbins_beam; ++i){
            x[i] = (i + 0.5) * dcos_beam;
            y[i] = aniso[e].w[i];
        }

        TGraphErrors* g = new TGraphErrors(
            nbins_beam, x.data(), y.data(), ex.data(), aniso[e].u_w.data());

        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.5);
        g->SetMarkerColor(kAnisoColor);
        g->SetLineColor(kAnisoColor);
        g->SetLineWidth(2);

        TVirtualPad* pad = c->cd(e + 1);
        stylePad(pad);
        pad->SetLeftMargin(0.18);
        pad->SetBottomMargin(0.16);

        g->SetTitle(";cos(#theta_{beam});W(#theta) / W(90^{#circ})");
        g->GetXaxis()->SetLimits(0.0, 1.0);
        g->SetMinimum(0.4);
        g->SetMaximum(2.5);
        g->GetXaxis()->SetTitleOffset(1.05);
        g->GetYaxis()->SetTitleOffset(1.35);
        g->Draw("AP");

        TLine* ref = new TLine(0.0, 1.0, 1.0, 1.0);
        ref->SetLineStyle(2);
        ref->SetLineColor(kGray + 1);
        ref->SetLineWidth(1);
        ref->Draw("SAME");

        TLatex lat; lat.SetNDC(); lat.SetTextFont(42);
        lat.SetTextSize(0.065); lat.SetTextAlign(33);
        lat.DrawLatex(0.94, 0.90,
                      Form("%.0f-%.0f MeV", energy_bins[e], energy_bins[e+1]));

        g->Draw("P SAME");
    }
    c->SaveAs(outname.c_str());
}

// ------------------------------------------------------------------------
// W(0)/W(90) vs neutron energy, single overlay figure.
// ------------------------------------------------------------------------
static void plotAnisotropyRatio(
    const std::vector<AnisotropyResult>& aniso,
    int nbins,
    int nbins_beam,
    const std::vector<double>& energy_bins,
    const std::string& outname)
{
    int bin_0 = nbins_beam - 1;

    std::vector<double> E_centers(nbins), ratio(nbins), u_ratio(nbins), ex(nbins);
    for(int e = 0; e < nbins; ++e){
        E_centers[e] = std::sqrt(energy_bins[e] * energy_bins[e+1]);
        ratio[e]     = aniso[e].w[bin_0];
        u_ratio[e]   = aniso[e].u_w[bin_0];
        // asymmetric in principle (log-spaced bins), TGraphErrors only takes
        // one value: use half the *linear* bin width as a simple width proxy
        ex[e]        = 0.5 * (energy_bins[e+1] - energy_bins[e]);
    }

    setPubStyle();

    TCanvas* c = new TCanvas("c_ratio", "", 800, 650);
    c->SetLogx();
    stylePad(gPad);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);

    TGraphErrors* g = new TGraphErrors(
        nbins, E_centers.data(), ratio.data(), ex.data(), u_ratio.data());

    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.);
    g->SetMarkerColor(kRatioColor);
    g->SetLineColor(kRatioColor);
    g->SetLineWidth(1);

    g->GetXaxis()->SetTitle("Neutron energy  E_{n}  (MeV)");
    g->GetYaxis()->SetTitle("W(0^{#circ}) / W(90^{#circ})");
    g->GetXaxis()->SetTitleOffset(1.2);
    g->GetYaxis()->SetTitleOffset(1.3);
    g->GetYaxis()->SetRangeUser(0.5, 2.2);
    g->GetXaxis()->SetMoreLogLabels();
    g->GetXaxis()->SetNoExponent();
    g->Draw("AP");

    TLine* line = new TLine(energy_bins.front(), 1.0, energy_bins.back(), 1.0);
    line->SetLineStyle(7);
    line->SetLineColor(kGray+1);
    line->SetLineWidth(2);
    line->Draw();

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

    g->Draw("P SAME");

    TLegend* leg = new TLegend(0.55, 0.74, 0.93, 0.93);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(g,    "^{238}U(n,f)  W(0^{#circ})/W(90^{#circ})", "lp");
    leg->AddEntry(line, "Isotropic",                                  "l");
    leg->AddEntry(band, "#pm5% band",                                 "f");
    leg->Draw();

    c->RedrawAxis();
    c->SaveAs((outname + ".pdf").c_str());

    TFile *fout = new TFile((outname + ".root").c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[ERROR] Cannot create output ROOT file: " << outname << ".root\n";
        return;
    }
    g->Write("anisotropy_ratio");
    line->Write("isotropic_line");
    band->Write("five_percent_band");
    fout->Close();

    std::cout << "[INFO] Saved: " << outname << ".pdf\n";
    std::cout << "[INFO] Saved: " << outname << ".root\n";
}

// ========================================================================
// Comparison of this-work anisotropy vs (any number of) EXFOR datasets.
//
// EXFOR CSV exports for DA (angular-distribution-ratio) quantities do not
// use a fixed column set: e.g. dataset 13709003 gives the y-error as
// "DATA-ERR (NO-DIM)", while 14660003 gives it as "ERR-T (NO-DIM)". The
// loader below matches by column-name PREFIX so it tolerates that, plus
// the resolution column being either a half-width (EN-RSL-HW) or a full
// width (EN-RSL) depending on the entry.
// ========================================================================

static const int kThisWorkColor = okabeIto(  0, 114, 178);  // navy blue — "this work", squares

// EXFOR accent colors: everything in the palette except navy (reserved
// for "this work") and yellow (too low-contrast to trust for data points)
static std::vector<int> exforPalette()
{
    return {
        okabeIto(128,   0,  32),  // burgundy
        okabeIto(  0, 158, 115),  // bluish green
        okabeIto(204, 121, 167),  // reddish purple
        okabeIto(230, 159,   0),  // burnt orange
        okabeIto(  0,   0,   0)   // black
    };
}

// minimal CSV line parser: handles double-quoted fields that themselves
// contain commas (EXFOR's Reacode column does exactly this)
static std::vector<std::string> parseCsvLine(const std::string& line)
{
    std::vector<std::string> fields;
    std::string cur;
    bool inQuotes = false;
    for(char c : line){
        if(c == '"'){ inQuotes = !inQuotes; continue; }
        if(c == ',' && !inQuotes){ fields.push_back(cur); cur.clear(); continue; }
        cur += c;
    }
    fields.push_back(cur);
    return fields;
}

struct ExforSource {
    std::string path;
    std::string label = "";   // leave empty to auto-build from author1/year1/DatasetID
};

struct ExforData {
    TGraphErrors* graph = nullptr;
    std::string   label;
};

// Loads one EXFOR DA-ratio CSV export. Energy EN (EV) -> MeV on X;
// DATA (NO-DIM) -> Y; best available *_ERR/*_ERR-T column -> Y error;
// EN-RSL(-HW) (EV) -> MeV, used as X half-width if present.
static ExforData loadExforAniso(const ExforSource& src)
{
    ExforData out;

    std::ifstream f(src.path);
    if(!f.is_open()){
        std::cerr << "[ERROR] Cannot open " << src.path << "\n";
        return out;
    }

    std::string header_line;
    std::getline(f, header_line);
    std::vector<std::string> header = parseCsvLine(header_line);

    auto findCol = [&](std::initializer_list<const char*> candidates) -> int {
        for(auto cand : candidates)
            for(size_t i = 0; i < header.size(); ++i)
                if(header[i].rfind(cand, 0) == 0)   // match by prefix
                    return (int)i;
        return -1;
    };

    int iID    = findCol({"DatasetID"});
    int iYear  = findCol({"year1"});
    int iAuth  = findCol({"author1"});
    int iEN    = findCol({"EN (EV)"});
    int iDATA  = findCol({"DATA (NO-DIM)"});
    // absolute (NO-DIM) error columns take priority; PER-CENT is a fallback
    // that must be converted to an absolute error using the DATA value
    int iERR_abs = findCol({"DATA-ERR (NO-DIM)", "ERR-T (NO-DIM)", "ERR (NO-DIM)"});
    int iERR_pct = findCol({"DATA-ERR (PER-CENT)", "ERR-T (PER-CENT)", "ERR (PER-CENT)"});
    int iERR     = (iERR_abs >= 0) ? iERR_abs : iERR_pct;
    bool errIsPercent = (iERR_abs < 0) && (iERR_pct >= 0);
    int iRSL   = findCol({"EN-RSL-HW (EV)", "EN-RSL (EV)"});
    bool rslIsHalfWidth = (iRSL >= 0) && header[iRSL].rfind("EN-RSL-HW", 0) == 0;

    if(iEN < 0 || iDATA < 0 || iERR < 0){
        std::cerr << "[ERROR] " << src.path << ": required columns not found\n";
        std::cerr << "        looked for EN (EV): " << (iEN>=0?"found":"MISSING")
                  << " | DATA (NO-DIM): " << (iDATA>=0?"found":"MISSING")
                  << " | error column: "  << (iERR>=0?"found":"MISSING") << "\n";
        std::cerr << "        header columns actually present in this file:\n";
        for(size_t i = 0; i < header.size(); ++i)
            std::cerr << "          [" << i << "] '" << header[i] << "'\n";
        return out;
    }

    std::vector<double> ex, ey, ex_err, ey_err;
    std::string line, ds_id, author, year;

    while(std::getline(f, line)){
        if(line.empty()) continue;
        auto fld = parseCsvLine(line);
        if((int)fld.size() <= std::max({iEN, iDATA, iERR})) continue;

        if(ds_id.empty() && iID >= 0 && (int)fld.size() > iID) ds_id  = fld[iID];
        if(author.empty() && iAuth >= 0 && (int)fld.size() > iAuth) author = fld[iAuth];
        if(year.empty()  && iYear >= 0 && (int)fld.size() > iYear) year   = fld[iYear];

        try{
            double E   = std::stod(fld[iEN]) * 1e-6;    // eV -> MeV
            double y   = std::stod(fld[iDATA]);
            double eyv = std::stod(fld[iERR]);
            if(errIsPercent) eyv = y * eyv / 100.0;     // % -> absolute
            double exv = 0.0;
            if(iRSL >= 0 && (int)fld.size() > iRSL && !fld[iRSL].empty()){
                double rsl = std::stod(fld[iRSL]) * 1e-6;
                exv = rslIsHalfWidth ? rsl : rsl / 2.0;
            }
            ex.push_back(E); ey.push_back(y);
            ex_err.push_back(exv); ey_err.push_back(eyv);
        } catch(...){ continue; }   // skip malformed/header-repeat rows
    }

    if(ex.empty()){
        std::cerr << "[WARN] " << src.path << ": no valid data rows parsed\n";
        return out;
    }

    out.graph = new TGraphErrors((int)ex.size(),
        ex.data(), ey.data(), ex_err.data(), ey_err.data());

    if(!src.label.empty()){
        out.label = src.label;
    } else {
        std::string entry = ds_id.size() >= 5 ? ds_id.substr(0, 5) : ds_id;
        out.label = (author.empty() ? "EXFOR" : author);
    }
    return out;
}

// Rough compatibility metric between "this work" and one EXFOR set:
// for each EXFOR point inside this-work's energy range, linearly
// interpolate this-work's value at that energy (in log E) and form a
// pull using the EXFOR point's own uncertainty. This is NOT a rigorous
// chi2 (this-work's own uncertainty at the interpolated point is not
// included), but it is a fast, honest-enough "how far off are we" number
// for a diagnostic panel label.
static bool chi2PerPointVsThis(TGraphErrors* ext, TGraphErrors* mine,
                               double& chi2_per_n, int& n_used)
{
    int nMine = mine->GetN();
    if(nMine < 2) return false;
    double* mx = mine->GetX();
    double* my = mine->GetY();

    double sum = 0.0; int n = 0;
    for(int i = 0; i < ext->GetN(); ++i){
        double x = ext->GetX()[i];
        if(x < mx[0] || x > mx[nMine-1]) continue;   // no extrapolation

        int k = (int)(std::lower_bound(mx, mx+nMine, x) - mx);
        if(k == 0) k = 1;
        double x0 = mx[k-1], x1 = mx[k], y0 = my[k-1], y1 = my[k];
        double t = (std::log(x) - std::log(x0)) / (std::log(x1) - std::log(x0));
        double yInterp = y0 + t * (y1 - y0);

        double ey = ext->GetEY()[i];
        if(ey <= 0.0) continue;
        double pull = (ext->GetY()[i] - yInterp) / ey;
        sum += pull*pull;
        ++n;
    }
    if(n == 0) return false;
    chi2_per_n = sum / n;
    n_used = n;
    return true;
}

// ------------------------------------------------------------------------
// Overlay: this work vs any number of EXFOR sources, single figure.
// ------------------------------------------------------------------------
static void plotAnisoVsExfor(
    TGraphErrors* g_this,
    const std::vector<ExforSource>& sources,
    const std::string& outname,
    double xmin = 1.0, double xmax = 1000.0)
{
    setPubStyle();

    std::vector<ExforData> exfor;
    for(auto& s : sources){
        ExforData d = loadExforAniso(s);
        if(d.graph) exfor.push_back(d);
    }
    if(!g_this){
        std::cerr << "[ERROR] plotAnisoVsExfor: this-work graph is null\n";
        return;
    }

    double ymin = 1e9, ymax = -1e9;
    auto updateRange = [&](TGraphErrors* g){
        for(int i = 0; i < g->GetN(); ++i){
            double v = g->GetY()[i], e = g->GetEY()[i];
            ymin = std::min(ymin, v - e);
            ymax = std::max(ymax, v + e);
        }
    };
    updateRange(g_this);
    for(auto& d : exfor) updateRange(d.graph);
    double margin = (ymax - ymin) * 0.15;
    ymin -= margin; ymax += margin;

    TCanvas* c = new TCanvas("c_aniso_exfor", "Anisotropy vs EXFOR", 900, 650);
    c->SetLogx();
    stylePad(gPad);

    TH1F* frame = c->DrawFrame(xmin, ymin, xmax, ymax);
    frame->GetXaxis()->SetTitle("E_{n} (MeV)");
    frame->GetYaxis()->SetTitle("W(0^{#circ}) / W(90^{#circ})");
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetYaxis()->SetTitleOffset(1.3);
    frame->GetXaxis()->SetMoreLogLabels();
    frame->GetXaxis()->SetNoExponent();

    TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray+1);
    line->SetLineWidth(1);
    line->Draw();

    TLegend* leg = new TLegend(0.55, 0.72 - 0.045*exfor.size(), 0.93, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    // EXFOR sources: open circles, one color each, drawn first (background)
    std::vector<int> pal = exforPalette();
    for(size_t k = 0; k < exfor.size(); ++k){
        int color = pal[k % pal.size()];
        exfor[k].graph->SetMarkerStyle(24);   // open circle
        exfor[k].graph->SetMarkerSize(1.);
        exfor[k].graph->SetMarkerColor(color);
        exfor[k].graph->SetLineColor(color);
        exfor[k].graph->SetLineWidth(1);
        exfor[k].graph->Draw("P SAME");
        leg->AddEntry(exfor[k].graph, exfor[k].label.c_str(), "lp");
    }

    // this work: filled navy square, drawn last (foreground)
    g_this->SetMarkerStyle(21);
    g_this->SetMarkerSize(1.);
    g_this->SetMarkerColor(kThisWorkColor);
    g_this->SetLineColor(kThisWorkColor);
    g_this->SetLineWidth(2);
    g_this->Draw("P SAME");
    leg->AddEntry(g_this, "This work", "lp");

    leg->AddEntry(line, "Isotropic", "l");
    leg->Draw();
    c->RedrawAxis();
    c->SaveAs(outname.c_str());
}


// ------------------------------------------------------------------------
// Per-point pull: (y_ext - y_this_interp) / sqrt(sigma_ext^2 + sigma_this_interp^2)
//
// Both this-work's value AND its uncertainty are linearly interpolated
// (in log E) to the EXFOR point's energy. Interpolating an uncertainty
// linearly is an approximation (the true point-to-point correlation of
// this-work's own systematic effects is not modeled), but it is the
// standard practical choice for this kind of compatibility check, and
// is far more honest than ignoring this-work's error entirely.
//
// A well-matched dataset should show pulls scattered around 0 with an
// RMS near 1; a systematically low/high or over/under-dispersed pull
// distribution flags disagreement or an underestimated uncertainty.
// ------------------------------------------------------------------------
struct PullSeries {
    TGraph* graph = nullptr;   // x = E_n (MeV), y = pull
    double  mean  = 0.0;
    double  rms   = 0.0;
    int     n     = 0;
};
 
static PullSeries computePulls(TGraphErrors* ext, TGraphErrors* mine)
{
    PullSeries out;
    int nMine = mine->GetN();
    if(nMine < 2) return out;
 
    double* mx  = mine->GetX();
    double* my  = mine->GetY();
    double* mey = mine->GetEY();
 
    std::vector<double> px, py;
    double sum = 0.0, sum2 = 0.0;
 
    for(int i = 0; i < ext->GetN(); ++i){
        double x = ext->GetX()[i];
        if(x < mx[0] || x > mx[nMine-1]) continue;   // no extrapolation
 
        int k = (int)(std::lower_bound(mx, mx+nMine, x) - mx);
        if(k == 0) k = 1;
        double x0 = mx[k-1], x1 = mx[k];
        double t  = (std::log(x) - std::log(x0)) / (std::log(x1) - std::log(x0));
 
        double yInterp  = my[k-1]  + t * (my[k]  - my[k-1]);
        double eyInterp = mey[k-1] + t * (mey[k] - mey[k-1]);
 
        double ey_ext = ext->GetEY()[i];
        double sigma  = std::sqrt(ey_ext*ey_ext + eyInterp*eyInterp);
        if(sigma <= 0.0) continue;
 
        double pull = (ext->GetY()[i] - yInterp) / sigma;
        px.push_back(x);
        py.push_back(pull);
        sum  += pull;
        sum2 += pull*pull;
    }
 
    if(px.empty()) return out;
 
    out.n     = (int)px.size();
    out.mean  = sum / out.n;
    out.rms   = std::sqrt(sum2 / out.n);
    out.graph = new TGraph(out.n, px.data(), py.data());
    return out;
}
 
// ------------------------------------------------------------------------
// Grid of pull plots, one panel per EXFOR source, to see which dataset
// is best matched (pulls near 0, RMS near 1) without cross-clutter.
// ------------------------------------------------------------------------
static void plotPullsVsExfor(
    TGraphErrors* g_this,
    const std::vector<ExforSource>& sources,
    const std::string& outname,
    double xmin = 1.1, double xmax = 1000.0)
{
    setPubStyle();
 
    std::vector<ExforData> exfor;
    for(auto& s : sources){
        ExforData d = loadExforAniso(s);
        if(d.graph) exfor.push_back(d);
    }
    if(exfor.empty() || !g_this) return;
 
    int ncol, nrow;
    gridLayout((int)exfor.size(), ncol, nrow);
 
    TCanvas* c = new TCanvas("c_pulls", "Pulls vs EXFOR", 480*ncol, 420*nrow);
    c->Divide(ncol, nrow, 0.0002, 0.0002);
 
    std::vector<int> pal = exforPalette();
 
    for(size_t k = 0; k < exfor.size(); ++k){
        int color = pal[k % pal.size()];
 
        PullSeries p = computePulls(exfor[k].graph, g_this);
 
        TVirtualPad* pad = c->cd((int)k + 1);
        stylePad(pad);
        pad->SetLogx();
 
        double ymax = 4.0;
        if(p.graph){
            for(int i = 0; i < p.graph->GetN(); ++i)
                ymax = std::max(ymax, std::abs(p.graph->GetY()[i]) * 1.2);
        }
 
        TH1F* frame = pad->DrawFrame(xmin, -ymax, xmax, ymax);
        frame->GetXaxis()->SetTitle("E_{n} (MeV)");
        frame->GetYaxis()->SetTitle("pull  = (data - this work) / #sigma");
        frame->GetXaxis()->SetTitleOffset(1.15);
        frame->GetYaxis()->SetTitleOffset(1.3);
        frame->GetXaxis()->SetMoreLogLabels();
        frame->GetXaxis()->SetNoExponent();
 
        TLine* l0 = new TLine(xmin, 0.0, xmax, 0.0);
        l0->SetLineColor(kGray+2); l0->SetLineWidth(1);
        l0->Draw();
        TLine* lp1 = new TLine(xmin,  1.0, xmax,  1.0);
        TLine* lm1 = new TLine(xmin, -1.0, xmax, -1.0);
        lp1->SetLineStyle(2); lp1->SetLineColor(kGray+1);
        lm1->SetLineStyle(2); lm1->SetLineColor(kGray+1);
        lp1->Draw(); lm1->Draw();
        TLine* lp2 = new TLine(xmin,  2.0, xmax,  2.0);
        TLine* lm2 = new TLine(xmin, -2.0, xmax, -2.0);
        lp2->SetLineStyle(3); lp2->SetLineColor(kGray);
        lm2->SetLineStyle(3); lm2->SetLineColor(kGray);
        lp2->Draw(); lm2->Draw();
 
        if(p.graph){
            p.graph->SetMarkerStyle(24);
            p.graph->SetMarkerSize(0.5);
            p.graph->SetMarkerColor(color);
            p.graph->SetLineColor(color);
            p.graph->Draw("P SAME");
        }
 
        TLatex lat; lat.SetNDC(); lat.SetTextFont(42);
        lat.SetTextSize(0.055); lat.SetTextAlign(13);
        lat.DrawLatex(0.06, 0.92, exfor[k].label.c_str());
 
        if(p.graph){
            lat.SetTextAlign(33);
            lat.SetTextSize(0.048);
            lat.DrawLatex(0.94, 0.92,
                Form("#LTpull#GT=%.2f  RMS=%.2f  (N=%d)", p.mean, p.rms, p.n));
        }
    }
    c->SaveAs(outname.c_str());
}
 
// ------------------------------------------------------------------------
// Overlay of all sources' pulls in one figure, for a direct side-by-side
// read of which dataset sits closest to 0 with the tightest scatter.
// ------------------------------------------------------------------------
static void plotPullsOverlay(
    TGraphErrors* g_this,
    const std::vector<ExforSource>& sources,
    const std::string& outname,
    double xmin = 1.1, double xmax = 1000.0)
{
    setPubStyle();
 
    std::vector<ExforData> exfor;
    for(auto& s : sources){
        ExforData d = loadExforAniso(s);
        if(d.graph) exfor.push_back(d);
    }
    if(exfor.empty() || !g_this) return;
 
    std::vector<PullSeries> pulls;
    double ymax = 4.0;
    for(auto& d : exfor){
        PullSeries p = computePulls(d.graph, g_this);
        if(p.graph)
            for(int i = 0; i < p.graph->GetN(); ++i)
                ymax = std::max(ymax, std::abs(p.graph->GetY()[i]) * 1.2);
        pulls.push_back(p);
    }
 
    TCanvas* c = new TCanvas("c_pulls_overlay", "Pulls overlay", 900, 650);
    c->SetLogx();
    stylePad(gPad);
 
    TH1F* frame = c->DrawFrame(xmin, -ymax, xmax, ymax);
    frame->GetXaxis()->SetTitle("E_{n} (MeV)");
    frame->GetYaxis()->SetTitle("pull");
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetYaxis()->SetTitleOffset(1.3);
    frame->GetXaxis()->SetMoreLogLabels();
    frame->GetXaxis()->SetNoExponent();
 
    TLine* l0 = new TLine(xmin, 0.0, xmax, 0.0);
    l0->SetLineColor(kGray+2); l0->SetLineWidth(1); l0->Draw();
    TLine* lp1 = new TLine(xmin, 1.0, xmax, 1.0);
    TLine* lm1 = new TLine(xmin,-1.0, xmax,-1.0);
    lp1->SetLineStyle(2); lp1->SetLineColor(kGray+1);
    lm1->SetLineStyle(2); lm1->SetLineColor(kGray+1);
    lp1->Draw(); lm1->Draw();
 
    std::vector<int> pal = exforPalette();
    TLegend* leg = new TLegend(0.55, 0.72, 0.93, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
 
    for(size_t k = 0; k < exfor.size(); ++k){
        if(!pulls[k].graph) continue;
        int color = pal[k % pal.size()];
        pulls[k].graph->SetMarkerStyle(21 + (int)(k % 2));  // open circle / open square
        pulls[k].graph->SetMarkerSize(0.4);
        pulls[k].graph->SetMarkerColor(color);
        pulls[k].graph->SetLineColor(color);
        pulls[k].graph->Draw("P SAME");
        leg->AddEntry(pulls[k].graph,
            Form("%s  (#LTpull#GT=%.2f, RMS=%.2f)",
                 exfor[k].label.c_str(), pulls[k].mean, pulls[k].rms), "p");
    }
    leg->Draw();
    c->RedrawAxis();
    c->SaveAs(outname.c_str());
}
static void plotAnisoVsExforIndividual(
    TGraphErrors* g_this,
    const std::vector<ExforSource>& sources,
    const std::string& outname,
    double xmin = 1.1, double xmax = 1000.0)
{
    setPubStyle();
 
    std::vector<ExforData> exfor;
    for(auto& s : sources){
        ExforData d = loadExforAniso(s);
        if(d.graph) exfor.push_back(d);
    }
    if(exfor.empty() || !g_this) return;
 
    int ncol, nrow;
    gridLayout((int)exfor.size(), ncol, nrow);
 
    TCanvas* c = new TCanvas("c_aniso_exfor_grid", "This work vs each EXFOR set",
                             480*ncol, 420*nrow);
    c->Divide(ncol, nrow, 0.0002, 0.0002);
 
    std::vector<int> pal = exforPalette();
 
    for(size_t k = 0; k < exfor.size(); ++k){
        int color = pal[k % pal.size()];
 
        double ymin = 1e9, ymax = -1e9;
        auto updateRange = [&](TGraphErrors* g){
            for(int i = 0; i < g->GetN(); ++i){
                double v = g->GetY()[i], e = g->GetEY()[i];
                ymin = std::min(ymin, v - e);
                ymax = std::max(ymax, v + e);
            }
        };
        updateRange(g_this);
        updateRange(exfor[k].graph);
        double margin = (ymax - ymin) * 0.20;
        ymin -= margin; ymax += margin;
 
        TVirtualPad* pad = c->cd((int)k + 1);
        stylePad(pad);
        pad->SetLogx();
 
        TH1F* frame = pad->DrawFrame(xmin, ymin, xmax, ymax);
        frame->GetXaxis()->SetTitle("E_{n} (MeV)");
        frame->GetYaxis()->SetTitle("W(0^{#circ})/W(90^{#circ})");
        frame->GetXaxis()->SetTitleOffset(1.15);
        frame->GetYaxis()->SetTitleOffset(1.35);
        frame->GetXaxis()->SetMoreLogLabels();
        frame->GetXaxis()->SetNoExponent();
 
        TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
        line->SetLineStyle(2);
        line->SetLineColor(kGray+1);
        line->Draw();
 
        exfor[k].graph->SetMarkerStyle(24);
        exfor[k].graph->SetMarkerSize(0.4);
        exfor[k].graph->SetMarkerColor(color);
        exfor[k].graph->SetLineColor(color);
        exfor[k].graph->Draw("P SAME");
 
        g_this->SetMarkerStyle(21);
        g_this->SetMarkerSize(0.2);
        g_this->SetMarkerColor(kThisWorkColor);
        g_this->SetLineColor(kThisWorkColor);
        g_this->Draw("P SAME");
 
        TLatex lat; lat.SetNDC(); lat.SetTextFont(42);
        lat.SetTextSize(0.055); lat.SetTextAlign(13);
        lat.DrawLatex(0.06, 0.92, exfor[k].label.c_str());
 
        double chi2n; int nUsed;
        if(chi2PerPointVsThis(exfor[k].graph, g_this, chi2n, nUsed)){
            lat.SetTextAlign(33);
            lat.SetTextSize(0.050);
            lat.DrawLatex(0.94, 0.92,
                Form("#chi^{2}/N #approx %.2f  (N=%d)", chi2n, nUsed));
        }
    }
    c->SaveAs(outname.c_str());
}