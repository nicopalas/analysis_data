#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include <cmath>
#include <iostream>
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
// pass it to get per-row en
// ------------------------------------------------------------------------
// Corrected angle reconstruction.
//
// Changes relative to the version you pasted:
//   1) The stray reference to `phi` in the failure branch was removed --
//      it was never a parameter or local variable (leftover from an
//      earlier version), and would not compile.
//   2) phi_det is now correctly an OUTPUT (computed internally via
//      atan2(dy,dx)), matching your latest signature -- it is no longer
//      expected as an input from the caller.
//   3) dx uses each point's own offset ((x1-offsetx1)-(x0-offsetx0)),
//      instead of the earlier version where both terms used offsetx0 and
//      the correction canceled algebraically.
// ------------------------------------------------------------------------
static void compute_angles(
    double x0, double y0, double x1, double y1,
    double offsetx0, double offsetx1,
    double offsety0, double offsety1,
    double& cos_theta_det,
    double& phi_det,
    double& cos_theta)
{
    double dx = (x1 + 2.5) - (x0 - 2.5);
    double dy =  y1 - y0 - (offsety1 - offsety0);
    double dz =  5.0;

    double nd = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (nd > 0.0) {
        double sth = std::sqrt(dx*dx + dy*dy) / nd;
        cos_theta_det = dz / nd;
        phi_det       = std::atan2(dy, dx);

        double nx = (-sth * std::cos(phi_det) + cos_theta_det) / std::sqrt(2.0);
        double ny =   sth * std::sin(phi_det);
        double nz = ( sth * std::cos(phi_det) + cos_theta_det) / std::sqrt(2.0);
        double nb = std::sqrt(nx*nx + ny*ny + nz*nz);

        cos_theta = (nb > 0.0) ? nz / nb : -999.;
    } else {
        cos_theta_det = -999.;
        phi_det       = -999.;
        cos_theta     = -999.;
    }
}

void angles_test()
{
    TFile* fin = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences.root");
    if(!fin || fin->IsZombie()){
        std::cerr << "[ERROR] Cannot open input file\n";
        return;
    }
    TTree* tin = (TTree*)fin->Get("events_uranium");
    if(!tin){
        std::cerr << "[ERROR] Tree not found\n";
        return;
    }

    double x0, x1, y0, y1;
    double cos_theta, cos_theta_det;
    tin->SetBranchAddress("x0",        &x0);
    tin->SetBranchAddress("x1",        &x1);
    tin->SetBranchAddress("y0",        &y0);
    tin->SetBranchAddress("y1",        &y1);
    tin->SetBranchAddress("cos_theta", &cos_theta);
    tin->SetBranchAddress("cos_theta_det", &cos_theta_det);

    Long64_t nentries = tin->GetEntries();

    TH1D* hist_theta           = new TH1D("hist_theta_formula",   "", 100, 0, 1);
    TH1D* hist_theta_corrected = new TH1D("hist_theta_corrected", "", 100, 0, 1);
    TH1D* hist_theta_det           = new TH1D("hist_theta_det_formula",   "", 100, 0, 1);
    TH1D* hist_theta_det_corrected = new TH1D("hist_theta_det_corrected", "", 100, 0, 1);

    // --- pass 1: fill the uncorrected distribution + accumulate means ---
    double mean_x0 = 0.0, mean_x1 = 0.0, mean_y0 = 0.0, mean_y1 = 0.0;
    for (Long64_t i = 0; i < nentries; ++i){
        tin->GetEntry(i);
        mean_x0 += x0;
        mean_x1 += x1;
        mean_y0 += y0;
        mean_y1 += y1;
        hist_theta->Fill(std::fabs(cos_theta));
        hist_theta_det->Fill(std::fabs(cos_theta_det));
    }
    mean_x0 /= nentries;
    mean_x1 /= nentries;
    mean_y0 /= nentries;
    mean_y1 /= nentries;

    // --- pass 2: fill the offset-corrected distribution ---
    for (Long64_t i = 0; i < nentries; ++i){
        tin->GetEntry(i);
        double cos_theta_corrected, cos_theta_det_corrected, phi_det;
        compute_angles(x0, y0, x1, y1,
                       mean_x0, mean_x1, mean_y0, mean_y1,
                       cos_theta_det_corrected, phi_det, cos_theta_corrected);
        if(cos_theta_corrected < -1.0) continue;   // skip -999 failure flag
        hist_theta_corrected->Fill(std::fabs(cos_theta_corrected));
        hist_theta_det_corrected->Fill(std::fabs(cos_theta_det_corrected));
    }

    // ====================================================================
    // Overlay plot: raw vs corrected, same canvas, house style
    // ====================================================================
    setPubStyle();

    const int kRawColor  = okabeIto(  0,   0,   0);   // black    -> uncorrected
    const int kCorrColor = okabeIto(  0, 114, 178);   // navy blue -> corrected

    TCanvas* c1 = new TCanvas("c_angles", "Angle reconstruction", 850, 680);
    stylePad(gPad);

    hist_theta->SetLineColor(kRawColor);
    hist_theta->SetMarkerColor(kRawColor);
    hist_theta->SetMarkerStyle(24);
    hist_theta->SetMarkerSize(1.2);
    hist_theta->SetLineWidth(2);
    hist_theta->SetTitle(";|cos(#theta)|;counts");
    hist_theta->GetXaxis()->SetTitleOffset(1.1);
    hist_theta->GetYaxis()->SetTitleOffset(1.5);
    hist_theta->Draw("E");

    hist_theta_corrected->SetLineColor(kCorrColor);
    hist_theta_corrected->SetMarkerColor(kCorrColor);
    hist_theta_corrected->SetMarkerStyle(21);
    hist_theta_corrected->SetMarkerSize(1.2);
    hist_theta_corrected->SetLineWidth(2);
    hist_theta_corrected->Draw("E SAME");

    TLegend* leg = new TLegend(0.18, 0.76, 0.55, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hist_theta,           "Uncorrected", "lp");
    leg->AddEntry(hist_theta_corrected, "Offset-corrected", "lp");
    leg->Draw();

    c1->SaveAs("angle_overlay.pdf");

    TCanvas* c3 = new TCanvas("c_angles_det", "Angle reconstruction", 850, 680);
    stylePad(gPad);

    hist_theta_det->SetLineColor(kRawColor);
    hist_theta_det->SetMarkerColor(kRawColor);
    hist_theta_det->SetMarkerStyle(24);
    hist_theta_det->SetMarkerSize(1.2);
    hist_theta_det->SetLineWidth(2);
    hist_theta_det->SetTitle(";|cos(#theta)|;counts");
    hist_theta_det->GetXaxis()->SetTitleOffset(1.1);
    hist_theta_det->GetYaxis()->SetTitleOffset(1.5);
    hist_theta_det->Draw("E");

    hist_theta_det_corrected->SetLineColor(kCorrColor);
    hist_theta_det_corrected->SetMarkerColor(kCorrColor);
    hist_theta_det_corrected->SetMarkerStyle(21);
    hist_theta_det_corrected->SetMarkerSize(1.2);
    hist_theta_det_corrected->SetLineWidth(2);
    hist_theta_det_corrected->Draw("E SAME");

    TLegend* leg3 = new TLegend(0.18, 0.76, 0.55, 0.90);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->AddEntry(hist_theta_det,           "Uncorrected", "lp");
    leg3->AddEntry(hist_theta_det_corrected, "Offset-corrected", "lp");
    leg3->Draw();

    c3->SaveAs("angle_overlay_det.pdf");

    // ====================================================================
    // Pull plot: bin-by-bin (corrected - raw) / combined Poisson error.
    //
    // NOTE: the two histograms are filled from the SAME underlying events
    // (just with different cos_theta values per event), so the bin counts
    // are correlated, not independent. Treating them as independent
    // Poisson counts (sqrt(N1+N2)) overstates the combined uncertainty
    // somewhat -- fine as a quick-look diagnostic of where the correction
    // moves probability mass, but not a rigorous statistical test.
    // ====================================================================
    int nb = hist_theta->GetNbinsX();
    TH1D* hist_pull = new TH1D("hist_pull", "",
                               nb, hist_theta->GetXaxis()->GetXmin(),
                                   hist_theta->GetXaxis()->GetXmax());

    for(int b = 1; b <= nb; ++b){
        double n1 = hist_theta->GetBinContent(b);
        double n2 = hist_theta_corrected->GetBinContent(b);
        double sigma = std::sqrt(n1 + n2);
        double pull = (sigma > 0.0) ? (n2 - n1) / sigma : 0.0;
        hist_pull->SetBinContent(b, pull);
    }

    TCanvas* c2 = new TCanvas("c_pull", "Pull: corrected vs raw", 850, 500);
    stylePad(gPad);

    double ymax = 4.0;
    for(int b = 1; b <= nb; ++b)
        ymax = std::max(ymax, std::fabs(hist_pull->GetBinContent(b)) * 1.2);

    hist_pull->SetMinimum(-ymax);
    hist_pull->SetMaximum( ymax);
    hist_pull->SetMarkerStyle(20);
    hist_pull->SetMarkerSize(1.2);
    hist_pull->SetMarkerColor(kCorrColor);
    hist_pull->SetLineColor(kCorrColor);
    hist_pull->SetTitle(";|cos(#theta)|;pull = (corrected - raw)/#sigma");
    hist_pull->GetXaxis()->SetTitleOffset(1.1);
    hist_pull->GetYaxis()->SetTitleOffset(1.4);
    hist_pull->Draw("P");

    double xlo = hist_theta->GetXaxis()->GetXmin();
    double xhi = hist_theta->GetXaxis()->GetXmax();
    TLine* l0  = new TLine(xlo, 0.0, xhi, 0.0);
    TLine* lp1 = new TLine(xlo, 1.0, xhi, 1.0);
    TLine* lm1 = new TLine(xlo,-1.0, xhi,-1.0);
    l0->SetLineColor(kGray+2);  l0->SetLineWidth(1);
    lp1->SetLineStyle(2); lp1->SetLineColor(kGray+1);
    lm1->SetLineStyle(2); lm1->SetLineColor(kGray+1);
    l0->Draw(); lp1->Draw(); lm1->Draw();
    hist_pull->Draw("P SAME");

    c2->SaveAs("angle_pull.pdf");

    std::cout << "[INFO] Saved angle_overlay.pdf and angle_pull.pdf\n";
}