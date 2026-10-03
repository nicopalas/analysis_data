// =============================================================================
//  gold_uranium_ratio.C
//
//  sigma_Au(n,f) / sigma_U(n,f) from the outputs of fission_xs.C
//  (h_cs_raw_gold, h_cs_raw_uranium in cs_toy_gold.root / cs_toy_uranium.root).
//
// -----------------------------------------------------------------------------
//  CHANGE WITH RESPECT TO THE PREVIOUS VERSION: NORMALISATION METHOD
//
//  OLD: the raw ratio was scaled so that ONE bin (the one containing
//  E = 73.9 MeV) matched the reference ratio at that energy exactly. All
//  the normalisation weight sat on a single point: a fluctuation, a bad
//  bin edge, or a typo in one reference value would silently propagate
//  into the scale of every other point.
//
//  NEW: same philosophy as computeScale() in fission_xs.C. For every gold
//  reference energy E_i that (a) has a valid interpolated uranium
//  reference value and (b) falls in a non-empty bin of the measured raw
//  ratio, compute
//
//      R_i  = ratio_ref(E_i) / ratio_raw(E_i)
//      u_R_i = R_i * sqrt( (u_ratio_ref/ratio_ref)^2 + (u_ratio_raw/ratio_raw)^2 )
//
//  and take the inverse-variance weighted mean:
//
//      scale   = ( sum_i R_i / u_R_i^2 ) / ( sum_i 1 / u_R_i^2 )
//      u_scale = 1 / sqrt( sum_i 1 / u_R_i^2 )              (correlated,
//                                                             applied to
//                                                             every bin)
//      chi2    = sum_i (R_i - scale)^2 / u_R_i^2,  ndf = n_used - 1
//
//  chi2/ndf tests whether the RAW ratio and the REFERENCE ratio agree in
//  shape over the whole overlap, not just at one point. A bad chi2/ndf is
//  a real warning that normalising by a single scale factor is not
//  appropriate (there is an energy-dependent discrepancy left over).
//
//  All points actually used in the fit are marked with an open ring in
//  the upper panel, instead of a single "normalisation point" marker.
//  Points outside the measured range (e.g. the 46.3 MeV gold reference
//  point, below the 60 MeV threshold of fission_xs.C) fall in an empty
//  bin and are automatically excluded, exactly like in computeScale().
//
//  If NO reference energy overlaps the measured range, no normalisation
//  is applied and the RAW ratio is plotted, with a clear warning -- same
//  fallback behaviour as fission_xs.C for uranium's own cross section.
//
// -----------------------------------------------------------------------------
//  OTHER SMALL FIXES ALONG THE WAY
//   - Bin-compatibility check now compares actual edges, not just the
//     number of bins.
//   - Graph x-positions use the geometric bin centre (sqrt(lo*hi)),
//     matching fission_xs.C's convention for log-spaced bins, instead of
//     ROOT's arithmetic GetBinCenter().
//   - dynamic_cast instead of C-style casts on TFile::Get() results.
//   - TFile* fout is deleted after Close().
//
//  STILL DUPLICATED, NOT YET SHARED (see gold_xs()/uranium_xs()):
//   - The gold and uranium reference tables are copied here rather than
//     read from a common header. If you correct a reference value in one
//     place, remember to update it here too.
// =============================================================================

#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/cross_section.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGaxis.h"

#include <vector>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>

// ============================================================================
//  House style — same palette and conventions as plots.h (anisotropy
//  analysis), kept local so this macro does not depend on that project.
// ============================================================================
static int rsColor(const char* hex) { return TColor::GetColor(hex); }

static const int kRsInk      = rsColor("#1B1B1B");  // axes, text
static const int kRsRefGray  = rsColor("#8C8C8C");  // reference lines
static const int kRsBandGray = rsColor("#E6E6E6");  // uncertainty bands
static const int kRsThisWork = rsColor("#1F4E79");  // deep blue   -> this work
static const int kRsRefData  = rsColor("#D1603D");  // terracotta  -> reference

static void rsSetStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameLineColor(kRsInk);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetEndErrorSize(0);
    gStyle->SetErrorX(0.);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
    TGaxis::SetMaxDigits(4);
}

// axis styling; s rescales sizes for pads shorter than the canvas
static void rsFrame(TH1* fr, const char* xt, const char* yt, double s, double yoff)
{
    TAxis* ax = fr->GetXaxis();
    TAxis* ay = fr->GetYaxis();
    ax->SetTitle(xt);
    ay->SetTitle(yt);
    for (TAxis* a : {ax, ay}) {
        a->SetLabelFont(42);
        a->SetTitleFont(42);
        a->SetLabelSize(0.042 * s);
        a->SetTitleSize(0.048 * s);
        a->SetAxisColor(kRsInk);
        a->SetLabelColor(kRsInk);
        a->SetTitleColor(kRsInk);
    }
    ax->SetLabelOffset(0.008 * s);
    ay->SetLabelOffset(0.010);
    ax->SetTitleOffset(1.05);
    ay->SetTitleOffset(yoff);
    ax->SetTickLength(0.030 * s);
    ay->SetTickLength(0.018);
    ax->SetMoreLogLabels(kFALSE);   // decades only: 10, 100, 1000
    ax->SetNoExponent(kTRUE);
}

static TLine* rsHLine(double x1, double x2, double y, int style = 2, int color = kRsRefGray)
{
    TLine* l = new TLine(x1, y, x2, y);
    l->SetLineStyle(style);
    l->SetLineColor(color);
    l->SetLineWidth(1);
    l->Draw();
    return l;
}

static TGraph* rsBand(double x1, double x2, double ylo, double yhi, int color = kRsBandGray)
{
    TGraph* b = new TGraph(4);
    b->SetPoint(0, x1, ylo);
    b->SetPoint(1, x2, ylo);
    b->SetPoint(2, x2, yhi);
    b->SetPoint(3, x1, yhi);
    b->SetFillColor(color);
    b->SetFillStyle(1001);
    b->SetLineColor(color);
    b->SetLineWidth(0);
    b->Draw("F SAME");
    return b;
}

static TLegend* rsLegend(double x1, double y1, double x2, double y2, double tsize)
{
    TLegend* lg = new TLegend(x1, y1, x2, y2);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->SetTextFont(42);
    lg->SetTextSize(tsize);
    lg->SetTextColor(kRsInk);
    lg->SetMargin(0.16);
    return lg;
}

// geometric bin centre, consistent with fission_xs.C's log-bin convention
static double geomCenter(double lo, double hi) { return std::sqrt(lo * hi); }

// ============================================================================
//  Result of the weighted-mean normalisation of the ratio
// ============================================================================
struct RatioNormalisation {
    bool                ok      = false;
    double              scale   = 1.0;
    double              u_scale = 0.0;
    double              chi2    = 0.0;
    int                 n_used  = 0;
    std::vector<double> E_used;    // energies actually used in the fit
    std::vector<double> R_used;    // per-point ratio ref/raw
    std::vector<double> uR_used;   // their uncertainties
};

// Weighted mean of R_i = ratio_ref(E_i) / ratio_raw(E_i) over every gold
// reference energy that (a) has a valid interpolated uranium reference and
// (b) lands in a non-empty bin of h_ratio_raw. Mirrors computeScale() in
// fission_xs.C, applied to the ratio instead of to a single cross section.
static RatioNormalisation computeRatioScale(TH1D* h_ratio_raw,
                                            const TGraphErrors* gr_ratio_ref)
{
    RatioNormalisation res;
    double sw = 0.0, swR = 0.0;

    for (int k = 0; k < gr_ratio_ref->GetN(); ++k) {
        const double E    = gr_ratio_ref->GetX()[k];
        const double rref  = gr_ratio_ref->GetY()[k];
        const double urref = gr_ratio_ref->GetEY()[k];
        if (rref <= 0.0) continue;

        const int    b   = h_ratio_raw->FindBin(E);
        const double x   = h_ratio_raw->GetBinContent(b);
        const double ux  = h_ratio_raw->GetBinError(b);
        if (x <= 0.0 || ux <= 0.0) continue;   // outside measured range / empty bin

        const double R  = rref / x;
        const double uR = R * std::sqrt(std::pow(urref / rref, 2) + std::pow(ux / x, 2));
        if (uR <= 0.0) continue;
        const double w = 1.0 / (uR * uR);

        sw  += w;
        swR += w * R;
        res.E_used.push_back(E);
        res.R_used.push_back(R);
        res.uR_used.push_back(uR);
        printf("[norm] E = %6.1f MeV  ratio_raw = %.4e  ratio_ref = %.4e  "
               "R = %.4e +/- %.4e\n", E, x, rref, R, uR);
    }
    if (res.E_used.empty()) return res;

    res.ok     = true;
    res.n_used = (int)res.E_used.size();
    res.scale  = swR / sw;
    res.u_scale = 1.0 / std::sqrt(sw);
    for (size_t i = 0; i < res.R_used.size(); ++i)
        res.chi2 += std::pow(res.R_used[i] - res.scale, 2) / std::pow(res.uR_used[i], 2);
    return res;
}

// ============================================================================
void gold_uranium_ratio(){

    const TString kOutDir = "/Users/nico/Desktop/Tese/Analysis/cross_section/output/";

    // ── open the two histograms produced by gold_xs() and uranium_xs() ────────
    TFile *f_gold = TFile::Open(kOutDir + "cs_toy_gold.root", "READ");
    if (!f_gold || f_gold->IsZombie()) { std::cerr << "Cannot open gold output file\n"; return; }
    TH1D *h_gold = dynamic_cast<TH1D*>(f_gold->Get("h_cs_raw_gold"));
    if (!h_gold) { std::cerr << "h_cs_raw_gold not found\n"; return; }
    h_gold = (TH1D*)h_gold->Clone("h_cs_raw_gold_clone");
    h_gold->SetDirectory(nullptr);
    f_gold->Close();

    TFile *f_u = TFile::Open(kOutDir + "cs_toy_uranium.root", "READ");
    if (!f_u || f_u->IsZombie()) { std::cerr << "Cannot open uranium output file\n"; return; }
    TH1D *h_u = dynamic_cast<TH1D*>(f_u->Get("h_cs_raw_uranium"));
    if (!h_u) { std::cerr << "h_cs_raw_uranium not found\n"; return; }
    h_u = (TH1D*)h_u->Clone("h_cs_raw_u_clone");
    h_u->SetDirectory(nullptr);
    f_u->Close();

    // ── check the two histograms really share the same binning ────────────────
    if (h_gold->GetNbinsX() != h_u->GetNbinsX()) {
        std::cerr << "Gold and uranium histograms have a different number of bins\n"; return;
    }
    const int nbins = h_gold->GetNbinsX();
    bool same_edges = true;
    for (int e = 1; e <= nbins + 1 && same_edges; ++e)
        same_edges = std::fabs(h_gold->GetXaxis()->GetBinLowEdge(e) -
                               h_u->GetXaxis()->GetBinLowEdge(e)) < 1e-9;
    if (!same_edges) { std::cerr << "Gold and uranium histograms have different bin edges\n"; return; }

    // ── raw bin-by-bin ratio Au/U (unscaled) ───────────────────────────────────
    // Note: the flux cancels EXACTLY here, since both sigma_raw histograms
    // were built from the same flux integration over the same binning
    // (see fission_xs.C). Acceptance and efficiency do NOT cancel: they
    // enter each numerator/denominator inside two different weighted sums
    // over two different count distributions, not as a common factor.
    TH1D *h_ratio_raw = (TH1D*)h_gold->Clone("h_ratio_raw");
    h_ratio_raw->SetDirectory(nullptr);
    h_ratio_raw->Reset();

    for (int e = 1; e <= nbins; ++e) {
        double g   = h_gold->GetBinContent(e);
        double ug  = h_gold->GetBinError(e);
        double u   = h_u->GetBinContent(e);
        double uu  = h_u->GetBinError(e);
        if (g <= 0.0 || u <= 0.0) continue;

        double r   = g / u;
        double ur  = r * std::sqrt((ug*ug)/(g*g) + (uu*uu)/(u*u));

        h_ratio_raw->SetBinContent(e, r);
        h_ratio_raw->SetBinError  (e, ur);
    }

    // ── reference gold data (Au-197), same as in gold_xs() ─────────────────────
    //    TODO: move into a shared reference_data.h so this table and the one
    //    in gold_xs() cannot drift apart.
    std::vector<double> E_ref_au   = {46.3, 66.6, 73.9, 94.1, 132.9, 144.6, 173.3};
    std::vector<double> sig_ref_au = {0.103, 0.81, 1.20, 2.81, 6.1,   8.1,   10.3}; // mb
    std::vector<double> u_ref_au   = {0.019, 0.12, 0.17, 0.39, 0.9,   1.2,   1.6};  // mb
    for (auto &s : sig_ref_au) s *= 1e-3;   // mb -> barn
    for (auto &u : u_ref_au)   u *= 1e-3;

    // ── reference uranium data (U-238(n,f)), full ENDF-style table ─────────────
    //    TODO: same as above, share with uranium_xs() once it uses the full
    //    table too (it currently only sees the 42-50 MeV subset, which is
    //    why it cannot normalise itself - see fission_xs.C).
    std::vector<double> E_ref_u = {
        0.50,0.52,0.54,0.57,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.94,0.96,0.98,
        1.00,1.10,1.25,1.40,1.60,1.80,2.00,2.20,2.40,2.60,2.80,3.00,3.60,4.00,
        4.50,4.70,5.00,5.30,5.50,5.80,6.00,6.20,6.50,7.00,7.50,7.75,8.00,8.50,
        9.00,10.00,11.00,11.50,12.00,13.00,14.00,14.50,15.00,16.00,17.00,18.00,
        19.00,20.00,21.00,22.00,23.00,24.00,25.00,26.00,27.00,28.00,29.00,30.00,
        32.00,34.00,36.00,38.00,40.00,42.00,44.00,46.00,48.00,50.00,52.00,54.00,
        56.00,58.00,60.00,64.00,68.00,72.00,76.00,80.00,84.00,88.00,92.00,96.00,
        100.00,104.00,108.00,112.00,116.00,120.00,128.00,136.00,144.00,152.00,
        160.00,168.00,176.00,184.00,192.00,200.00,300.00,400.00,500.00,600.00,
        700.00,800.00,900.00,1000.00
    };
    std::vector<double> sig_ref_u = {
        0.00026980,0.00067119,0.00059598,0.00065116,0.00116531,0.00129202,0.00184360,
        0.00264695,0.00450784,0.00675016,0.01378015,0.01690077,0.01550547,0.01588961,
        0.01410576,0.02895309,0.03320903,0.18619023,0.41835630,0.48206710,0.53553997,
        0.54749355,0.54624349,0.54151562,0.53708920,0.52481337,0.54752327,0.55451845,
        0.55921914,0.55980956,0.54833662,0.55240729,0.54820060,0.56779697,0.61392195,
        0.68694589,0.82159525,0.95004294,0.99851630,0.99569271,1.01920402,1.01472289,
        1.01452311,1.01117324,1.00882950,1.00505332,0.98752794,1.03156514,1.14994058,
        1.19263970,1.24300610,1.32225423,1.32703520,1.32162220,1.35796397,1.40707482,
        1.51901925,1.55512900,1.60213378,1.54300822,1.57505132,1.58626925,1.56289817,
        1.62440482,1.61076480,1.65532032,1.70277087,1.69981759,1.65676440,1.63970023,
        1.66526066,1.65841277,1.65557833,1.68089904,1.64272747,1.61694393,1.64917876,
        1.61998653,1.63154717,1.62447590,1.58763664,1.57865896,1.52995539,1.50445837,
        1.51984211,1.49515595,1.49600313,1.44827525,1.40759496,1.39989479,1.42732347,
        1.38255921,1.42826710,1.35680338,1.36360367,1.33475089,1.32909427,1.29432732,
        1.28907594,1.32032323,1.31550749,1.29939456,1.33658872,1.34154384,1.31700244,
        1.32466100,1.44100000,1.49700000,1.46100000,1.45800000,1.47500000,1.48200000,
        1.46800000,1.46700001
    };
    std::vector<double> u_ref_u = {
        0.00002049,0.00010438,0.00002249,0.00002441,0.00003694,0.00002922,0.00004329,
        0.00004879,0.00008038,0.00010401,0.00017266,0.00021365,0.00023027,0.00021990,
        0.00015131,0.00029119,0.00032341,0.00134738,0.00262918,0.00302871,0.00290432,
        0.00301741,0.00311783,0.00300113,0.00344929,0.00290484,0.00298915,0.00307938,
        0.00360237,0.00358805,0.00362927,0.00400621,0.00365471,0.00415490,0.00463600,
        0.00555936,0.00612914,0.00688432,0.00719061,0.00859903,0.00773614,0.00729852,
        0.00729594,0.00815660,0.00806416,0.00962545,0.00864000,0.00748159,0.00674434,
        0.00570974,0.00834980,0.01067722,0.01146856,0.01263919,0.01082061,0.01487632,
        0.01428978,0.02590652,0.01836317,0.02226853,0.01699061,0.02641331,0.02127290,
        0.02271004,0.02358701,0.01967756,0.02992725,0.02426689,0.02619886,0.03118292,
        0.02826156,0.02867969,0.02690437,0.03005004,0.02926416,0.03718286,0.03321005,
        0.03569087,0.04041503,0.03035458,0.03109424,0.02486835,0.02395078,0.02856810,
        0.02707488,0.03238298,0.03572723,0.03841870,0.03488576,0.03670170,0.03425680,
        0.03842244,0.04389951,0.04586037,0.04611334,0.04018981,0.03825088,0.03790471,
        0.03312034,0.03876146,0.03685174,0.04049074,0.03526862,0.05257003,0.03563506,
        0.03218347,0.07185272,0.10454333,0.08123140,0.08067812,0.08530600,0.06620452,
        0.06558076,0.06740687
    };
    // sig_ref_u, u_ref_u already in barn

    std::vector<double> ex_ref_au(E_ref_au.size(), 0.0);
    std::vector<double> ex_ref_u(E_ref_u.size(), 0.0);

    TGraphErrors *gr_ref_au = new TGraphErrors(
        (int)E_ref_au.size(), E_ref_au.data(), sig_ref_au.data(),
        ex_ref_au.data(), u_ref_au.data());
    gr_ref_au->SetName("g_ref_au");

    TGraphErrors *gr_ref_u = new TGraphErrors(
        (int)E_ref_u.size(), E_ref_u.data(), sig_ref_u.data(),
        ex_ref_u.data(), u_ref_u.data());
    gr_ref_u->SetName("g_ref_u");

    // uranium uncertainty interpolated the same way as its value
    TGraph gr_ref_u_err((int)E_ref_u.size(), E_ref_u.data(), u_ref_u.data());
    const double u_ref_u_Emin = E_ref_u.front(), u_ref_u_Emax = E_ref_u.back();

    // ── reference ratio Au/U at the gold reference energies ────────────────────
    //    All 7 gold energies lie inside [0.5, 1000] MeV, so the uranium
    //    interpolation never extrapolates here; the range check below is a
    //    safeguard if this table is edited later.
    TGraphErrors *gr_ratio_ref = new TGraphErrors();
    gr_ratio_ref->SetName("g_ratio_ref");
    gr_ratio_ref->SetTitle("Reference ratio;E_{n} (MeV);Ratio");
    for (size_t k = 0; k < E_ref_au.size(); ++k) {
        const double E = E_ref_au[k];
        if (E < u_ref_u_Emin || E > u_ref_u_Emax) {
            std::cerr << "[ref] " << E << " MeV outside the uranium reference "
                      << "range, skipped for the ratio\n";
            continue;
        }
        const double u_val = gr_ref_u->Eval(E);
        const double u_err = gr_ref_u_err.Eval(E);
        const double a_val = sig_ref_au[k];
        const double a_err = u_ref_au[k];
        if (u_val <= 0.0 || a_val <= 0.0) continue;
        const double r  = a_val / u_val;
        const double ur = r * std::sqrt(std::pow(a_err/a_val, 2) + std::pow(u_err/u_val, 2));
        const int ip = gr_ratio_ref->GetN();
        gr_ratio_ref->SetPoint(ip, E, r);
        gr_ratio_ref->SetPointError(ip, 0.0, ur);
    }

    // ── weighted-mean normalisation over every overlapping point ──────────────
    const RatioNormalisation norm = computeRatioScale(h_ratio_raw, gr_ratio_ref);

    TH1D *h_ratio_norm = (TH1D*)h_ratio_raw->Clone("h_ratio_norm");
    h_ratio_norm->SetDirectory(nullptr);

    if (norm.ok) {
        h_ratio_norm->Scale(norm.scale);   // scales contents and errors
        const int ndf = norm.n_used - 1;
        printf("[norm] scale = %.4e +/- %.4e (%.1f%%, correlated)  from %d point(s)",
               norm.scale, norm.u_scale, 100.0 * norm.u_scale / norm.scale, norm.n_used);
        if (ndf > 0) printf("  chi2/ndf = %.2f/%d", norm.chi2, ndf);
        printf("\n");
    } else {
        std::cerr << "[norm] WARNING: no gold reference energy falls in a "
                     "non-empty bin of the measured ratio. NO normalisation "
                     "applied; the RAW ratio is plotted.\n";
    }
    const double scale_for_plot   = norm.ok ? norm.scale   : 1.0;
    const double u_scale_rel_plot = norm.ok ? (norm.u_scale / norm.scale) : 0.0;

    // ── graphs of this work: file version (symmetric) + drawing version ──────
    //    x position and x errors use the geometric bin centre, matching
    //    fission_xs.C's convention for the underlying log-spaced bins.
    std::vector<double> x_r, y_r, ex_r, ey_r, exl_r, exh_r;
    for (int e = 1; e <= nbins; ++e) {
        if (h_ratio_norm->GetBinContent(e) <= 0.0) continue;
        const double lo = h_ratio_norm->GetXaxis()->GetBinLowEdge(e);
        const double hi = h_ratio_norm->GetXaxis()->GetBinUpEdge(e);
        const double xc = geomCenter(lo, hi);
        x_r.push_back(xc);
        y_r.push_back(h_ratio_norm->GetBinContent(e));
        ex_r.push_back(0.5 * (hi - lo));
        exl_r.push_back(xc - lo);
        exh_r.push_back(hi - xc);
        ey_r.push_back(h_ratio_norm->GetBinError(e));
    }
    if (x_r.empty()) { std::cerr << "Empty ratio, nothing to plot\n"; return; }

    TGraphErrors *gr_ratio = new TGraphErrors(
        (int)x_r.size(), x_r.data(), y_r.data(), ex_r.data(), ey_r.data());
    gr_ratio->SetName("g_ratio_au_u");
    gr_ratio->SetTitle(norm.ok ? "#sigma_{Au}/#sigma_{U} (normalised);E_{n} (MeV);Ratio"
                                : "#sigma_{Au}/#sigma_{U} (raw);E_{n} (MeV);Ratio");

    TGraphAsymmErrors *gd_ratio = new TGraphAsymmErrors(
        (int)x_r.size(), x_r.data(), y_r.data(), exl_r.data(), exh_r.data(),
        ey_r.data(), ey_r.data());

    // points that actually entered the weighted-mean fit, marked with a ring
    TGraph *g_used = new TGraph();
    for (int k = 0; k < norm.n_used; ++k) {
        const int b = h_ratio_norm->FindBin(norm.E_used[k]);
        g_used->SetPoint(k, geomCenter(h_ratio_norm->GetXaxis()->GetBinLowEdge(b),
                                       h_ratio_norm->GetXaxis()->GetBinUpEdge(b)),
                         h_ratio_norm->GetBinContent(b));
    }

    // ── lower panel: this work / reference at the reference energies ─────────
    //    Residuals of the same fit performed above: points used in the
    //    weighted mean should scatter around 1 within their own errors,
    //    consistent with the reported chi2/ndf.
    TGraphErrors *gr_cmp = new TGraphErrors();
    gr_cmp->SetName("g_ratio_this_over_ref");
    for (int k = 0; k < gr_ratio_ref->GetN(); ++k) {
        const double E    = gr_ratio_ref->GetX()[k];
        const double yref = gr_ratio_ref->GetY()[k];
        const double eref = gr_ratio_ref->GetEY()[k];
        const int    b    = h_ratio_norm->FindBin(E);
        const double y    = h_ratio_norm->GetBinContent(b);
        const double ey   = h_ratio_norm->GetBinError(b);
        if (y <= 0.0 || yref <= 0.0) continue;
        const double q  = y / yref;
        const double eq = q * std::sqrt(std::pow(ey/y, 2) + std::pow(eref/yref, 2));
        const int ip = gr_cmp->GetN();
        gr_cmp->SetPoint(ip, E, q);
        gr_cmp->SetPointError(ip, 0.0, eq);
    }

    // ── ranges ────────────────────────────────────────────────────────────────
    double xmin = 1e30, xmax = 0.0, ymin = 1e30, ymax = 0.0;
    for (size_t i = 0; i < x_r.size(); ++i) {
        xmin = std::min(xmin, x_r[i] - exl_r[i]);
        xmax = std::max(xmax, x_r[i] + exh_r[i]);
        ymin = std::min(ymin, (y_r[i] - ey_r[i] > 0.0) ? y_r[i] - ey_r[i] : y_r[i]);
        ymax = std::max(ymax, y_r[i] + ey_r[i]);
    }
    for (int k = 0; k < gr_ratio_ref->GetN(); ++k) {
        const double v = gr_ratio_ref->GetY()[k], ev = gr_ratio_ref->GetEY()[k];
        xmin = std::min(xmin, gr_ratio_ref->GetX()[k]);
        xmax = std::max(xmax, gr_ratio_ref->GetX()[k]);
        ymin = std::min(ymin, (v - ev > 0.0) ? v - ev : v);
        ymax = std::max(ymax, v + ev);
    }
    xmin = std::max(xmin * 0.9, 1e-3);
    xmax = xmax * 1.1;
    const double ylo  = ymin / 2.5;
    // leave the top ~30% of the upper frame (in log) for the legend
    const double yTop = ylo * std::pow(ymax * 1.5 / ylo, 1.0 / (1.0 - 0.30));

    double dq = 1.5 * u_scale_rel_plot;
    for (int k = 0; k < gr_cmp->GetN(); ++k)
        dq = std::max(dq, std::fabs(gr_cmp->GetY()[k] - 1.0) + gr_cmp->GetEY()[k]);
    dq = std::max(0.2, 1.25 * dq);

    // ── drawing ───────────────────────────────────────────────────────────────
    rsSetStyle();

    const double fUp = 0.70, fLo = 0.30;
    const double sUp = 1.0 / fUp, sLo = 1.0 / fLo;
    const double kLeft = 0.14, kRight = 0.04;

    TCanvas *c = new TCanvas("c_ratio", "Au/U cross section ratio", 900, 820);
    TPad *pUp = new TPad("pad_ratio_up", "", 0.0, fLo, 1.0, 1.0);
    TPad *pLo = new TPad("pad_ratio_lo", "", 0.0, 0.0, 1.0, fLo);
    pUp->Draw();
    pLo->Draw();

    // ===== upper panel: sigma_Au / sigma_U ====================================
    pUp->cd();
    pUp->SetLeftMargin(kLeft);
    pUp->SetRightMargin(kRight);
    pUp->SetTopMargin(0.045 * sUp);
    pUp->SetBottomMargin(0.015);
    pUp->SetTickx(1); pUp->SetTicky(1);
    pUp->SetLogx();
    pUp->SetLogy();

    TH1F *fr = pUp->DrawFrame(xmin, ylo, xmax, yTop);
    rsFrame(fr, "E_{n} (MeV)", "#sigma_{Au} / #sigma_{U}", sUp, 1.25);
    fr->GetXaxis()->SetLabelSize(0.0);
    fr->GetXaxis()->SetTitleSize(0.0);

    // reference: terracotta open squares
    gr_ratio_ref->SetMarkerStyle(21);
    gr_ratio_ref->SetMarkerSize(1.3);
    gr_ratio_ref->SetMarkerColor(kRsRefData);
    gr_ratio_ref->SetLineColor(kRsRefData);
    gr_ratio_ref->SetLineWidth(1);
    gr_ratio_ref->Draw("P");

    // this work: deep-blue filled circles, horizontal bar = energy bin
    gd_ratio->SetMarkerStyle(20);
    gd_ratio->SetMarkerSize(1.0);
    gd_ratio->SetMarkerColor(kRsThisWork);
    gd_ratio->SetLineColor(kRsThisWork);
    gd_ratio->SetLineWidth(1);
    gd_ratio->Draw("P");

    // points that entered the weighted-mean fit: open rings
    g_used->SetMarkerStyle(24);
    g_used->SetMarkerSize(1.8);
    g_used->SetMarkerColor(kRsInk);
    if (g_used->GetN() > 0) g_used->Draw("P");

    // same style on the file graph, so it looks right if redrawn from the file
    gr_ratio->SetMarkerStyle(20);
    gr_ratio->SetMarkerSize(1.0);
    gr_ratio->SetMarkerColor(kRsThisWork);
    gr_ratio->SetLineColor(kRsThisWork);

    const double yL = 1.0 - pUp->GetTopMargin() - 0.02;
    const double rowH = 0.052 * sUp;
    TLegend *leg = rsLegend(kLeft + 0.03, yL - 3 * rowH, 0.62, yL, 0.040 * sUp);
    leg->AddEntry(gd_ratio,     "This work",                             "pe");
    leg->AddEntry(gr_ratio_ref, "Reference  (^{197}Au / ^{238}U)",       "pe");
    leg->AddEntry(g_used,       Form("Used in normalisation (%d pts)", norm.n_used), "p");
    leg->Draw();

    TLatex lat; lat.SetNDC(); lat.SetTextFont(42); lat.SetTextColor(kRsInk);
    lat.SetTextSize(0.036 * sUp); lat.SetTextAlign(33);
    lat.DrawLatex(1.0 - kRight - 0.03, yL,
                  "^{197}Au(n,f) / ^{238}U(n,f)");
    lat.SetTextColor(kRsRefGray);
    if (norm.ok) {
        const int ndf = norm.n_used - 1;
        if (ndf > 0)
            lat.DrawLatex(1.0 - kRight - 0.03, yL - rowH,
                          Form("scale #pm%.1f%%  (#chi^{2}/ndf = %.1f/%d)",
                               100.0 * u_scale_rel_plot, norm.chi2, ndf));
        else
            lat.DrawLatex(1.0 - kRight - 0.03, yL - rowH,
                          Form("scale #pm%.1f%% (1 point)", 100.0 * u_scale_rel_plot));
    } else {
        lat.DrawLatex(1.0 - kRight - 0.03, yL - rowH, "not normalised (no overlap)");
    }

    pUp->RedrawAxis();

    // ===== lower panel: this work / reference ================================
    pLo->cd();
    pLo->SetLeftMargin(kLeft);
    pLo->SetRightMargin(kRight);
    pLo->SetTopMargin(0.03);
    pLo->SetBottomMargin(0.12 * sLo);
    pLo->SetTickx(1); pLo->SetTicky(1);
    pLo->SetLogx();

    TH1F *fr2 = pLo->DrawFrame(xmin, 1.0 - dq, xmax, 1.0 + dq);
    rsFrame(fr2, "E_{n} (MeV)", "This work / ref.", sLo, 1.25);
    fr2->GetYaxis()->SetNdivisions(504);
    fr2->GetYaxis()->CenterTitle();

    // grey band: correlated normalisation uncertainty of this work
    rsBand(xmin, xmax, 1.0 - u_scale_rel_plot, 1.0 + u_scale_rel_plot);
    rsHLine(xmin, xmax, 1.0, 1, kRsRefGray);

    gr_cmp->SetMarkerStyle(20);
    gr_cmp->SetMarkerSize(1.0);
    gr_cmp->SetMarkerColor(kRsThisWork);
    gr_cmp->SetLineColor(kRsThisWork);
    gr_cmp->SetLineWidth(1);
    gr_cmp->Draw("P");

    pLo->RedrawAxis();

    c->cd();
    c->SaveAs(kOutDir + "cs_toy_ratio.pdf");

    // ── save ──────────────────────────────────────────────────────────────────
    TFile *fout = TFile::Open(kOutDir + "cs_toy_ratio.root", "RECREATE");
    if (!fout || fout->IsZombie()) { std::cerr << "Cannot create cs_toy_ratio.root\n"; return; }
    c->Write();
    gr_ratio->Write();
    gr_ratio_ref->Write();
    gr_cmp->Write();
    h_ratio_raw->Write();
    h_ratio_norm->Write();
    fout->Close();
    delete fout;

    std::cout << "[DONE] cs_toy_ratio.root + cs_toy_ratio.pdf\n";
}