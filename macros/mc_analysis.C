#include "../include/types.h"
#include "../include/constants.h"
#include "../include/config.h"
#include "../include/utils.h"
#include "../include/acceptance.h"
#include "../include/efficiency.h"
#include "../include/efficiency_no_overlap.h"
#include "../include/anisotropy.h"
#include "../include/plotting.h"
#include "../include/fit_anisotropy.h"

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
#include "TStyle.h"
#include "TColor.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <limits>
#include <algorithm>
#include <fstream>
#include <sstream>

// ============================================================================
//  MC analysis for n + 238U (target_i == 2, ppac0 == 8)
//            and    n + 197Au (target_i == 1, ppac0 == 7)
//
//  Each target goes through EXACTLY the same pipeline as uranium_analysis():
//    * counts in (cos theta_beam, cos theta_det) with the acceptance cuts;
//    * efficiency with computeEfficiency();
//    * W(|cos theta|) with anisotropy();
//    * Legendre fit W(cos) = A0 (1 + a2 P2 + a4 P4) with legendre_fit();
//    * W(0)/W(90) from the fit, with the full covariance.
//
//  Two beam-angle definitions are propagated through the SAME pipeline:
//    reco (lab) : cos_theta        (PPAC reconstruction only)
//    boost (CM) : cos_theta_cm_i   (reconstruction + boost correction)
//
//  Fragment-charge window (true Z of the two fragments, from the MC):
//    U  : 80 <= Z1+Z2 <= 92   (Z_CN = 92, up to 12 charges lost before scission)
//    Au : 67 <= Z1+Z2 <= 79   (Z_CN = 79, same 12-charge window)
//    and Z1, Z2 >= 4 for both.
//
//  Per-target figures (suffix _u238 / _au197): every figure of the data
//  analysis plus the MC diagnostics (boost shift, CM vs lab, LMT).
//
//  Target comparisons (suffix _au_vs_u):
//    lmt_fraction_au_vs_u_mc     <LMT> vs E_n for both targets (+ Fatyga
//                                p + 238U), Au/U ratio below
//    boost_ratio_au_vs_u_mc      A_boost / A_reco vs E_n, both targets
//    anisotropy_cm_au_vs_u_mc    boost-corrected W(0)/W(90) vs E_n, both
//  with A = W(0 deg)/W(90 deg) from the Legendre fit.
// ============================================================================

static const int kGoldColor = hexColor("#B8860B");   // dark goldenrod -> Au

// |cos| -> bin index in [0,1]
static inline int cosBin(double c, int nbins)
{
    double a = std::fabs(c);
    if (!(a >= 0. && a <= 1.)) return -1;     // also rejects NaN
    int b = (int)(a * nbins);
    if (b >= nbins) b = nbins - 1;
    return b;
}

struct McCounts {
    Vec3D n, u_n;
};

static McCounts makeCounts(int nbins_e)
{
    McCounts m;
    m.n   = Vec3D(nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    m.u_n = Vec3D(nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    return m;
}

static void poissonErrors(McCounts& m)
{
    for (size_t e = 0; e < m.n.size(); ++e)
        for (int j = 0; j < nbins_beam; ++j)
            for (int i = 0; i < nbins_det; ++i)
                m.u_n[e][j][i] = (m.n[e][j][i] > 0.) ? std::sqrt(m.n[e][j][i]) : 1.0;
}

static void saveBoth(TCanvas* c, const std::string& outname)
{
    c->SaveAs((outname + ".pdf").c_str());
    c->SaveAs((outname + ".png").c_str());
}

static bool goodFit(const LegendreResult& r)
{
    return r.valid && std::isfinite(r.anisotropy) && std::isfinite(r.u_anisotropy)
        && r.anisotropy > 0.0 && r.u_anisotropy >= 0.0;
}

// ============================================================================
//  Linear momentum transfer statistics
// ============================================================================
struct LmtBin {
    double   lo = 0.0, hi = 0.0, xc = 0.0;
    Long64_t n  = 0;
    double   mean    = std::numeric_limits<double>::quiet_NaN();
    double   sem     = std::numeric_limits<double>::quiet_NaN();
    double   rms     = std::numeric_limits<double>::quiet_NaN();
    double   q50     = std::numeric_limits<double>::quiet_NaN();
    double   fAbove1 = std::numeric_limits<double>::quiet_NaN();   // share with LMT > 1
};

static bool okLmtBin(const LmtBin& b)
{
    return b.n >= 20 && std::isfinite(b.mean) && std::isfinite(b.sem) && std::isfinite(b.rms);
}

// Step-wise band mean +- 1 sigma (event-by-event standard deviation) over
// the energy bins: upper edge left -> right, lower edge right -> left.
static TGraph* sigmaBand(const std::vector<const LmtBin*>& ok, int fillColor)
{
    const int n = (int)ok.size();
    TGraph* g = new TGraph(4 * n);
    for (int i = 0; i < n; ++i) {
        const LmtBin& b = *ok[i];
        g->SetPoint(2 * i,             b.lo, b.mean + b.rms);
        g->SetPoint(2 * i + 1,         b.hi, b.mean + b.rms);
        g->SetPoint(4 * n - 2 - 2 * i, b.hi, b.mean - b.rms);
        g->SetPoint(4 * n - 1 - 2 * i, b.lo, b.mean - b.rms);
    }
    g->SetFillColor(fillColor);
    g->SetFillStyle(1001);
    g->SetLineColor(fillColor);
    g->SetLineWidth(0);
    return g;
}

static std::vector<double> logEdges(double lo, double hi, int n)
{
    std::vector<double> e(n + 1);
    for (int i = 0; i <= n; ++i) e[i] = lo * std::pow(hi / lo, double(i) / n);
    return e;
}

static int edgeBin(const std::vector<double>& edges, double x)
{
    if (!(x >= edges.front() && x < edges.back())) return -1;
    return (int)(std::upper_bound(edges.begin(), edges.end(), x) - edges.begin()) - 1;
}

static double quantileSorted(const std::vector<float>& v, double p)
{
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    const double h = p * (v.size() - 1);
    const size_t i = (size_t)h;
    const size_t j = std::min(i + 1, v.size() - 1);
    return v[i] + (h - i) * (v[j] - v[i]);
}

static std::vector<LmtBin> lmtStats(const std::vector<float>& E,
                                    const std::vector<float>& F,
                                    const std::vector<double>& edges)
{
    const int nb = (int)edges.size() - 1;
    std::vector<std::vector<float>> per(nb);
    for (size_t i = 0; i < E.size(); ++i) {
        const int b = edgeBin(edges, E[i]);
        if (b >= 0) per[b].push_back(F[i]);
    }

    std::vector<LmtBin> out(nb);
    for (int b = 0; b < nb; ++b) {
        LmtBin& o = out[b];
        o.lo = edges[b];
        o.hi = edges[b + 1];
        o.xc = std::sqrt(o.lo * o.hi);
        o.n  = (Long64_t)per[b].size();
        if (o.n < 2) continue;

        double s = 0.0, s2 = 0.0;
        Long64_t nAbove = 0;
        for (float v : per[b]) { s += v; s2 += double(v) * v; if (v > 1.0f) ++nAbove; }
        o.mean = s / o.n;
        const double var = std::max(0.0, s2 / o.n - o.mean * o.mean) * o.n / (o.n - 1.0);
        o.rms  = std::sqrt(var);
        o.sem  = o.rms / std::sqrt((double)o.n);
        o.fAbove1 = double(nAbove) / o.n;

        std::sort(per[b].begin(), per[b].end());
        o.q50 = quantileSorted(per[b], 0.50);
    }
    return out;
}

static void writeLmtGraphs(const std::vector<LmtBin>& bins, const std::string& suffix)
{
    std::vector<double> x, exl, exh, ym, eym, ers;
    for (const auto& b : bins) {
        if (b.n < 2) continue;
        x.push_back(b.xc);
        exl.push_back(b.xc - b.lo);
        exh.push_back(b.hi - b.xc);
        ym.push_back(b.mean);
        eym.push_back(b.sem);
        ers.push_back(b.rms);
    }
    if (x.empty()) return;
    const int n = (int)x.size();

    TGraphAsymmErrors gm(n, x.data(), ym.data(), exl.data(), exh.data(), eym.data(), eym.data());
    gm.SetTitle("LMT fraction: mean #pm s.e.m.;E_{n} (MeV);LMT fraction");
    gm.Write(("lmt_mean_" + suffix).c_str());

    TGraphAsymmErrors gs(n, x.data(), ym.data(), exl.data(), exh.data(), ers.data(), ers.data());
    gs.SetTitle("LMT fraction: mean #pm 1#sigma (event spread);E_{n} (MeV);LMT fraction");
    gs.Write(("lmt_mean_sigma_" + suffix).c_str());
}

// ============================================================================
//  Per-target containers
// ============================================================================
struct TargetSpec {
    std::string key;          // file-name suffix
    std::string name;         // legend label, e.g. "n + ^{238}U"
    std::string shortName;    // for ratio axes, e.g. "U"
    std::string tag;          // figure tag
    int    target_i;
    int    ppac0;
    int    zSumMin, zSumMax;  // window on Z1 + Z2
    int    color;
    int    marker;
    std::string acceptance_file;
    bool   compareWithRefs;   // Fatyga (p + 238U) only for uranium
    std::vector<double> ebEff;                          // efficiency energy bins
    AnalysisConfig (*makeConfig)(const std::vector<double>&, const std::string&);
};

struct TargetResult {
    TargetSpec spec;
    Vec2D      acceptance;

    McCounts eff_counts, ani_reco, ani_cm;
    std::vector<std::vector<double>> sSum, sSum2, sN;
    std::vector<float> lmtE, lmtV;
    Long64_t nTarget = 0, nZcut = 0, nsel = 0, nnoCM = 0, nnoLMT = 0;

    std::vector<EfficiencyResult> eff;
    std::vector<AnisotropyResult> ani_r, ani_c;
    std::vector<LegendreResult>   fit_r, fit_c;
    std::vector<LmtBin>           lmtFine, lmtAniso;
};

// ============================================================================
//  Two-panel layout (main + ratio), shared by the ratio figures.
//  Sizes: s = 0.8 / f_pad keeps absolute text size = single-panel figures.
// ============================================================================
struct RatioLayout {
    TCanvas* c   = nullptr;
    TPad*    top = nullptr;
    TPad*    bot = nullptr;
    double   sT  = 1.0, sB = 1.0;
};

static RatioLayout makeRatioLayout(const std::string& base, bool logx)
{
    RatioLayout L;
    const double fBot = 0.36, fTop = 1.0 - fBot;
    L.sT = 0.80 / fTop;
    L.sB = 0.80 / fBot;

    L.c = new TCanvas(uniqueName(base).c_str(), "", 820, 800);
    L.c->SetFillColor(kWhite);
    L.c->SetBorderMode(0);

    L.top = new TPad(uniqueName(base + "_top").c_str(), "", 0.0, fBot, 1.0, 1.0);
    L.bot = new TPad(uniqueName(base + "_bot").c_str(), "", 0.0, 0.0,  1.0, fBot);
    for (TPad* p : {L.top, L.bot}) {
        p->SetLeftMargin(0.15);
        p->SetRightMargin(0.04);
        p->SetTickx(1);
        p->SetTicky(1);
        if (logx) p->SetLogx();
        p->SetFillColor(kWhite);
        p->SetBorderMode(0);
    }
    L.top->SetTopMargin(figureTag().empty() ? 0.05 : 0.085);
    L.top->SetBottomMargin(0.022);
    L.bot->SetTopMargin(0.035);
    L.bot->SetBottomMargin(0.34);
    L.top->Draw();
    L.bot->Draw();
    return L;
}

// ============================================================================
//  MC-only diagnostic: <|cos theta| - |cos theta*|> vs |cos theta*|
// ============================================================================
static void plotBoostShift(
    const std::vector<double>& eEdges,
    const std::vector<std::vector<double>>& sSum,
    const std::vector<std::vector<double>>& sSum2,
    const std::vector<std::vector<double>>& sN,
    const std::string& outname)
{
    setPubStyle();

    const int nE = (int)eEdges.size() - 1;
    const int nS = sSum.empty() ? 0 : (int)sSum[0].size();
    if (nE <= 0 || nS <= 0) return;

    const std::vector<int> col = energyPalette(nE);
    const double w = 1.0 / nS;

    std::vector<std::vector<double>> x(nE), y(nE), ey(nE);
    double ymin = 0.0, ymax = 0.0;
    for (int e = 0; e < nE; ++e)
        for (int k = 0; k < nS; ++k) {
            if (sN[e][k] < 5) continue;
            const double m   = sSum[e][k] / sN[e][k];
            const double var = std::max(0.0, sSum2[e][k] / sN[e][k] - m * m);
            const double sem = std::sqrt(var / sN[e][k]);
            x[e].push_back((k + 0.5) * w + dodge(e, nE, 0.5 * w));
            y[e].push_back(m);
            ey[e].push_back(sem);
            ymin = std::min(ymin, m - sem);
            ymax = std::max(ymax, m + sem);
        }
    if (ymax - ymin <= 0.0) { ymin = -0.01; ymax = 0.01; }
    const double mg = 0.08 * (ymax - ymin);

    TCanvas* c = new TCanvas(uniqueName("c_boost").c_str(), "", 820, 640);
    stylePad(c);
    reserveTag(c);

    const int    ncols = nE > 12 ? 3 : (nE > 4 ? 2 : 1);
    const int    rows  = (nE + ncols - 1) / ncols;
    const double rowH  = 0.046;
    const double ylo   = ymin - mg;
    const double yTop  = withHeadroom(ylo, ymax + mg, legendFrac(c, rows, rowH));

    TH1F* fr = c->DrawFrame(0.0, ylo, 1.0, yTop);
    styleFrame(fr, "|cos#theta*|", "#LT|cos#theta| #minus |cos#theta*|#GT");
    drawHLine(0.0, 1.0, 0.0);

    TLegend* lg = topLegend(c, rows, rowH, 0.034, 0.0, ncols == 1 ? 0.55 : 1.0);
    lg->SetNColumns(ncols);

    for (int e = 0; e < nE; ++e) {
        if (x[e].empty()) continue;
        TGraphErrors* gr = new TGraphErrors(
            (int)x[e].size(), x[e].data(), y[e].data(), nullptr, ey[e].data());
        const int mk = kSeriesMarkers[e % 6];
        styleGraph(gr, col[e], mk, 0.95 * markerScale(mk), 2);
        gr->Draw("LP");
        lg->AddEntry(gr, energyRange(eEdges[e], eEdges[e+1]).c_str(), "lp");
    }
    lg->Draw();

    drawTag();
    c->RedrawAxis();
    saveBoth(c, outname);
    delete c;
}

// ============================================================================
//  One target: boost-corrected and reconstructed W(0)/W(90) on the same
//  panel, CM/lab ratio below (same convention as the Au-vs-U comparison).
//  CM and lab come from the SAME events, so the quadrature error on the
//  ratio is a conservative upper bound.
// ============================================================================
static void plotCMvsLab(
    const std::vector<double>& eEdges,
    const std::vector<LegendreResult>& fitC,
    const std::vector<LegendreResult>& fitL,
    const std::string& outname)
{
    setPubStyle();

    const int nE = static_cast<int>(eEdges.size()) - 1;
    if (nE <= 0) return;

    std::vector<double> xC, yC, eC, xL, yL, eL, xR, yR, eR, xRl, xRh;
    for (int e = 0; e < nE; ++e) {
        const double E1 = eEdges[e], E2 = eEdges[e + 1];
        const double xc = std::sqrt(E1 * E2);
        const double f  = std::pow(E2 / E1, 0.07);

        const bool okC = goodFit(fitC[e]);
        const bool okL = goodFit(fitL[e]);
        if (okC) { xC.push_back(xc / f); yC.push_back(fitC[e].anisotropy); eC.push_back(fitC[e].u_anisotropy); }
        if (okL) { xL.push_back(xc * f); yL.push_back(fitL[e].anisotropy); eL.push_back(fitL[e].u_anisotropy); }

        if (okC && okL) {
            const double C = fitC[e].anisotropy, uC = fitC[e].u_anisotropy;
            const double L = fitL[e].anisotropy, uL = fitL[e].u_anisotropy;
            const double R  = C / L;
            const double uR = R * std::sqrt(std::pow(uL / L, 2) + std::pow(uC / C, 2));
            if (std::isfinite(R) && std::isfinite(uR)) {
                xR.push_back(xc);  yR.push_back(R);  eR.push_back(uR);
                xRl.push_back(xc - E1);  xRh.push_back(E2 - xc);
            }
        }
    }
    if (xC.empty() && xL.empty()) return;

    const double xmin = eEdges.front(), xmax = eEdges.back();
    RatioLayout P = makeRatioLayout("c_cm_vs_lab", true);

    // ── main panel ─────────────────────────────────────────────────────────
    P.top->cd();
    double ymin = 1.0, ymax = 1.0;
    for (size_t i = 0; i < yC.size(); ++i) { ymin = std::min(ymin, yC[i] - eC[i]); ymax = std::max(ymax, yC[i] + eC[i]); }
    for (size_t i = 0; i < yL.size(); ++i) { ymin = std::min(ymin, yL[i] - eL[i]); ymax = std::max(ymax, yL[i] + eL[i]); }
    const double m    = 0.08 * std::max(ymax - ymin, 0.02);
    const double rowH = 0.052 * P.sT;
    const double ylo  = std::max(0.0, ymin - m);
    const double yTop = withHeadroom(ylo, ymax + m, legendFrac(P.top, 1, rowH));

    TH1F* fr = P.top->DrawFrame(xmin, ylo, xmax, yTop);
    styleFrame(fr, "", "W(0^{#circ}) / W(90^{#circ})", P.sT, 1.40);
    fr->GetXaxis()->SetLabelSize(0.0);
    fr->GetXaxis()->SetTitleSize(0.0);
    fr->GetYaxis()->SetNdivisions(505);
    drawHLine(xmin, xmax, 1.0, 7, kRefGray, 2);

    TGraphErrors* gCM = nullptr;
    TGraphErrors* gLab = nullptr;
    if (!xC.empty()) {
        gCM = new TGraphErrors((int)xC.size(), xC.data(), yC.data(), nullptr, eC.data());
        styleGraph(gCM, kAnisoColor, 20, 1.15, 2);
        gCM->Draw("P");
    }
    if (!xL.empty()) {
        gLab = new TGraphErrors((int)xL.size(), xL.data(), yL.data(), nullptr, eL.data());
        styleGraph(gLab, kFitColor, 21, 1.05, 2);
        gLab->Draw("P");
    }

    TLegend* leg = topLegend(P.top, 1, rowH, 0.036 * P.sT, 0.0, 1.0);
    leg->SetNColumns(2);
    leg->SetMargin(0.14);
    if (gCM)  leg->AddEntry(gCM,  "Boost-corrected (CM)", "pe");
    if (gLab) leg->AddEntry(gLab, "Reconstructed (lab)",  "pe");
    leg->Draw();

    padHeader("", figureTag(), 0.040 * P.sT);
    P.top->RedrawAxis();

    // ── ratio panel ────────────────────────────────────────────────────────
    P.bot->cd();
    double rMin = 0.98, rMax = 1.02;
    for (size_t i = 0; i < yR.size(); ++i) {
        rMin = std::min(rMin, yR[i] - eR[i]);
        rMax = std::max(rMax, yR[i] + eR[i]);
    }
    const double rm = 0.15 * (rMax - rMin);
    rMin -= rm;
    rMax += rm;

    TH1F* fr2 = P.bot->DrawFrame(xmin, rMin, xmax, rMax);
    styleFrame(fr2, "E_{n} (MeV)", "CM / lab", P.sB, 1.40, 1.05);
    fr2->GetYaxis()->SetNdivisions(504);
    fr2->GetYaxis()->CenterTitle();
    logXLabels(P.bot, fr2, xmin, xmax, 0.045 * P.sB, 0.012 * P.sB);

    drawBand(xmin, xmax, 0.98, 1.02);
    drawHLine(xmin, xmax, 1.0, 7, kRefGray, 2);

    if (!xR.empty()) {
        TGraphAsymmErrors* gR = new TGraphAsymmErrors(
            (int)xR.size(), xR.data(), yR.data(), xRl.data(), xRh.data(), eR.data(), eR.data());
        styleGraph(gR, kInk, 20, 1.0, 2);
        gR->Draw("P");
    }
    P.bot->RedrawAxis();

    P.c->cd();
    P.c->Modified();
    P.c->Update();
    saveBoth(P.c, outname);
    delete P.c;
}

// ============================================================================
//  LMT figures for one target
// ============================================================================
//    points  -> mean (horizontal bar = energy bin, vertical bar = s.e.m.)
//    dashed  -> median
//    band    -> mean +- 1 sigma, event-by-event standard deviation (NOT an
//               uncertainty on the mean)
//    grey    -> full momentum transfer (LMT = 1)
static void plotLMT(const std::vector<LmtBin>& bins, const std::string& outname)
{
    setPubStyle();

    std::vector<const LmtBin*> ok;
    for (const auto& b : bins) if (okLmtBin(b) && std::isfinite(b.q50)) ok.push_back(&b);
    if (ok.empty()) { std::cerr << "[WARN] plotLMT: nothing to draw\n"; return; }

    const double xmin = bins.front().lo, xmax = bins.back().hi;
    const int n = (int)ok.size();

    TGraphAsymmErrors* gMean = new TGraphAsymmErrors(n);
    TGraph* gMed  = new TGraph(2 * n);
    TGraph* gBand = sigmaBand(ok, kAnisoFill);

    double ymin = 1.0, ymax = 1.0;
    for (int i = 0; i < n; ++i) {
        const LmtBin& b = *ok[i];
        gMean->SetPoint(i, b.xc, b.mean);
        gMean->SetPointError(i, b.xc - b.lo, b.hi - b.xc, b.sem, b.sem);
        gMed->SetPoint(2 * i,     b.lo, b.q50);
        gMed->SetPoint(2 * i + 1, b.hi, b.q50);
        ymin = std::min({ymin, b.mean - b.rms, b.q50});
        ymax = std::max({ymax, b.mean + b.rms, b.q50});
    }
    const double m = 0.08 * std::max(ymax - ymin, 0.02);

    TCanvas* c = new TCanvas(uniqueName("c_lmt").c_str(), "", 820, 640);
    stylePad(c);
    reserveTag(c);
    c->SetLogx();

    const int    rows = 2;
    const double rowH = 0.052;
    const double ylo  = ymin - m;
    const double yTop = withHeadroom(ylo, ymax + m, legendFrac(c, rows, rowH));

    TH1F* fr = c->DrawFrame(xmin, ylo, xmax, yTop);
    styleFrame(fr, "E_{n} (MeV)", "LMT fraction");
    logXLabels(c, fr, xmin, xmax);

    gBand->Draw("F");
    TLine* ref = drawHLine(xmin, xmax, 1.0, 7, kRefGray, 2);

    gMed->SetLineColor(kAnisoColor);
    gMed->SetLineStyle(2);
    gMed->SetLineWidth(2);
    gMed->Draw("L");

    styleGraph(gMean, kAnisoColor, 20, 1.1, 2);
    gMean->Draw("P");

    TLegend* lg = topLegend(c, rows, rowH, 0.036, 0.0, 0.95);
    lg->SetNColumns(2);
    lg->AddEntry(gMean, "Mean #pm s.e.m.",      "pe");
    lg->AddEntry(gMed,  "Median",               "l");
    lg->AddEntry(gBand, "#pm1#sigma (events)",  "f");
    lg->AddEntry(ref,   "Full transfer",        "l");
    lg->Draw();

    drawTag();
    c->RedrawAxis();
    saveBoth(c, outname);
    delete c;
}

// sequential palette white -> light blue -> deep blue -> ink (house colours)
static void setHousePalette()
{
    Double_t stops[5] = {0.00, 0.15, 0.45, 0.75, 1.00};
    Double_t r[5]     = {1.00, 0.81, 0.36, 0.12, 0.05};
    Double_t g[5]     = {1.00, 0.88, 0.55, 0.31, 0.12};
    Double_t b[5]     = {1.00, 0.93, 0.72, 0.47, 0.20};
    TColor::CreateGradientColorTable(5, stops, r, g, b, 255);
    gStyle->SetNumberContours(255);
}

// Density map, each energy column normalised to its maximum; mean overlaid.
static void plotLMTMap(const std::vector<float>& E, const std::vector<float>& F,
                       const std::vector<LmtBin>& bins, const std::string& outname)
{
    setPubStyle();
    if (E.empty() || bins.empty()) return;

    const double xmin = bins.front().lo, xmax = bins.back().hi;
    const std::vector<double> xe = logEdges(xmin, xmax, 48);

    std::vector<float> tmp;
    tmp.reserve(F.size());
    for (size_t i = 0; i < F.size(); ++i)
        if (E[i] >= xmin && E[i] < xmax) tmp.push_back(F[i]);
    if (tmp.size() < 100) return;
    const size_t iLo = (size_t)(0.002 * (tmp.size() - 1));
    const size_t iHi = (size_t)(0.998 * (tmp.size() - 1));
    std::nth_element(tmp.begin(), tmp.begin() + iLo, tmp.end());
    const double qLo = tmp[iLo];
    std::nth_element(tmp.begin(), tmp.begin() + iHi, tmp.end());
    const double qHi = tmp[iHi];
    std::vector<float>().swap(tmp);

    double yLoD = std::min(qLo, 1.0), yHiD = std::max(qHi, 1.0);
    const double mg = 0.05 * std::max(yHiD - yLoD, 0.02);
    yLoD -= mg;
    yHiD += mg;

    TCanvas* c = new TCanvas(uniqueName("c_lmt_map").c_str(), "", 860, 640);
    stylePad(c);
    c->SetRightMargin(0.17);
    reserveTag(c);
    c->SetLogx();

    const double rowH = 0.052;
    const double yTop = withHeadroom(yLoD, yHiD, legendFrac(c, 1, rowH));
    const int    nyD  = 100;
    const double dy   = (yHiD - yLoD) / nyD;
    const int    ny   = (int)std::ceil((yTop - yLoD) / dy);

    TH2D* h = new TH2D(uniqueName("h_lmt_map").c_str(), "",
                       (int)xe.size() - 1, xe.data(), ny, yLoD, yLoD + ny * dy);
    h->SetDirectory(nullptr);
    for (size_t i = 0; i < E.size(); ++i) h->Fill(E[i], F[i]);

    for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
        double mx = 0.0;
        for (int iy = 1; iy <= ny; ++iy) mx = std::max(mx, h->GetBinContent(ix, iy));
        if (mx <= 0.0) continue;
        for (int iy = 1; iy <= ny; ++iy)
            h->SetBinContent(ix, iy, h->GetBinContent(ix, iy) / mx);
    }
    h->SetMinimum(1e-4);
    h->SetMaximum(1.0);

    setHousePalette();
    h->SetContour(255);

    styleFrame(h, "E_{n} (MeV)", "LMT fraction");
    TAxis* az = h->GetZaxis();
    az->SetTitle("Normalised density");
    az->SetTitleFont(42);
    az->SetLabelFont(42);
    az->SetTitleSize(0.045);
    az->SetLabelSize(0.040);
    az->SetTitleOffset(1.30);
    az->SetNdivisions(505);
    az->SetAxisColor(kInk);
    az->SetLabelColor(kInk);
    az->SetTitleColor(kInk);

    h->Draw("COLZ");
    logXLabels(c, h, xmin, xmax);

    TLine* ref = drawHLine(xmin, xmax, 1.0, 7, kRefGray, 2);

    std::vector<double> mx, my;
    for (const auto& b : bins)
        if (okLmtBin(b)) { mx.push_back(b.xc); my.push_back(b.mean); }
    TGraph* gMean = nullptr;
    if (!mx.empty()) {
        TGraph* gHalo = new TGraph((int)mx.size(), mx.data(), my.data());
        gHalo->SetMarkerStyle(20);
        gHalo->SetMarkerSize(1.45);
        gHalo->SetMarkerColor(kWhite);
        gHalo->SetLineColor(kWhite);
        gHalo->SetLineWidth(4);
        gHalo->Draw("LP");

        gMean = new TGraph((int)mx.size(), mx.data(), my.data());
        styleGraph(gMean, kFitColor, 20, 1.0, 2);
        gMean->Draw("LP");
    }

    TLegend* lg = topLegend(c, 1, rowH, 0.036, 0.0, 0.80);
    lg->SetNColumns(2);
    if (gMean) lg->AddEntry(gMean, "Mean", "lp");
    lg->AddEntry(ref, "Full transfer", "l");
    lg->Draw();

    drawTag();

    c->Update();
    if (auto* pal = dynamic_cast<TPaletteAxis*>(h->GetListOfFunctions()->FindObject("palette"))) {
        const double x1 = 1.0 - c->GetRightMargin() + 0.015;
        pal->SetX1NDC(x1);
        pal->SetX2NDC(x1 + 0.030);
        pal->SetY1NDC(c->GetBottomMargin());
        pal->SetY2NDC(1.0 - c->GetTopMargin());
        c->Modified();
        c->Update();
    }

    c->RedrawAxis();
    saveBoth(c, outname);
    delete c;
}

// ============================================================================
//  LMT reference data (Fatyga et al., PRC 32, 1496 (1985): p + 238U)
//  File: E_MeV, value, error [, unit]   (unit = frac | MeVc), '#' = comment
// ============================================================================
struct LmtRefSource {
    std::string path;
    std::string label;
    double      projMass = 938.272;     // MeV/c^2 (proton)
};

struct LmtRef {
    TGraphErrors* graph = nullptr;
    std::string   label;
};

static LmtRef loadLmtRef(const LmtRefSource& src)
{
    LmtRef out;
    out.label = src.label;

    std::ifstream f(src.path);
    if (!f.is_open()) {
        std::cerr << "[WARN] LMT reference file not found: " << src.path << "\n";
        return out;
    }

    std::vector<double> x, y, ey;
    std::string line;
    while (std::getline(f, line)) {
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        std::replace(line.begin(), line.end(), ',', ' ');
        std::replace(line.begin(), line.end(), ';', ' ');

        std::istringstream is(line);
        std::vector<std::string> tok;
        std::string s;
        while (is >> s) tok.push_back(s);
        if (tok.size() < 3) continue;

        try {
            const double E = std::stod(tok[0]);
            double v = std::stod(tok[1]);
            double e = std::stod(tok[2]);
            const std::string unit = (tok.size() > 3) ? tok[3] : "frac";
            if (unit == "MeVc" || unit == "MeV/c") {
                const double pb = std::sqrt(E * E + 2.0 * E * src.projMass);
                v /= pb;
                e /= pb;
            }
            if (!(E > 0.0) || !std::isfinite(v) || !std::isfinite(e)) continue;
            x.push_back(E);
            y.push_back(v);
            ey.push_back(std::fabs(e));
        } catch (...) { continue; }
    }

    if (x.empty()) {
        std::cerr << "[WARN] " << src.path << ": no numeric rows\n";
        return out;
    }
    out.graph = new TGraphErrors((int)x.size(), x.data(), y.data(), nullptr, ey.data());
    return out;
}

static void compareLmtWithRef(const std::vector<float>& E, const std::vector<float>& F,
                              const LmtRef& ref, const std::string& who, double w = 0.10)
{
    if (!ref.graph) return;
    std::printf("\nLMT comparison, %s vs %s  (MC window: E_ref/%.2f .. E_ref*%.2f)\n",
                who.c_str(), ref.label.c_str(), 1.0 + w, 1.0 + w);
    std::printf("  %8s  %16s  %24s  %8s  %10s\n",
                "E (MeV)", "reference", "MC mean +- sem  (rms)", "N_MC", "diff/sigma");

    for (int i = 0; i < ref.graph->GetN(); ++i) {
        const double Er = ref.graph->GetX()[i];
        const double yr = ref.graph->GetY()[i];
        const double er = ref.graph->GetEY()[i];
        const double lo = Er / (1.0 + w), hi = Er * (1.0 + w);

        double s = 0.0, s2 = 0.0;
        Long64_t n = 0;
        for (size_t k = 0; k < E.size(); ++k)
            if (E[k] >= lo && E[k] < hi) { s += F[k]; s2 += double(F[k]) * F[k]; ++n; }

        if (n < 2) {
            std::printf("  %8.1f  %7.3f +- %5.3f  %24s  %8lld  %10s\n",
                        Er, yr, er, "no MC events", (long long)n, "--");
            continue;
        }
        const double m   = s / n;
        const double rms = std::sqrt(std::max(0.0, s2 / n - m * m) * n / (n - 1.0));
        const double sem = rms / std::sqrt((double)n);
        const double sig = std::sqrt(er * er + sem * sem);
        std::printf("  %8.1f  %7.3f +- %5.3f  %7.4f +- %6.4f (%5.3f)  %8lld  %10.2f\n",
                    Er, yr, er, m, sem, rms, (long long)n,
                    sig > 0.0 ? (m - yr) / sig : 0.0);
    }
}

// One target vs published reference. The published errors are uncertainties
// on the mean; the shaded band is the MC event spread (+-1 sigma).
static void plotLMTvsRef(const std::vector<LmtBin>& bins,
                         const std::vector<LmtRef>& refs,
                         const std::string& mcLabel,
                         const std::string& outname)
{
    setPubStyle();

    std::vector<const LmtBin*> ok;
    for (const auto& b : bins) if (okLmtBin(b)) ok.push_back(&b);
    if (ok.empty()) { std::cerr << "[WARN] plotLMTvsRef: no MC bins\n"; return; }

    double xmin = bins.front().lo, xmax = bins.back().hi;
    double ymin = 1.0, ymax = 1.0;

    const int n = (int)ok.size();
    TGraphAsymmErrors* gMean = new TGraphAsymmErrors(n);
    for (int i = 0; i < n; ++i) {
        const LmtBin& b = *ok[i];
        gMean->SetPoint(i, b.xc, b.mean);
        gMean->SetPointError(i, b.xc - b.lo, b.hi - b.xc, b.sem, b.sem);
        ymin = std::min(ymin, b.mean - b.rms);
        ymax = std::max(ymax, b.mean + b.rms);
    }
    TGraph* gBand = sigmaBand(ok, kAnisoFill);

    int nRef = 0;
    for (const auto& r : refs) {
        if (!r.graph) continue;
        ++nRef;
        for (int i = 0; i < r.graph->GetN(); ++i) {
            const double x = r.graph->GetX()[i], y = r.graph->GetY()[i], e = r.graph->GetEY()[i];
            xmin = std::min(xmin, x / 1.25);
            xmax = std::max(xmax, x * 1.25);
            ymin = std::min(ymin, y - e);
            ymax = std::max(ymax, y + e);
        }
    }
    const double m = 0.08 * std::max(ymax - ymin, 0.02);

    TCanvas* c = new TCanvas(uniqueName("c_lmt_ref").c_str(), "", 820, 640);
    stylePad(c);
    reserveTag(c);
    c->SetLogx();

    const int    rows = 3 + nRef;
    const double rowH = 0.050;
    const double ylo  = ymin - m;
    const double yTop = withHeadroom(ylo, ymax + m, legendFrac(c, rows, rowH));

    TH1F* fr = c->DrawFrame(xmin, ylo, xmax, yTop);
    styleFrame(fr, "E_{n}, E_{p} (MeV)", "#LTp_{#parallel}#GT / p_{beam}");
    logXLabels(c, fr, xmin, xmax);

    gBand->Draw("F");
    TLine* ref1 = drawHLine(xmin, xmax, 1.0, 7, kRefGray, 2);

    styleGraph(gMean, kAnisoColor, 20, 1.1, 2);
    gMean->Draw("P");

    TLegend* lg = topLegend(c, rows, rowH, 0.034, 0.0, 0.75);
    lg->AddEntry(gMean, (mcLabel + ": mean #pm s.e.m.").c_str(), "pe");
    lg->AddEntry(gBand, (mcLabel + ": #pm1#sigma (events)").c_str(), "f");

    const std::vector<int> pal = exforPalette();
    int k = 0;
    for (const auto& r : refs) {
        if (!r.graph) continue;
        styleGraph(r.graph, pal[k % pal.size()], kOpenMarkers[k % 6], 1.4, 2);
        r.graph->Draw("P");
        lg->AddEntry(r.graph, r.label.c_str(), "pe");
        ++k;
    }
    lg->AddEntry(ref1, "Full transfer", "l");
    lg->Draw();

    drawTag();
    c->RedrawAxis();
    saveBoth(c, outname);
    delete c;
}

// ============================================================================
//  TARGET COMPARISONS
// ============================================================================

// dodge factor in log E, same for every bin: a fixed fraction of the axis
static double logDodge(int i, int n, double xmin, double xmax, double frac = 0.014)
{
    if (n <= 1) return 1.0;
    return std::pow(xmax / xmin, frac * (i - 0.5 * (n - 1)));
}

// ----------------------------------------------------------------------------
//  <LMT> vs E_n for every target (+ reference data), ratio T[1]/T[0] below.
//  Points: mean +- s.e.m.; translucent bands: +-1 sigma event spread.
//  Both targets use the same energy bins, so the ratio is bin by bin; the
//  samples are independent, so the quadrature error is exact.
// ----------------------------------------------------------------------------
static void plotLmtTargets(const std::vector<const TargetResult*>& T,
                           const std::vector<LmtRef>& refs,
                           const std::string& outname)
{
    setPubStyle();
    if (T.empty() || T[0]->lmtFine.empty()) return;

    const int nT = (int)T.size();
    const int nb = (int)T[0]->lmtFine.size();
    for (const auto* t : T)
        if ((int)t->lmtFine.size() != nb) { std::cerr << "[WARN] plotLmtTargets: bin mismatch\n"; return; }

    double xmin = T[0]->lmtFine.front().lo, xmax = T[0]->lmtFine.back().hi;
    int nRef = 0;
    for (const auto& r : refs) {
        if (!r.graph) continue;
        ++nRef;
        for (int i = 0; i < r.graph->GetN(); ++i) {
            xmin = std::min(xmin, r.graph->GetX()[i] / 1.25);
            xmax = std::max(xmax, r.graph->GetX()[i] * 1.25);
        }
    }

    // ── series: mean +- s.e.m. (points) and +-1 sigma (translucent band) ────
    double ymin = 1.0, ymax = 1.0;
    std::vector<TGraphErrors*> g(nT, nullptr);
    std::vector<TGraph*>       band(nT, nullptr);
    for (int t = 0; t < nT; ++t) {
        std::vector<double> x, y, ey;
        std::vector<const LmtBin*> ok;
        const double f = logDodge(t, nT, xmin, xmax);
        for (const auto& b : T[t]->lmtFine) {
            if (!okLmtBin(b)) continue;
            ok.push_back(&b);
            x.push_back(b.xc * f);
            y.push_back(b.mean);
            ey.push_back(b.sem);
            ymin = std::min(ymin, b.mean - b.rms);
            ymax = std::max(ymax, b.mean + b.rms);
        }
        if (x.empty()) continue;
        g[t] = new TGraphErrors((int)x.size(), x.data(), y.data(), nullptr, ey.data());
        styleGraph(g[t], T[t]->spec.color, T[t]->spec.marker, 1.1, 2);
        band[t] = sigmaBand(ok, TColor::GetColorTransparent(T[t]->spec.color, 0.22));
    }
    for (const auto& r : refs) {
        if (!r.graph) continue;
        for (int i = 0; i < r.graph->GetN(); ++i) {
            ymin = std::min(ymin, r.graph->GetY()[i] - r.graph->GetEY()[i]);
            ymax = std::max(ymax, r.graph->GetY()[i] + r.graph->GetEY()[i]);
        }
    }

    // ratio T[1] / T[0]
    std::vector<double> xR, yR, eR, xRl, xRh;
    if (nT >= 2)
        for (int b = 0; b < nb; ++b) {
            const LmtBin& A = T[1]->lmtFine[b];
            const LmtBin& B = T[0]->lmtFine[b];
            if (!okLmtBin(A) || !okLmtBin(B) || B.mean == 0.0 || A.mean == 0.0) continue;
            const double R  = A.mean / B.mean;
            const double uR = R * std::sqrt(std::pow(A.sem / A.mean, 2) + std::pow(B.sem / B.mean, 2));
            xR.push_back(B.xc);  yR.push_back(R);  eR.push_back(uR);
            xRl.push_back(B.xc - B.lo);  xRh.push_back(B.hi - B.xc);
        }

    RatioLayout P = makeRatioLayout("c_lmt_targets", true);

    // ── main panel ─────────────────────────────────────────────────────────
    P.top->cd();
    // legend in two columns: [mean | +-1 sigma] per target, then refs + line
    const int    rows = nT + (nRef + 2) / 2;
    const double rowH = 0.050 * P.sT;
    const double m    = 0.08 * std::max(ymax - ymin, 0.02);
    const double ylo  = ymin - m;
    const double yTop = withHeadroom(ylo, ymax + m, legendFrac(P.top, rows, rowH));

    TH1F* fr = P.top->DrawFrame(xmin, ylo, xmax, yTop);
    styleFrame(fr, "", "#LTp_{#parallel}#GT / p_{beam}", P.sT, 1.40);
    fr->GetXaxis()->SetLabelSize(0.0);
    fr->GetXaxis()->SetTitleSize(0.0);
    fr->GetYaxis()->SetNdivisions(505);

    for (int t = 0; t < nT; ++t) if (band[t]) band[t]->Draw("F");
    TLine* ref1 = drawHLine(xmin, xmax, 1.0, 7, kRefGray, 2);

    TLegend* lg = topLegend(P.top, rows, rowH, 0.034 * P.sT, 0.0, 1.0);
    lg->SetNColumns(2);
    lg->SetMargin(0.12);
    for (int t = 0; t < nT; ++t) {
        if (!g[t]) continue;
        g[t]->Draw("P");
        lg->AddEntry(g[t],    (T[t]->spec.name + " (MC), mean").c_str(), "pe");
        lg->AddEntry(band[t], (T[t]->spec.name + ", #pm1#sigma").c_str(), "f");
    }
    const std::vector<int> pal = exforPalette();
    int k = 0;
    for (const auto& r : refs) {
        if (!r.graph) continue;
        TGraphErrors* gr = (TGraphErrors*)r.graph->Clone();
        styleGraph(gr, pal[k % pal.size()], kOpenMarkers[k % 6], 1.4, 2);
        gr->Draw("P");
        lg->AddEntry(gr, r.label.c_str(), "pe");
        ++k;
    }
    lg->AddEntry(ref1, "Full transfer", "l");
    lg->Draw();

    padHeader("", figureTag(), 0.040 * P.sT);
    P.top->RedrawAxis();

    // ── ratio panel ────────────────────────────────────────────────────────
    P.bot->cd();
    double rMin = 1.0, rMax = 1.0;
    for (size_t i = 0; i < yR.size(); ++i) {
        rMin = std::min(rMin, yR[i] - eR[i]);
        rMax = std::max(rMax, yR[i] + eR[i]);
    }
    const double rm = 0.15 * std::max(rMax - rMin, 0.02);
    rMin -= rm;
    rMax += rm;

    const std::string rTitle = (nT >= 2) ? T[1]->spec.shortName + " / " + T[0]->spec.shortName : "";
    TH1F* fr2 = P.bot->DrawFrame(xmin, rMin, xmax, rMax);
    styleFrame(fr2, nRef > 0 ? "E_{n}, E_{p} (MeV)" : "E_{n} (MeV)", rTitle.c_str(), P.sB, 1.40, 1.05);
    fr2->GetYaxis()->SetNdivisions(504);
    fr2->GetYaxis()->CenterTitle();
    logXLabels(P.bot, fr2, xmin, xmax, 0.045 * P.sB, 0.012 * P.sB);
    drawHLine(xmin, xmax, 1.0, 7, kRefGray, 2);

    if (!xR.empty()) {
        TGraphAsymmErrors* gR = new TGraphAsymmErrors(
            (int)xR.size(), xR.data(), yR.data(), xRl.data(), xRh.data(), eR.data(), eR.data());
        styleGraph(gR, kInk, 20, 1.0, 2);
        gR->Draw("P");
    }
    P.bot->RedrawAxis();

    P.c->cd();
    P.c->Modified();
    P.c->Update();
    saveBoth(P.c, outname);
    delete P.c;
}

// ----------------------------------------------------------------------------
//  Generic "one value per energy bin, one series per target" figure, used for
//  the boost ratio and for the CM anisotropy comparison.
// ----------------------------------------------------------------------------
struct TargetSeries {
    std::string label;
    int color  = kInk;
    int marker = 20;
    std::vector<double> y, ey;
    std::vector<char>   ok;
};

static void plotTargetsVsEnergy(const std::vector<TargetSeries>& S,
                                const std::vector<double>& edges,
                                const std::string& ytitle,
                                double bandLo, double bandHi,
                                const std::string& bandLabel,
                                const std::string& lineLabel,
                                const std::string& outname)
{
    setPubStyle();
    const int nE = (int)edges.size() - 1;
    const int nS = (int)S.size();
    if (nE <= 0 || nS == 0) return;

    const double xmin = edges.front(), xmax = edges.back();

    double ymin = std::min(1.0, bandLo), ymax = std::max(1.0, bandHi);
    std::vector<TGraphErrors*> g(nS, nullptr);
    for (int s = 0; s < nS; ++s) {
        std::vector<double> x, y, ey;
        const double f = logDodge(s, nS, xmin, xmax);
        for (int e = 0; e < nE && e < (int)S[s].y.size(); ++e) {
            if (!S[s].ok[e]) continue;
            x.push_back(std::sqrt(edges[e] * edges[e + 1]) * f);
            y.push_back(S[s].y[e]);
            ey.push_back(S[s].ey[e]);
            ymin = std::min(ymin, S[s].y[e] - S[s].ey[e]);
            ymax = std::max(ymax, S[s].y[e] + S[s].ey[e]);
        }
        if (x.empty()) continue;
        g[s] = new TGraphErrors((int)x.size(), x.data(), y.data(), nullptr, ey.data());
        styleGraph(g[s], S[s].color, S[s].marker, 1.15, 2);
    }
    const double m = 0.08 * std::max(ymax - ymin, 0.02);

    TCanvas* c = new TCanvas(uniqueName("c_targets").c_str(), "", 820, 640);
    stylePad(c);
    reserveTag(c);
    c->SetLogx();

    const int    nEntries = nS + 2;
    const int    rows = (nEntries + 1) / 2;
    const double rowH = 0.052;
    const double ylo  = ymin - m;
    const double yTop = withHeadroom(ylo, ymax + m, legendFrac(c, rows, rowH));

    TH1F* fr = c->DrawFrame(xmin, ylo, xmax, yTop);
    styleFrame(fr, "E_{n} (MeV)", ytitle.c_str());
    logXLabels(c, fr, xmin, xmax);

    TGraph* band = drawBand(xmin, xmax, bandLo, bandHi);
    TLine*  line = drawHLine(xmin, xmax, 1.0, 7, kRefGray, 2);

    TLegend* lg = topLegend(c, rows, rowH, 0.036, 0.0, 1.0);
    lg->SetNColumns(2);
    for (int s = 0; s < nS; ++s) {
        if (!g[s]) continue;
        g[s]->Draw("P");
        lg->AddEntry(g[s], S[s].label.c_str(), "pe");
    }
    lg->AddEntry(line, lineLabel.c_str(), "l");
    lg->AddEntry(band, bandLabel.c_str(), "f");
    lg->Draw();

    drawTag();
    c->RedrawAxis();
    saveBoth(c, outname);
    delete c;
}

// ============================================================================
//  Per-target analysis, output and figures
// ============================================================================
static void analyseTarget(TargetResult& R,
                          const std::vector<double>& ebAni,
                          const std::vector<double>& ebLmt,
                          bool fit_a4)
{
    const std::vector<double>& ebEff = R.spec.ebEff;
    const int nEff = (int)ebEff.size() - 1;
    const int nAni = (int)ebAni.size() - 1;
    const std::string& key = R.spec.key;

    poissonErrors(R.eff_counts);
    poissonErrors(R.ani_reco);
    poissonErrors(R.ani_cm);

    // efficiency (same routine as the data analysis)
    R.eff.assign(nEff, EfficiencyResult());
    for (int e = 0; e < nEff; ++e)
        R.eff[e] = computeEfficiency(nbins_det - 1, nbins_beam, nbins_det,
                                     R.eff_counts.n, R.eff_counts.u_n,
                                     R.acceptance, e);

    // anisotropy: two beam-angle definitions, same eff and acceptance
    AnalysisConfig cfg = R.spec.makeConfig(ebAni, "mc");   // makeUraniumConfig / makeGoldConfig

    R.ani_r.assign(nAni, AnisotropyResult());
    R.ani_c.assign(nAni, AnisotropyResult());
    R.fit_r.assign(nAni, LegendreResult());
    R.fit_c.assign(nAni, LegendreResult());

    std::cout << "\n=== " << key << " ===\n";
    for (int e = 0; e < nAni; ++e) {
        const double Ec = std::sqrt(ebAni[e] * ebAni[e+1]);
        int e_eff = findBin(ebEff, Ec);
        if (e_eff < 0)     e_eff = 0;
        if (e_eff >= nEff) e_eff = nEff - 1;

        R.ani_r[e] = anisotropy(nbins_beam, nbins_det, R.ani_reco.n, R.ani_reco.u_n,
                                R.acceptance, e, R.eff[e_eff].eps, R.eff[e_eff].u_eps, cfg);
        R.ani_c[e] = anisotropy(nbins_beam, nbins_det, R.ani_cm.n, R.ani_cm.u_n,
                                R.acceptance, e, R.eff[e_eff].eps, R.eff[e_eff].u_eps, cfg);

        R.fit_r[e] = legendre_fit(nbins_beam, nbins_det, R.ani_reco.n, R.acceptance, e,
                                  R.eff[e_eff].eps, R.eff[e_eff].u_eps, cfg, fit_a4);
        R.fit_c[e] = legendre_fit(nbins_beam, nbins_det, R.ani_cm.n, R.acceptance, e,
                                  R.eff[e_eff].eps, R.eff[e_eff].u_eps, cfg, fit_a4);

        std::cout << "  ebin " << e << "  E=" << Ec << " MeV";
        if (R.fit_r[e].valid) std::cout << "  A(reco)=" << R.fit_r[e].anisotropy << " +/- " << R.fit_r[e].u_anisotropy;
        else                  std::cout << "  A(reco)=  --";
        if (R.fit_c[e].valid) std::cout << "   A(boost)=" << R.fit_c[e].anisotropy << " +/- " << R.fit_c[e].u_anisotropy;
        else                  std::cout << "   A(boost)=  --";
        std::cout << "\n";
    }

    // LMT
    if (!R.lmtV.empty()) {
        R.lmtFine  = lmtStats(R.lmtE, R.lmtV, ebLmt);
        R.lmtAniso = lmtStats(R.lmtE, R.lmtV, ebAni);
        std::cout << "  LMT fraction per anisotropy bin:\n";
        for (const auto& b : R.lmtAniso) {
            std::cout << "    [" << b.lo << ", " << b.hi << "] MeV  N=" << b.n;
            if (b.n >= 2)
                std::cout << "  mean=" << b.mean << " +/- " << b.sem
                          << "  sigma=" << b.rms
                          << "  median=" << b.q50
                          << "  LMT>1: " << 100.0 * b.fAbove1 << "%";
            std::cout << "\n";
        }
    }
}

static void saveTargetRoot(const TargetResult& R, const std::string& outdir, int nAni)
{
    const int nEff = (int)R.spec.ebEff.size() - 1;
    TFile* fout = TFile::Open((outdir + "output_mc_analysis_" + R.spec.key + ".root").c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) { std::cerr << "Error creating output file for " << R.spec.key << "\n"; return; }

    for (int e = 0; e < nEff; ++e) {
        TH1D* h = new TH1D(Form("heff_ebin%d", e), "", nbins_det, 0, 1);
        for (int i = 0; i < nbins_det; ++i) {
            h->SetBinContent(i+1, R.eff[e].eps[i]);
            h->SetBinError  (i+1, R.eff[e].u_eps[i]);
        }
        h->Write();
    }
    for (int e = 0; e < nAni; ++e) {
        TH1D* h = new TH1D(Form("hc_reco_ebin%d", e), "", nbins_beam, 0, 1);
        TH1D* h2 = new TH1D(Form("hc_cm_ebin%d", e), "", nbins_beam, 0, 1);
        for (int j = 0; j < nbins_beam; ++j) {
            double sr = 0., sc = 0.;
            for (int i = 0; i < nbins_det; ++i) { sr += R.ani_reco.n[e][j][i]; sc += R.ani_cm.n[e][j][i]; }
            h ->SetBinContent(j+1, sr);  h ->SetBinError(j+1, std::sqrt(std::max(0.0, sr)));
            h2->SetBinContent(j+1, sc);  h2->SetBinError(j+1, std::sqrt(std::max(0.0, sc)));
        }
        h->Write();
        h2->Write();
    }
    if (!R.lmtFine.empty())  writeLmtGraphs(R.lmtFine,  "fine");
    if (!R.lmtAniso.empty()) writeLmtGraphs(R.lmtAniso, "anisobins");
    fout->Close();
}

static void targetFigures(const TargetResult& R,
                          const std::vector<double>& ebAni,
                          const std::vector<LmtRef>& refs,
                          const std::string& outdir)
{
    const std::vector<double>& ebEff = R.spec.ebEff;
    const int nEff = (int)ebEff.size() - 1;
    const int nAni = (int)ebAni.size() - 1;
    const std::string k = "_" + R.spec.key;
    setFigureTag(R.spec.tag);

    plotEfficiency(R.eff, nEff, nbins_det, ebEff, outdir + "efficiency_mc" + k + ".pdf");
    plotEfficiencyResolution(R.eff, nEff, nbins_det, ebEff, outdir + "efficiency_resolution_mc" + k + ".pdf");

    plotAnisotropy(R.ani_c, nAni, nbins_beam, ebAni, outdir + "anisotropy_mc_cm"   + k + ".pdf");
    plotAnisotropy(R.ani_r, nAni, nbins_beam, ebAni, outdir + "anisotropy_mc_reco" + k + ".pdf");

    plotAnisotropyFit(R.fit_c, ebAni, outdir + "anisotropy_fit_mc_cm"   + k + ".pdf", 1);
    plotAnisotropyFit(R.fit_r, ebAni, outdir + "anisotropy_fit_mc_reco" + k + ".pdf", 1);

    plotAnisotropyRatioFit(R.fit_c, ebAni, outdir + "anisotropy_ratio_fit_mc_cm"   + k, "Boost-corrected (CM)");
    plotAnisotropyRatioFit(R.fit_r, ebAni, outdir + "anisotropy_ratio_fit_mc_reco" + k, "Reconstructed (lab)");

    plotBoostShift(ebAni, R.sSum, R.sSum2, R.sN, outdir + "boost_axis_mc" + k);
    plotCMvsLab(ebAni, R.fit_c, R.fit_r, outdir + "anisotropy_mc_cm_vs_lab" + k);

    if (!R.lmtFine.empty()) {
        plotLMT(R.lmtFine, outdir + "lmt_fraction_mc" + k);
        plotLMTMap(R.lmtE, R.lmtV, R.lmtFine, outdir + "lmt_fraction_map_mc" + k);

        if (R.spec.compareWithRefs && !refs.empty()) {
            for (const auto& r : refs) compareLmtWithRef(R.lmtE, R.lmtV, r, R.spec.name + " (MC)");
            plotLMTvsRef(R.lmtFine, refs, R.spec.name + ", MC", outdir + "lmt_fraction_vs_fatyga1985_mc" + k);
        }
    }
}

// ============================================================================
void mc_analysis()
{
    const Bool_t wasBatch = gROOT->IsBatch();
    gROOT->SetBatch(kTRUE);

    const std::string infile =
        "/Users/nico/Desktop/Tese/Analysis/montecarlo/data/mc_complete.root";
    const std::string outdir =
        "/Users/nico/Desktop/Tese/Analysis/montecarlo/output/";
    const std::string acceptance_file =
        "/Users/nico/Desktop/Tese/Analysis/cross_section/acceptance_coincidence.csv";

    // ── binning ────────────────────────────────────────────────────────────
    //  anisotropy and LMT bins are shared, so the targets can be compared bin
    //  by bin; efficiency bins follow each target's data analysis (config.h),
    //  capped at the MC range.
    const std::vector<double> energy_bins_eff_u  = {10, 100, 500, 880};
    const std::vector<double> energy_bins_eff_au = {60, 400, 880};
    const std::vector<double> energy_bins_aniso  = {60, 320, 580, 900};
    const std::vector<double> energy_bins_lmt    =
        logEdges(energy_bins_aniso.front(), energy_bins_aniso.back(), 18);
    const int nbins_aniso = (int)energy_bins_aniso.size() - 1;

    const bool   fit_a4  = false;
    const int    kNShift = 20;
    const double kAmpMin = 0.15;

    // LMT selection: target_i and ppac0 of each target (+ Z window). Set to
    // true to also require the fission-pair cuts (coinc_type, amplitudes).
    const bool kLmtUseFissionPairCuts = true;

    const std::vector<LmtRefSource> lmt_refs_src = {
        { "/Users/nico/Desktop/Tese/Analysis/montecarlo/data/lmt_fatyga1985.csv",
          "p + ^{238}U, Fatyga #it{et al.} (1985)", 938.272 }
    };

    // ── targets ────────────────────────────────────────────────────────────
    //  NOTE: both use the same acceptance file. If the Au target sits at a
    //  different position (ppac0 == 7), point its entry to its own file.
    TargetSpec u;
    u.key = "u238";   u.name = "n + ^{238}U";   u.shortName = "U";
    u.tag = "^{238}U(n,f), MC";
    u.target_i = 2;   u.ppac0 = 8;   u.zSumMin = 80;   u.zSumMax = 92;
    u.color = kAnisoColor;   u.marker = 20;
    u.acceptance_file = acceptance_file;
    u.compareWithRefs = true;
    u.ebEff = energy_bins_eff_u;
    u.makeConfig = makeUraniumConfig;

    TargetSpec au;
    au.key = "au197"; au.name = "n + ^{197}Au"; au.shortName = "Au";
    au.tag = "^{197}Au(n,f), MC";
    au.target_i = 1;  au.ppac0 = 7;  au.zSumMin = 67;  au.zSumMax = 79;
    au.color = kGoldColor;   au.marker = 21;
    au.acceptance_file = acceptance_file;
    au.compareWithRefs = false;
    au.ebEff = energy_bins_eff_au;
    au.makeConfig = makeGoldConfig;

    const std::vector<TargetSpec> specs = {u, au};    // U first: ratios are Au / U

    std::vector<TargetResult> T;
    for (const auto& s : specs) {
        TargetResult r;
        r.spec = s;
        Vec2D dOmega_fine;
        if (!loadAcceptanceCSV(s.acceptance_file, dOmega_fine)) {
            std::cerr << "Failed to load acceptance CSV for " << s.key << "\n";
            return;
        }
        r.acceptance = rebin(dOmega_fine);
        r.eff_counts = makeCounts((int)s.ebEff.size() - 1);
        r.ani_reco   = makeCounts(nbins_aniso);
        r.ani_cm     = makeCounts(nbins_aniso);
        r.sSum .assign(nbins_aniso, std::vector<double>(kNShift, 0.0));
        r.sSum2.assign(nbins_aniso, std::vector<double>(kNShift, 0.0));
        r.sN   .assign(nbins_aniso, std::vector<double>(kNShift, 0.0));
        T.push_back(r);
    }

    // ── open tree ──────────────────────────────────────────────────────────
    TFile* fin = TFile::Open(infile.c_str());
    if (!fin || fin->IsZombie()) { std::cerr << "cannot open " << infile << "\n"; return; }
    TTree* t = (TTree*)fin->Get("CoincTree");
    if (!t) { std::cerr << "CoincTree not found\n"; fin->Close(); return; }
    std::cout << "Entries: " << t->GetEntries() << "\n";

    Double_t neutronE, cos_theta, cos_theta_det, cos_theta_cm_i;
    Int_t    coinc_type, target_i;
    Bool_t   hasCM;
    Int_t    ppac0;
    Double_t amp0_c1, amp1_c1, amp0_c2, amp1_c2;
    Int_t    z1, z2;

    t->SetBranchAddress("neutronE",       &neutronE);
    t->SetBranchAddress("cos_theta",      &cos_theta);
    t->SetBranchAddress("cos_theta_det",  &cos_theta_det);
    t->SetBranchAddress("cos_theta_cm_i", &cos_theta_cm_i);
    t->SetBranchAddress("coinc_type",     &coinc_type);
    t->SetBranchAddress("target_i",       &target_i);
    t->SetBranchAddress("hasCM",          &hasCM);
    t->SetBranchAddress("ppac0",          &ppac0);
    t->SetBranchAddress("amp0_c1", &amp0_c1);
    t->SetBranchAddress("amp1_c1", &amp1_c1);
    t->SetBranchAddress("amp0_c2", &amp0_c2);
    t->SetBranchAddress("amp1_c2", &amp1_c2);
    t->SetBranchAddress("Z1", &z1);
    t->SetBranchAddress("Z0", &z2);

    Double_t lmt_d = std::numeric_limits<double>::quiet_NaN();
    Float_t  lmt_f = std::numeric_limits<float>::quiet_NaN();
    bool haveLMT = false, lmtIsFloat = false;
    if (TLeaf* lf = t->GetLeaf("LMT_fraction")) {
        haveLMT    = true;
        lmtIsFloat = std::string(lf->GetTypeName()) == "Float_t";
        if (lmtIsFloat) t->SetBranchAddress("LMT_fraction", &lmt_f);
        else            t->SetBranchAddress("LMT_fraction", &lmt_d);
    } else {
        std::cerr << "[WARN] branch LMT_fraction not found: LMT figures skipped\n";
    }

    // ── loop (one pass, events dispatched to their target) ─────────────────
    for (Long64_t k = 0, n = t->GetEntries(); k < n; ++k) {
        t->GetEntry(k);

        for (auto& R : T) {
            const TargetSpec& S = R.spec;
            if (target_i != S.target_i) continue;
            ++R.nTarget;

            const int zs = z1 + z2;
            if (zs < S.zSumMin || zs > S.zSumMax || z1 < 4 || z2 < 4) { ++R.nZcut; continue; }

            const bool pairCuts = coinc_type == 1 && ppac0 == S.ppac0 &&
                                  !(amp0_c1 < kAmpMin || amp1_c1 < kAmpMin ||
                                    amp0_c2 < kAmpMin || amp1_c2 < kAmpMin);

            if (haveLMT && ppac0 == S.ppac0 && (pairCuts || !kLmtUseFissionPairCuts)) {
                const double lmt = lmtIsFloat ? (double)lmt_f : lmt_d;
                if (std::isfinite(lmt)) { R.lmtE.push_back((float)neutronE); R.lmtV.push_back((float)lmt); }
                else                    ++R.nnoLMT;
            }

            if (!pairCuts) continue;
            ++R.nsel;

            const bool cmOK = hasCM && std::fabs(cos_theta_cm_i) <= 1.0;
            if (!cmOK) ++R.nnoCM;

            const int d = cosBin(cos_theta_det, nbins_det);
            if (d < 0) continue;

            const int b_reco = cosBin(cos_theta,      nbins_beam);
            const int b_cm   = cosBin(cos_theta_cm_i, nbins_beam);

            const int nEffT = (int)S.ebEff.size() - 1;
            const int e_eff = findBin(S.ebEff, neutronE);
            if (e_eff >= 0 && e_eff < nEffT && b_reco >= 0)
                R.eff_counts.n[e_eff][b_reco][d] += 1.;

            const int e_ani = findBin(energy_bins_aniso, neutronE);
            if (e_ani >= 0 && e_ani < nbins_aniso && cmOK && b_reco >= 0 && b_cm >= 0) {
                R.ani_reco.n[e_ani][b_reco][d] += 1.;
                R.ani_cm  .n[e_ani][b_cm  ][d] += 1.;

                const int ks = cosBin(cos_theta_cm_i, kNShift);
                const double dc = std::fabs(cos_theta) - std::fabs(cos_theta_cm_i);
                R.sSum [e_ani][ks] += dc;
                R.sSum2[e_ani][ks] += dc * dc;
                R.sN   [e_ani][ks] += 1.;
            }
        }
    }
    fin->Close();

    for (const auto& R : T) {
        std::cout << R.spec.key << ": " << R.nTarget << " entries on target, "
                  << R.nZcut << " outside Z1+Z2 in [" << R.spec.zSumMin << ", " << R.spec.zSumMax << "], "
                  << R.nsel << " selected fission pairs (" << R.nnoCM << " without CM angle)";
        if (haveLMT) std::cout << ", LMT: " << R.lmtV.size() << " values (" << R.nnoLMT << " non-finite)";
        std::cout << "\n";
    }

    // ── analysis, output, per-target figures ───────────────────────────────
    std::vector<LmtRef> refs;
    for (const auto& src : lmt_refs_src) {
        LmtRef r = loadLmtRef(src);
        if (r.graph) refs.push_back(r);
    }

    for (auto& R : T) {
        analyseTarget(R, energy_bins_aniso, energy_bins_lmt, fit_a4);
        saveTargetRoot(R, outdir, nbins_aniso);
        targetFigures(R, energy_bins_aniso, refs, outdir);
    }

    // ── comparisons Au vs U ────────────────────────────────────────────────
    setFigureTag("MC");
    std::vector<const TargetResult*> TP;
    for (const auto& R : T) TP.push_back(&R);

    // LMT: U first, so the ratio panel is Au / U
    if (haveLMT)
        plotLmtTargets(TP, refs, outdir + "lmt_fraction_au_vs_u_mc");

    // boost correction and boost-corrected anisotropy
    std::vector<TargetSeries> sBoost, sAniso;
    for (const auto& R : T) {
        TargetSeries b, a;
        b.label = a.label = R.spec.name + " (MC)";
        b.color = a.color = R.spec.color;
        b.marker = a.marker = R.spec.marker;
        for (int e = 0; e < nbins_aniso; ++e) {
            const LegendreResult& C = R.fit_c[e];
            const LegendreResult& L = R.fit_r[e];
            const bool okC = goodFit(C), okL = goodFit(L);

            a.y.push_back(okC ? C.anisotropy : 0.0);
            a.ey.push_back(okC ? C.u_anisotropy : 0.0);
            a.ok.push_back(okC);

            double r = 0.0, ur = 0.0;
            bool ok = okC && okL;
            if (ok) {
                r  = C.anisotropy / L.anisotropy;
                ur = r * std::sqrt(std::pow(C.u_anisotropy / C.anisotropy, 2) +
                                   std::pow(L.u_anisotropy / L.anisotropy, 2));
                ok = std::isfinite(r) && std::isfinite(ur);
            }
            b.y.push_back(r);
            b.ey.push_back(ur);
            b.ok.push_back(ok);
        }
        sBoost.push_back(b);
        sAniso.push_back(a);
    }

    plotTargetsVsEnergy(sBoost, energy_bins_aniso,
                        "A_{boost} / A_{reco}", 0.98, 1.02,
                        "#pm2%", "No boost effect",
                        outdir + "boost_ratio_au_vs_u_mc");

    plotTargetsVsEnergy(sAniso, energy_bins_aniso,
                        "W(0^{#circ}) / W(90^{#circ})  (boost-corrected)", 0.95, 1.05,
                        "#pm5% band", "Isotropic",
                        outdir + "anisotropy_cm_au_vs_u_mc");

    setFigureTag("");
    gROOT->SetBatch(wasBatch);
}