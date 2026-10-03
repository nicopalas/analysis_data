// acceptance_mc.C
//
// Geometric acceptance of a PPAC pair for back-to-back fission fragments.
//
// Works in the DETECTOR frame, where the 45-degree tilt is undone and the
// PPAC planes are perpendicular to the local z axis.
//
// FRAME CONSTRUCTION
// ------------------
// The PPAC volumes are rotated +45 deg about Y in the lab, so
//     v_lab = R_y(45) v_local   =>   v_local = R_y(-45) v_lab
// Explicitly, R_y(-45) (X,Y,Z) = ( (X-Z)/sqrt2 , Y , (X+Z)/sqrt2 ).
//
// BEAM AXIS
// ---------
// The beam runs along +z_lab, so in the detector frame its direction is
//     b_hat = R_y(-45) (0,0,1) = (-1, 0, +1)/sqrt2
// i.e. the beam axis is the line  x_det = -z_det , tilted -45 deg about Y
// with respect to the detector z axis.  The PPACs are stacked ALONG THE
// BEAM (z_ppac_lab = z0 + i*d_PPACs, with x_ppac_lab = 0), so every PPAC
// centre lies on that line and
//     x_ppac = -z_ppac                        <-- never a free parameter
// This is enforced below instead of being hard-coded as four independent
// numbers, so the geometry cannot drift out of consistency.
//
// CATHODE PLACEMENT
// -----------------
// CreatePPAC puts the cathodes at local (0, 0, +-(gas_gap + mylar)).  A
// displacement that is already purely local-z maps to detector-frame z with
// NO x component:
//     cathode_lab = (0,0,Delta) + R_y(45)(0,0,+-c)
//     => detector  ( -Delta/sqrt2 , 0 , +Delta/sqrt2 +- c )
// So both cathodes of a PPAC share the x of their PPAC centre and differ
// only in z.  A cathode centre is therefore NOT on the beam axis: it is off
// it by c = 0.32017 cm.  If CreatePPAC instead stacks the cathodes along the
// beam direction, set cathodes_on_beam_axis = true below (costs ~0.7 % of
// absolute coincidence probability).
//
// CATHODE ORDERING  (this is what was wrong before)
// -------------------------------------------------
// The internal layout along the PPAC normal is  Y, anode, X, and ALL PPACs
// are added with the SAME rotation -- the back one is not flipped.  Hence,
// for BOTH PPACs, in detector-frame z:
//     z_Y = z_ppac - c        z_X = z_ppac + c
// The forward fragment (travelling +z_det) meets front-Y then front-X; the
// backward fragment (travelling -z_det) meets back-X then back-Y.  The
// ARRIVAL ORDER is reversed but the PLANE IDENTITY is not: labelling the
// planes by arrival order mirrors the back PPAC and gives unequal lever arms
// (5+2c for X, 5-2c for Y), a spurious +-12.8 % azimuthal distortion of the
// reconstructed direction.  With the correct assignment both lever arms are
// exactly d_PPACs/sqrt2 and the reconstruction closes to machine precision.
//
// COINCIDENCE CONDITION
// ---------------------
// A fragment must cross BOTH cathode planes of its PPAC inside the active
// area (a signal is needed on the X plane and on the Y plane), and the same
// for the back-to-back partner in the opposite PPAC: four intersections.
//
// Geometry from create_ntof_geo.C:
//   target  TGeoEltu semi-axes (7.8*sqrt2/2, 7.8/2) cm, U thickness 0.411e-4 cm
//           (the sqrt2 elongation in local x is the footprint of a round beam
//            on a foil tilted 45 deg; in the detector frame the foil is flat
//            at z_det = 0, so sampling (ox,oy) in that ellipse is correct)
//   PPAC    TGeoBBox half-size 10 x 10 cm, gas_gap 0.32, mylar 1.7e-4
//
// TRUE-vs-RECONSTRUCTED ANGLE CHECK
// ---------------------------------
// dfx,dfy,dfz is the TRUE emission direction in the detector frame (dfz =
// cth by construction).  The reconstructed direction comes from the four
// cathode intersection points.  For an ideal back-to-back pair the two
// cathodes of a PPAC lie on the SAME straight track, so with the correct
// plane assignment the reconstruction is exact and these histograms collapse
// to delta functions at zero -- that is the intended closure test.  To study
// a real reconstruction bias, switch on strip_pitch below (0 = off), which
// digitises the measured coordinates.
//
// FIGURES
// -------
// angle_true_vs_reco.{pdf,png} and costheta_true_vs_reco.{pdf,png}: one
// figure each, beam frame | detector frame, with true vs reconstructed
// distributions, pulls per bin (exact, from bin migrations) and the
// reco - true map.  See the recofig block below.

#include "TRandom3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TStyle.h"
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

// ============================================================================
//  recofig: true vs reconstructed angle, one figure per variable (theta, cos)
//
//  Two columns (beam frame | detector frame) sharing x within each column:
//    top     distributions: true (filled, deep blue) and reconstructed
//            (terracotta, dashed), the two accents of the result figures;
//    middle  pulls per bin, bars in the reconstructed colour, with the
//            +-1 sigma and +-2 sigma bands;
//    bottom  reco - true vs true, viridis, shared log colour scale.
//
//  PULLS.  True and reconstructed histograms are filled with the SAME events,
//  so they are not independent and (N_reco - N_true)/sqrt(N_reco + N_true) is
//  wrong (too small).  An event changes bin i of the difference only if it
//  migrates, i.e. its true and reco values fall in different bins, so
//      N_reco - N_true = n_in - n_out,     sigma = sqrt(n_in + n_out),
//  with n_in / n_out counted in the event loop.  With strip_pitch = 0 nothing
//  migrates and every pull is exactly zero (the closure test), which the
//  figure states instead of showing an empty panel.
//
//  Layout: exact pixel geometry, Helvetica pixel fonts, hand-drawn frames,
//  ticks, tick labels, legend and colour bar (as in the efficiency and cut
//  figures), batch mode so the canvas is exactly kW x kH.
// ============================================================================
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TColor.h"
#include "TAxis.h"
#include "TROOT.h"
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

namespace recofig {

// ---- fonts (Helvetica, precision 3: size in pixels) ------------------------------
const int kFont = 43, kFontBold = 63;
const int kLabelPx  = 40, kTitlePx = 46, kHeadPx = 38, kLegendPx = 40;
const int kFigTitPx = 50, kSubPx   = 34, kNotePx = 32, kBarLabPx = 36;

// ---- geometry, in pixels (checked against measured text widths) ----------------
const int kPxLeft   = 240;   // y labels up to 101 px ("-0.01") + rotated y titles
const int kPxRight  = 220;   // colour bar, its labels and title
const int kPxColW   = 560;
const int kPxColGap =  80;   // "180" | "0" under adjacent columns
const int kPxTop    = 170;   // title row + legend row
const int kPxHead   =  70;   // column headers
const int kPxDistH  = 400;   // distributions
const int kPxGap1   =  20;   // distributions -> pulls (shared x, no labels between)
const int kPxPullH  = 200;   // pulls
const int kPxGap2   =  40;   // pulls -> residual map
const int kPxMapH   = 400;   // residual map
const int kPxBot    = 150;   // x labels + x title

const int kTickMaj = 20, kTickMin = 10, kLabelGap = 14;

const int kDistTop = kPxTop + kPxHead,   kDistBot = kDistTop + kPxDistH;
const int kPullTop = kDistBot + kPxGap1, kPullBot = kPullTop + kPxPullH;
const int kMapTop  = kPullBot + kPxGap2, kMapBot  = kMapTop  + kPxMapH;
const int kW = kPxLeft + 2 * kPxColW + kPxColGap + kPxRight;   // 1660
const int kH = kMapBot + kPxBot;                                // 1450
const int kFramesRightPx = kPxLeft + 2 * kPxColW + kPxColGap;  // 1440

inline int    colLeftPx(int c) { return kPxLeft + c * (kPxColW + kPxColGap); }
inline double ndcX(double px)  { return px / kW; }
inline double ndcY(double px)  { return 1.0 - px / kH; }

// ---- colours -------------------------------------------------------------------------
// true: deep blue #1F4E79 on a pale blue fill; reconstructed and pulls:
// terracotta #D1603D (its complement); map: viridis, as the other 2D figures.
int cInk = 1, cTrueFill = 0, cTrueLine = 1, cReco = 1, cBand1 = 0, cBand2 = 0, cRef = 1, cNote = 1;
const int kStyleReco = 11;   // dashed

inline void setStyle()
{
    cInk      = TColor::GetColor("#1B1B1B");
    cTrueFill = TColor::GetColor("#D9E4F1");
    cTrueLine = TColor::GetColor("#1F4E79");
    cReco     = TColor::GetColor("#D1603D");
    cBand1    = TColor::GetColor("#E3E6EA");   // +-1 sigma
    cBand2    = TColor::GetColor("#F2F3F5");   // +-2 sigma
    cRef      = TColor::GetColor("#8C8C8C");
    cNote     = TColor::GetColor("#6E6E6E");
    gStyle->SetLineStyleString(kStyleReco, "24 12");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kViridis);
    gStyle->SetNumberContours(255);
}

// ---- primitives -------------------------------------------------------------------------
inline TLine* newLine(double x1, double y1, double x2, double y2, int col, int style, int width)
{
    TLine* l = new TLine(x1, y1, x2, y2);
    l->SetLineColor(col); l->SetLineStyle(style); l->SetLineWidth(width);
    l->SetBit(kCanDelete); l->Draw();
    return l;
}

inline TBox* newBox(double x1, double y1, double x2, double y2, int fill, int line = -1)
{
    TBox* b = new TBox(x1, y1, x2, y2);
    b->SetFillColor(fill); b->SetFillStyle(1001); b->SetLineWidth(0);
    b->SetBit(kCanDelete); b->Draw();
    if (line >= 0) {
        TBox* o = new TBox(x1, y1, x2, y2);
        o->SetFillStyle(0); o->SetLineColor(line); o->SetLineWidth(2);
        o->SetBit(kCanDelete); o->Draw();
    }
    return b;
}

// text at canvas pixels (from the left, from the top); gPad = canvas
inline void text(double xpx, double ypx, const std::string& s, int px, int align, int col,
                 double angle = 0.0, int font = kFont)
{
    TLatex t;
    t.SetNDC();
    t.SetTextFont(font); t.SetTextSize(px); t.SetTextAlign(align);
    t.SetTextAngle(angle); t.SetTextColor(col);
    t.DrawLatex(ndcX(xpx), ndcY(ypx), s.c_str());
}

// tick label: %g, typographic minus
inline std::string num(double v)
{
    if (std::fabs(v) < 1e-12) return "0";
    const std::string s = Form("%g", std::round(std::fabs(v) * 1e9) / 1e9);
    return v < 0 ? "#minus" + s : s;
}

// smallest 1, 2, 2.5, 5 x 10^k >= v
inline double niceStep(double v)
{
    const double d = std::pow(10.0, std::floor(std::log10(v)));
    for (double m : {1.0, 2.0, 2.5, 5.0, 10.0}) if (m * d >= v * 0.9999) return m * d;
    return 10.0 * d;
}
inline int minorsFor(double step)
{
    const double m = step / std::pow(10.0, std::floor(std::log10(step) + 1e-9));
    return (std::fabs(m - 2.0) < 1e-6) ? 4 : 5;
}

// ---- axes -------------------------------------------------------------------------------
struct AxisSpec {
    double lo, hi;        // range
    double anchor, major; // majors at anchor + n*major
    int    nMinor;        // minor intervals per major
    bool   edgeLabels;    // label the values that fall on the frame edges
};

inline std::vector<double> majors(const AxisSpec& a)
{
    std::vector<double> v;
    const int n0 = (int)std::ceil ((a.lo - a.anchor) / a.major - 1e-9);
    const int n1 = (int)std::floor((a.hi - a.anchor) / a.major + 1e-9);
    for (int n = n0; n <= n1; ++n) {
        const double x = a.anchor + n * a.major;
        const bool edge = std::fabs(x - a.lo) < 1e-9 * a.major || std::fabs(x - a.hi) < 1e-9 * a.major;
        if (edge && !a.edgeLabels) continue;
        v.push_back(x);
    }
    return v;
}

// a transparent full-width pad whose frame is column c, rows [topPx, botPx]
inline TPad* makeFramePad(const std::string& name, int c, int topPx, int botPx,
                          const AxisSpec& X, const AxisSpec& Y, bool logz = false)
{
    const int padTop = topPx - 2, padBot = botPx + 2;
    const double h = padBot - padTop;
    TPad* p = new TPad(name.c_str(), "", 0.0, ndcY(padBot), 1.0, ndcY(padTop));
    p->SetLeftMargin  (double(colLeftPx(c)) / kW);
    p->SetRightMargin (double(kW - colLeftPx(c) - kPxColW) / kW);
    p->SetTopMargin   (2.0 / h);
    p->SetBottomMargin(2.0 / h);
    p->SetFillStyle(0); p->SetBorderMode(0);
    p->SetFrameLineColor(cInk); p->SetFrameLineWidth(2);
    p->SetTickx(0); p->SetTicky(0);
    if (logz) p->SetLogz();
    p->Draw();
    p->cd();
    TH1F* fr = p->DrawFrame(X.lo, Y.lo, X.hi, Y.hi);
    for (TAxis* a : {fr->GetXaxis(), fr->GetYaxis()}) {
        a->SetLabelSize(0); a->SetTickLength(0); a->SetTitle(""); a->SetAxisColor(cInk);
    }
    return p;
}

// frame + ticks inside on the four sides; gPad = the frame's pad
inline void drawFrameAndTicks(const AxisSpec& X, const AxisSpec& Y, int wpx, int hpx)
{
    const double ux = (X.hi - X.lo) / wpx, uy = (Y.hi - Y.lo) / hpx;
    newLine(X.lo, Y.lo, X.hi, Y.lo, cInk, 1, 2);
    newLine(X.lo, Y.hi, X.hi, Y.hi, cInk, 1, 2);
    newLine(X.lo, Y.lo, X.lo, Y.hi, cInk, 1, 2);
    newLine(X.hi, Y.lo, X.hi, Y.hi, cInk, 1, 2);

    auto ticks = [&](const AxisSpec& a, bool isX) {
        const double st = a.major / a.nMinor;
        const int n0 = (int)std::ceil ((a.lo - a.anchor) / st - 1e-9);
        const int n1 = (int)std::floor((a.hi - a.anchor) / st + 1e-9);
        for (int n = n0; n <= n1; ++n) {
            const double v = a.anchor + n * st;
            if (std::fabs(v - a.lo) < 1e-9 * st || std::fabs(v - a.hi) < 1e-9 * st) continue;
            const int len = (n % a.nMinor == 0) ? kTickMaj : kTickMin;
            if (isX) {
                newLine(v, Y.lo, v, Y.lo + len * uy, cInk, 1, 2);
                newLine(v, Y.hi, v, Y.hi - len * uy, cInk, 1, 2);
            } else {
                newLine(X.lo, v, X.lo + len * ux, v, cInk, 1, 2);
                newLine(X.hi, v, X.hi - len * ux, v, cInk, 1, 2);
            }
        }
    };
    ticks(X, true);
    ticks(Y, false);
}

// tick labels; gPad = canvas
inline void xLabels(int c, const AxisSpec& X, int framBotPx)
{
    for (double v : majors(X)) {
        const double px = colLeftPx(c) + (v - X.lo) / (X.hi - X.lo) * kPxColW;
        text(px, framBotPx + kLabelGap, num(v), kLabelPx, 23, cInk);
    }
}
inline void yLabels(const AxisSpec& Y, int topPx, int botPx)
{
    for (double v : majors(Y)) {
        const double py = botPx - (v - Y.lo) / (Y.hi - Y.lo) * (botPx - topPx);
        text(kPxLeft - kLabelGap, py, num(v), kLabelPx, 32, cInk);
    }
}

// step outline of h * scale over its non-empty range; closed to y = 0 for filling
inline TGraph* stepGraph(const TH1* h, double scale, bool closeToZero)
{
    int b0 = 1, b1 = h->GetNbinsX();
    while (b0 <= b1 && h->GetBinContent(b0) == 0) ++b0;
    while (b1 >= b0 && h->GetBinContent(b1) == 0) --b1;
    TGraph* g = new TGraph();
    if (b0 > b1) return g;
    if (closeToZero) g->SetPoint(g->GetN(), h->GetBinLowEdge(b0), 0.0);
    for (int b = b0; b <= b1; ++b) {
        const double v = h->GetBinContent(b) * scale;
        g->SetPoint(g->GetN(), h->GetBinLowEdge(b),     v);
        g->SetPoint(g->GetN(), h->GetBinLowEdge(b + 1), v);
    }
    if (closeToZero) g->SetPoint(g->GetN(), h->GetBinLowEdge(b1 + 1), 0.0);
    g->SetBit(kCanDelete);
    return g;
}

// log colour bar, drawn by hand on the canvas (gPad = canvas)
inline double niceUp(double v)
{
    const double d = std::pow(10.0, std::floor(std::log10(v)));
    for (double m : {1.0, 2.0, 5.0, 10.0}) if (m * d >= v * 0.9999) return m * d;
    return 10.0 * d;
}
inline void drawColourBar(double zmin, double zmax, int x1, int x2, int top, int bot)
{
    const int ncol = gStyle->GetNumberOfColors(), nbox = 255;
    const double l0 = std::log10(zmin), l1 = std::log10(zmax), h = bot - top;
    for (int k = 0; k < nbox; ++k) {
        const double ya = bot - h * k / nbox;
        const double yb = bot - h * (k + 1) / nbox - (k + 1 < nbox ? 0.6 : 0.0);
        newBox(ndcX(x1), ndcY(ya), ndcX(x2), ndcY(yb),
               gStyle->GetColorPalette(std::min(ncol - 1, (int)((k + 0.5) * ncol / nbox))));
    }
    newLine(ndcX(x1), ndcY(top), ndcX(x2), ndcY(top), cInk, 1, 2);
    newLine(ndcX(x1), ndcY(bot), ndcX(x2), ndcY(bot), cInk, 1, 2);
    newLine(ndcX(x1), ndcY(top), ndcX(x1), ndcY(bot), cInk, 1, 2);
    newLine(ndcX(x2), ndcY(top), ndcX(x2), ndcY(bot), cInk, 1, 2);
    auto ypx = [&](double v) { return bot - (std::log10(v) - l0) / (l1 - l0) * h; };
    std::vector<double> lab;
    for (int k = -1; k <= 10; ++k)
        for (int m = 1; m <= 9; ++m) {
            const double v = m * std::pow(10.0, k);
            if (v < zmin * 0.9999 || v > zmax * 1.0001) continue;
            newLine(ndcX(x2), ndcY(ypx(v)), ndcX(x2 + (m == 1 ? 10 : 5)), ndcY(ypx(v)), cInk, 1, 2);
            if (m == 1) lab.push_back(v);
        }
    if (lab.size() < 3)
        for (int k = -1; k <= 10; ++k)
            for (double m : {2.0, 5.0}) {
                const double v = m * std::pow(10.0, k);
                if (v >= zmin * 0.9999 && v <= zmax * 1.0001) lab.push_back(v);
            }
    for (double v : lab) {
        const int k = (int)std::lround(std::log10(v));
        const bool dec = std::fabs(v - std::pow(10.0, k)) < 1e-6 * v;
        text(x2 + 16, ypx(v), (dec && k >= 2) ? Form("10^{%d}", k) : Form("%g", v), kBarLabPx, 12, cInk);
    }
}

// ---- one column of the figure -------------------------------------------------------
struct Column {
    TH1D *hTrue, *hReco;   // distributions (display binning)
    TH1D *hIn,   *hOut;    // migrations into / out of each bin
    TH2D *hMap;            // reco - true vs true
    AxisSpec X, Ymap;
    std::string header, xTitle;
};

// ---- the figure ---------------------------------------------------------------------------
inline void drawFigure(const Column cols[2], const std::string& binUnit,
                       const std::string& mapTitle, const std::string& subtitle,
                       const std::string& outbase)
{
    setStyle();

    // distributions: common y, scaled by 10^k (k multiple of 3) for short labels
    double dmax = 0.0;
    for (int c = 0; c < 2; ++c)
        dmax = std::max({dmax, cols[c].hTrue->GetMaximum(), cols[c].hReco->GetMaximum()});
    if (dmax <= 0.0) dmax = 1.0;
    const int    kExp  = (dmax >= 1e3) ? 3 * (int)std::floor(std::log10(dmax) / 3.0) : 0;
    const double scale = std::pow(10.0, -kExp);
    const double yStep = niceStep(1.08 * dmax * scale / 4.0);
    const AxisSpec Ydist = {0.0, std::ceil(1.08 * dmax * scale / yStep - 1e-9) * yStep,
                            0.0, yStep, minorsFor(yStep), true};

    // pulls: exact, from migrations; common range
    std::vector<std::vector<double>> pull(2);
    double pmax = 0.0;
    for (int c = 0; c < 2; ++c) {
        const Column& k = cols[c];
        for (int b = 1; b <= k.hTrue->GetNbinsX(); ++b) {
            const double d = k.hReco->GetBinContent(b) - k.hTrue->GetBinContent(b);
            const double s = std::sqrt(k.hIn->GetBinContent(b) + k.hOut->GetBinContent(b));
            const double p = (s > 0.0) ? d / s : 0.0;
            pull[c].push_back(p);
            pmax = std::max(pmax, std::fabs(p));
        }
    }
    // labels at 0 and +-step only, with the range = 1.5 step: the labels sit a
    // third of the way in from the edges and never meet the panels around
    double pStep = 1000.0;
    for (double s : {3.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0})
        if (1.5 * s >= 1.15 * pmax) { pStep = s; break; }
    const double R = std::max(5.0, 1.5 * pStep);
    const int pMinor = (pStep == 3.0) ? 3 : (pStep == 20.0 ? 4 : 5);
    const AxisSpec Ypull = {-R, R, 0.0, pStep, pMinor, false};

    // residual maps: common log colour scale
    double zmax = 2.0;
    for (int c = 0; c < 2; ++c) zmax = std::max(zmax, cols[c].hMap->GetMaximum());
    zmax = niceUp(zmax);
    const double zmin = 0.9;

    const Bool_t wasBatch = gROOT->IsBatch();
    gROOT->SetBatch(kTRUE);
    TCanvas* cv = new TCanvas(("c_" + outbase).c_str(), "", kW, kH);
    cv->SetCanvasSize(kW, kH);
    cv->SetFillColor(kWhite);
    cv->SetBorderMode(0);

    for (int c = 0; c < 2; ++c) {
        const Column& k = cols[c];

        // -- distributions --
        cv->cd();
        makeFramePad(Form("pd_%d", c), c, kDistTop, kDistBot, k.X, Ydist);
        TGraph* gf = stepGraph(k.hTrue, scale, true);
        gf->SetFillColor(cTrueFill); gf->SetFillStyle(1001);
        if (gf->GetN() > 2) gf->Draw("F");
        TGraph* gt = stepGraph(k.hTrue, scale, true);
        gt->SetLineColor(cTrueLine); gt->SetLineWidth(3);
        if (gt->GetN() > 1) gt->Draw("L");
        TGraph* gr = stepGraph(k.hReco, scale, false);
        gr->SetLineColor(cReco); gr->SetLineWidth(4); gr->SetLineStyle(kStyleReco);
        if (gr->GetN() > 1) gr->Draw("L");
        drawFrameAndTicks(k.X, Ydist, kPxColW, kPxDistH);

        // -- pulls --
        cv->cd();
        makeFramePad(Form("pp_%d", c), c, kPullTop, kPullBot, k.X, Ypull);
        newBox(k.X.lo, -2.0, k.X.hi, 2.0, cBand2);
        newBox(k.X.lo, -1.0, k.X.hi, 1.0, cBand1);
        for (int b = 1; b <= k.hTrue->GetNbinsX(); ++b) {
            const double p = std::max(-R, std::min(R, pull[c][b - 1]));
            if (p != 0.0)
                newBox(k.hTrue->GetBinLowEdge(b), std::min(0.0, p),
                       k.hTrue->GetBinLowEdge(b + 1), std::max(0.0, p), cReco);
        }
        newLine(k.X.lo, 0.0, k.X.hi, 0.0, cInk, 1, 2);
        drawFrameAndTicks(k.X, Ypull, kPxColW, kPxPullH);

        // -- residual map --
        cv->cd();
        makeFramePad(Form("pm_%d", c), c, kMapTop, kMapBot, k.X, k.Ymap, true);
        k.hMap->SetMinimum(zmin);
        k.hMap->SetMaximum(zmax);
        k.hMap->SetStats(0);
        k.hMap->Draw("COL SAME");
        newLine(k.X.lo, 0.0, k.X.hi, 0.0, cRef, 1, 1);
        drawFrameAndTicks(k.X, k.Ymap, kPxColW, kPxMapH);
    }

    // ---- text, legend, colour bar: on the canvas ------------------------------------
    cv->cd();
    for (int c = 0; c < 2; ++c) {
        const int x0 = colLeftPx(c);
        const double yh = kDistTop - 0.5 * kPxHead;
        text(x0,      yh, Form("(%c)", 'a' + c), kHeadPx, 12, cInk, 0.0, kFontBold);
        text(x0 + 60, yh, cols[c].header,         kHeadPx, 12, cInk);
        xLabels(c, cols[c].X, kMapBot);
        text(x0 + 0.5 * kPxColW, kMapBot + 98, cols[c].xTitle, kTitlePx, 22, cInk);
    }
    if (pmax == 0.0)
        for (int c = 0; c < 2; ++c)
            text(colLeftPx(c) + 0.5 * kPxColW, kPullTop + 0.25 * kPxPullH,
                 "All pulls = 0 (exact closure)", kNotePx, 22, cNote);

    yLabels(Ydist, kDistTop, kDistBot);
    yLabels(Ypull, kPullTop, kPullBot);
    yLabels(cols[0].Ymap, kMapTop, kMapBot);
    const std::string yDist = "Events / " + binUnit + (kExp > 0 ? Form(" (10^{%d})", kExp) : "");
    text(62, 0.5 * (kDistTop + kDistBot), yDist,    kTitlePx, 22, cInk, 90.0);
    text(62, 0.5 * (kPullTop + kPullBot), "Pull",   kTitlePx, 22, cInk, 90.0);
    text(62, 0.5 * (kMapTop  + kMapBot),  mapTitle, kTitlePx, 22, cInk, 90.0);

    drawColourBar(zmin, zmax, kFramesRightPx + 24, kFramesRightPx + 58, kMapTop, kMapBot);
    text(kFramesRightPx + 180, 0.5 * (kMapTop + kMapBot), "Events", kTitlePx, 22, cInk, 90.0);

    // title row
    text(kPxLeft, 50, "True vs reconstructed angle", kFigTitPx, 12, cInk);
    text(kFramesRightPx, 50, subtitle, kSubPx, 32, cNote);

    // legend row: true | reconstructed | +-1 sigma | +-2 sigma  (widths measured)
    {
        const double y = 122;
        double x = kPxLeft;
        newBox(ndcX(x), ndcY(y + 15), ndcX(x + 60), ndcY(y - 15), cTrueFill, cTrueLine);
        text(x + 74, y, "True", kLegendPx, 12, cInk);                x += 74 + 81 + 48;
        newLine(ndcX(x), ndcY(y), ndcX(x + 60), ndcY(y), cReco, kStyleReco, 4);
        text(x + 74, y, "Reconstructed", kLegendPx, 12, cInk);       x += 74 + 258 + 48;
        newBox(ndcX(x), ndcY(y + 15), ndcX(x + 60), ndcY(y - 15), cBand1);
        text(x + 74, y, "#pm1#sigma", kLegendPx, 12, cInk);          x += 74 + 69 + 48;
        newBox(ndcX(x), ndcY(y + 15), ndcX(x + 60), ndcY(y - 15), cBand2);
        text(x + 74, y, "#pm2#sigma", kLegendPx, 12, cInk);
    }

    cv->Modified();
    cv->Update();
    cv->SaveAs((outbase + ".pdf").c_str());
    cv->SaveAs((outbase + ".png").c_str());
    delete cv;
    gROOT->SetBatch(wasBatch);
}

} // namespace recofig

void acceptance_mc()
{
    // ── configuration ──────────────────────────────────────────────────────
    const Long64_t nevents = 50000000;      // 5e7: omega = counts/nevents

    const double target_a         = 7.8*TMath::Sqrt(2)/2.0;  // 5.5154 cm semi-axis x
    const double target_b         = 7.8/2.0;                 // 3.9    cm semi-axis y
    const double target_thickness = 0.41105294e-4;           // cm, U target

    const double ppac_half = 10.0;          // half active area [cm]

    // internal PPAC layout (CreatePPAC), along the PPAC normal
    const double gas_gap         = 0.32;
    const double mylar_thickness = 1.7e-4;
    const double cath_off        = gas_gap + mylar_thickness;   // 0.32017 cm

    // true  -> cathode centres displaced along the BEAM direction
    // false -> cathode centres displaced along the PPAC NORMAL (CreatePPAC)
    const bool cathodes_on_beam_axis = false;

    // optional strip digitisation of the measured coordinates [cm]; 0 = off.
    // Only affects the reconstructed direction, never the acceptance.
    const double strip_pitch = 0.0;

    const double inv_sqrt2 = 1.0/TMath::Sqrt(2.0);

    const double d_PPACs = 5.0*TMath::Sqrt(2.0);   // lab-z spacing between PPACs
    const double step    = d_PPACs*inv_sqrt2;      // 5.0 cm in the detector frame

    // PPAC centres sit ON the beam axis (x = -z), target midway between them
    const double z_front = +0.5*step,  x_front = -z_front;   // (-2.5, +2.5)
    const double z_back  = -0.5*step,  x_back  = -z_back;    // (+2.5, -2.5)

    // Cathode planes: SAME ordering in both PPACs (back is not flipped)
    const double zX_front = z_front - cath_off;   // +2.17983
    const double zY_front = z_front + cath_off;   // +2.82017
    const double zX_back  = z_back  - cath_off;   // -2.82017
    const double zY_back  = z_back  + cath_off;   // -2.17983

    // Active-area centres in x
    const double xY_front = cathodes_on_beam_axis ? -zY_front : x_front;
    const double xX_front = cathodes_on_beam_axis ? -zX_front : x_front;
    const double xY_back  = cathodes_on_beam_axis ? -zY_back  : x_back;
    const double xX_back  = cathodes_on_beam_axis ? -zX_back  : x_back;

    // Per-coordinate lever arms.  Equal to `step` for the nominal layout, but
    // computed rather than assumed so a mirrored or re-aligned geometry stays
    // correctly reconstructed.
    const double Lx = zX_front - zX_back;
    const double Ly = zY_front - zY_back;

    std::cout << "detector-frame geometry\n"
              << "  front PPAC centre  (x,z) = (" << x_front << ", " << z_front << ")\n"
              << "  back  PPAC centre  (x,z) = (" << x_back  << ", " << z_back  << ")\n"
              << "  front  Y / X planes  z   = " << zY_front << " / " << zX_front << "\n"
              << "  back   Y / X planes  z   = " << zY_back  << " / " << zX_back  << "\n"
              << "  lever arms  Lx / Ly      = " << Lx << " / " << Ly << "\n"
              << "  cathode centres on beam axis: "
              << (cathodes_on_beam_axis ? "yes" : "no") << "\n"
              << "  strip pitch              = " << strip_pitch << " cm\n\n";

    const int    nbins_beam = 100;          // matches cos_theta_center in the CSV
    const int    nbins_det  = 100;           // matches cos_theta_det_center
    const double dcos_beam  = 1.0/nbins_beam;
    const double dcos_det   = 1.0/nbins_det;

    // ── bookkeeping ────────────────────────────────────────────────────────
    std::vector<std::vector<double>> cell_counts(
        nbins_beam, std::vector<double>(nbins_det, 0.0));

    Long64_t counts_forward = 0, counts_backward = 0;
    Long64_t coincidence    = 0, missed          = 0;
    Long64_t coinc_midplane = 0;   // would pass a single-plane-per-PPAC test
    Long64_t lost_2cathode  = 0;   // passes mid-plane but fails both cathodes

    TRandom3 rng(0);

    // Acceptance numerator and denominator are both binned in the TRUE angles,
    // so the ratio is a genuine acceptance and does not inherit any
    // reconstruction bias.
    TH2D* hist_thetas  = new TH2D("theta_det_beam", ";|cos#theta_{beam}|;cos#theta_{det}",
                                  nbins_beam, 0, 1, nbins_det, 0, 1);
    TH2D* hist_emitted = new TH2D("emitted",        ";|cos#theta_{beam}|;cos#theta_{det}",
                                  nbins_beam, 0, 1, nbins_det, 0, 1);
    hist_thetas->Sumw2();
    hist_emitted->Sumw2();

    // cos_theta_beam is genuinely signed (negative whenever theta_det > 45 deg),
    // so the axis must cover [-1,1] or those events go silently to underflow.
    TH2D* eff_beam = new TH2D("eff_beam", ";cos#theta;#phi",
                              200, -1, 1, 100, -TMath::Pi(), TMath::Pi());
    TH2D* eff_det  = new TH2D("eff_det",  ";cos#theta_{det};#phi_{det} (deg)",
                              100, 0, 1, 100, -180, 180);

    TH1D* hist_cos_theta_det     = new TH1D("cos_theta_det",     "", nbins_det,  0, 1);
    TH1D* hist_cos_theta         = new TH1D("cos_theta",         "", nbins_beam, 0, 1);
    TH1D* hist_cos_theta_emitted = new TH1D("cos_theta_emitted", "", nbins_beam, 0, 1);
    hist_cos_theta->Sumw2();
    hist_cos_theta_emitted->Sumw2();

    TH1D* hist_x_front = new TH1D("x_front", ";x_{front} (cm)", 200, -20, 20);
    TH1D* hist_x_back  = new TH1D("x_back",  ";x_{back} (cm)",  200, -20, 20);
    TH1D* hist_y_back  = new TH1D("y_back",  ";y_{back} (cm)",  200, -20, 20);

    // ── true vs reconstructed angle check ──────────────────────────────────
    // Signed cos is used throughout: folding with Abs() before ACos() makes
    // true and reco land on opposite sides of 90 deg and fakes a residual.
    const int    nbins_theta   = 1000;
    const double dtheta_range  = 5.0;    // deg; widen if the 2D map overflows

    TH1D* hist_theta_true_beam = new TH1D("theta_true_beam",
        ";#theta_{beam} (deg);counts", nbins_theta, 0, 180);
    TH1D* hist_theta_reco_beam = new TH1D("theta_reco_beam",
        ";#theta_{beam} (deg);counts", nbins_theta, 0, 180);
    TH1D* hist_theta_true_det  = new TH1D("theta_true_det",
        ";#theta_{det} (deg);counts",  nbins_theta, 0, 90);
    TH1D* hist_theta_reco_det  = new TH1D("theta_reco_det",
        ";#theta_{det} (deg);counts",  nbins_theta, 0, 90);

    TH2D* hist_dtheta_vs_theta_beam = new TH2D("dtheta_vs_theta_beam",
        ";#theta_{true,beam} (deg);#theta_{reco}-#theta_{true} (deg)",
        nbins_theta, 0, 180, 120, -dtheta_range, dtheta_range);
    TH2D* hist_dtheta_vs_theta_det = new TH2D("dtheta_vs_theta_det",
        ";#theta_{true,det} (deg);#theta_{reco}-#theta_{true} (deg)",
        nbins_theta, 0, 90, 120, -dtheta_range, dtheta_range);

    // same check, but in cos(theta) instead of theta (deg)
    const int    nbins_costheta  = 200;
    const double dcostheta_range = 0.02;   // widen if the 2D map overflows

    TH1D* hist_costheta_true_beam = new TH1D("costheta_true_beam",
        ";cos#theta_{beam};counts", nbins_costheta, -1, 1);
    TH1D* hist_costheta_reco_beam = new TH1D("costheta_reco_beam",
        ";cos#theta_{beam};counts", nbins_costheta, -1, 1);
    TH1D* hist_costheta_true_det  = new TH1D("costheta_true_det",
        ";cos#theta_{det};counts",  nbins_costheta, 0, 1);
    TH1D* hist_costheta_reco_det  = new TH1D("costheta_reco_det",
        ";cos#theta_{det};counts",  nbins_costheta, 0, 1);

    TH2D* hist_dcostheta_vs_costheta_beam = new TH2D("dcostheta_vs_costheta_beam",
        ";cos#theta_{true,beam};cos#theta_{reco}-cos#theta_{true}",
        nbins_costheta, -1, 1, 120, -dcostheta_range, dcostheta_range);
    TH2D* hist_dcostheta_vs_costheta_det = new TH2D("dcostheta_vs_costheta_det",
        ";cos#theta_{true,det};cos#theta_{reco}-cos#theta_{true}",
        nbins_costheta, 0, 1, 120, -dcostheta_range, dcostheta_range);

    // ── display histograms for the figures + migration counters ──────────────
    // Binning of the figures: 1 deg in theta, 0.02 in cos(theta).  n_out / n_in
    // count the events whose true and reco values fall in DIFFERENT bins; they
    // give the exact uncertainty of N_reco - N_true = n_in - n_out (see recofig).
    struct Disp { TH1D *t, *r, *in, *out; };
    auto makeDisp = [](const char* n, int nb, double lo, double hi) {
        Disp d;
        d.t   = new TH1D(Form("disp_%s_true", n), "", nb, lo, hi);
        d.r   = new TH1D(Form("disp_%s_reco", n), "", nb, lo, hi);
        d.in  = new TH1D(Form("disp_%s_in",   n), "", nb, lo, hi);
        d.out = new TH1D(Form("disp_%s_out",  n), "", nb, lo, hi);
        for (TH1D* h : {d.t, d.r, d.in, d.out}) h->SetDirectory(nullptr);
        return d;
    };
    auto fillDisp = [](Disp& d, double vTrue, double vReco) {
        d.t->Fill(vTrue);
        d.r->Fill(vReco);
        if (d.t->FindFixBin(vTrue) != d.r->FindFixBin(vReco)) {
            d.out->Fill(vTrue);
            d.in ->Fill(vReco);
        }
    };
    Disp dThetaBeam = makeDisp("theta_beam", 180,  0., 180.);
    Disp dThetaDet  = makeDisp("theta_det",   90,  0.,  90.);
    Disp dCosBeam   = makeDisp("cos_beam",   100, -1.,   1.);
    Disp dCosDet    = makeDisp("cos_det",     50,  0.,   1.);

    // ── event loop ─────────────────────────────────────────────────────────
    for (Long64_t i = 0; i < nevents; ++i) {

        if (i % 5000000 == 0) std::cout << "  event " << i << std::endl;

        // uniform over the elliptical footprint: uniform in the unit disk,
        // then scale.  In the detector frame the foil is flat at z_det = 0.
        const double u_r   = TMath::Sqrt(rng.Uniform(0., 1.));
        const double u_phi = TMath::TwoPi() * rng.Uniform(0., 1.);
        const double ox = target_a * u_r * TMath::Cos(u_phi);
        const double oy = target_b * u_r * TMath::Sin(u_phi);
        const double oz = (0.5 - rng.Uniform(0., 1.)) * target_thickness;

        // isotropic into the forward hemisphere; the partner is exactly
        // back-to-back, so one event covers the full sphere as a pair
        const double cth = rng.Uniform(0., 1.);
        const double sth = TMath::Sqrt(1. - cth*cth);
        const double phd = TMath::TwoPi() * rng.Uniform(0., 1.);

        const double dfx = TMath::Cos(phd)*sth;
        const double dfy = TMath::Sin(phd)*sth;
        const double dfz = cth;

        if (dfz <= 0.) continue;            // parallel to the planes (guard)

        // TRUE emission direction in the beam frame.  (dfx,dfy,dfz) is already
        // a unit vector in the detector frame, so no renormalisation is needed.
        const double cos_theta_true_beam = (-dfx + dfz) * inv_sqrt2;
        // cos_theta_true_det is just dfz == cth, already at hand.

        // intersection of a ray from (ox,oy,oz) with the plane z = zp
        auto cross = [&](double zp, double vx, double vy, double vz,
                         double& X, double& Y) -> bool {
            if (vz == 0.) return false;
            const double t = (zp - oz) / vz;
            if (t < 0.) return false;
            X = ox + t*vx;
            Y = oy + t*vy;
            return true;
        };
        auto inArea = [&](double X, double Y, double xoff) -> bool {
            return TMath::Abs(X - xoff) <= ppac_half
                && TMath::Abs(Y)        <= ppac_half;
        };

        // ---- forward fragment through the two front cathodes --------------
        double xfX, yfX, xfY, yfY;
        const bool okfX = cross(zX_front,  dfx,  dfy,  dfz, xfX, yfX);
        const bool okfY = cross(zY_front,  dfx,  dfy,  dfz, xfY, yfY);
        const bool hit_front = okfX && okfY
                            && inArea(xfX, yfX, xX_front)
                            && inArea(xfY, yfY, xY_front);

        // ---- backward fragment through the two back cathodes --------------
        double xbX, ybX, xbY, ybY;
        const bool okbX = cross(zX_back,  -dfx, -dfy, -dfz, xbX, ybX);
        const bool okbY = cross(zY_back,  -dfx, -dfy, -dfz, xbY, ybY);
        const bool hit_back = okbX && okbY
                           && inArea(xbX, ybX, xX_back)
                           && inArea(xbY, ybY, xY_back);

        if (!okfX || !okfY || !okbX || !okbY) continue;

        // measured point per PPAC: each coordinate from its OWN cathode plane.
        // Plane identity, not arrival order: X always at z_ppac + c, Y always
        // at z_ppac - c, in both PPACs.
        const double xf = xfX, yf = (yfY);
        const double xb = xbX, yb = (ybY);

        // diagnostic: single-plane-per-PPAC test, for comparison
        double xfm, yfm, xbm, ybm;
        cross(z_front,  dfx,  dfy,  dfz, xfm, yfm);
        cross(z_back,  -dfx, -dfy, -dfz, xbm, ybm);
        const bool hit_mid = inArea(xfm, yfm, x_front) && inArea(xbm, ybm, x_back);

        // Direction from per-coordinate slopes.  Both fragments share a vertex
        // and are exactly back-to-back, so the vertex cancels and this
        // reproduces the emission direction exactly (machine precision) when
        // strip_pitch = 0 -- the internal closure check.
        const double tx = (xf - xb) / Lx;
        const double ty = (yf - yb) / Ly;
        const double nn = TMath::Sqrt(tx*tx + ty*ty + 1.0);

        const double cos_theta_det = 1.0 / nn;
        const double phi_det       = TMath::ATan2(ty, tx);

        // back to the beam frame: v_lab = R_y(45) v_local
        //   => cos_theta_beam = (-vx + vz)/sqrt2
        const double vx = tx / nn;
        const double vy = ty / nn;
        const double vz = 1.0 / nn;

        const double nx = ( vx + vz) * inv_sqrt2;
        const double ny =   vy;
        const double nz = (-vx + vz) * inv_sqrt2;

        const double cos_theta = nz;                     // already normalised
        const double phi       = TMath::ATan2(ny, nx);

        // acceptance denominator, in TRUE angles
        hist_emitted->Fill(TMath::Abs(cos_theta_true_beam), cth);
        hist_cos_theta_emitted->Fill(TMath::Abs(cos_theta_true_beam));

        if (hit_front) ++counts_forward;
        if (hit_back)  ++counts_backward;
        if (!hit_front || !hit_back) ++missed;
        if (hit_mid)   ++coinc_midplane;
        if (hit_mid && !(hit_front && hit_back)) ++lost_2cathode;

        if (hit_front && hit_back) {
            ++coincidence;

            // clamp rather than skip: cos = 1 exactly is reachable, and a
            // `continue` here would also drop the event from every histogram
            // filled further down.
            const int bin_beam =
                int(TMath::Abs(cos_theta_true_beam) / dcos_beam);
            const int bin_det  =
                int(cth / dcos_det);

            cell_counts[bin_beam][bin_det] += 1.;

            hist_thetas->Fill(TMath::Abs(cos_theta_true_beam), cth);
            hist_cos_theta->Fill(TMath::Abs(cos_theta_true_beam));
            hist_cos_theta_det->Fill(cth);
            hist_x_front->Fill(xf-x_front);
            hist_x_back->Fill(xb-x_back);
            hist_y_back->Fill(yb);
            eff_beam->Fill(cos_theta, phi);
            eff_det->Fill(cos_theta_det, phi_det*TMath::RadToDeg());

            // ── true vs reconstructed angle, only where we have a full
            //    coincidence (i.e. an actual reconstructed point pair) ──────
            const double theta_true_beam_deg =
                TMath::ACos(cos_theta_true_beam) * TMath::RadToDeg();
            const double theta_reco_beam_deg =
                TMath::ACos(cos_theta)           * TMath::RadToDeg();
            const double theta_true_det_deg  = TMath::ACos(cth)           * TMath::RadToDeg();
            const double theta_reco_det_deg  = TMath::ACos(cos_theta_det) * TMath::RadToDeg();

            hist_theta_true_beam->Fill(theta_true_beam_deg);
            hist_theta_reco_beam->Fill(theta_reco_beam_deg);
            hist_theta_true_det->Fill(theta_true_det_deg);
            hist_theta_reco_det->Fill(theta_reco_det_deg);

            hist_dtheta_vs_theta_beam->Fill(theta_true_beam_deg,
                                             theta_reco_beam_deg - theta_true_beam_deg);
            hist_dtheta_vs_theta_det->Fill(theta_true_det_deg,
                                            theta_reco_det_deg  - theta_true_det_deg);

            // ── same check in cos(theta) instead of theta (deg) ─────────────
            hist_costheta_true_beam->Fill(cos_theta_true_beam);
            hist_costheta_reco_beam->Fill(cos_theta);
            hist_costheta_true_det->Fill(cth);
            hist_costheta_reco_det->Fill(cos_theta_det);

            hist_dcostheta_vs_costheta_beam->Fill(cos_theta_true_beam,
                                                   cos_theta - cos_theta_true_beam);
            hist_dcostheta_vs_costheta_det->Fill(cth,
                                                  cos_theta_det - cth);

            // ── figures: distributions and migrations ───────────────────────
            fillDisp(dThetaBeam, theta_true_beam_deg, theta_reco_beam_deg);
            fillDisp(dThetaDet,  theta_true_det_deg,  theta_reco_det_deg);
            fillDisp(dCosBeam,   cos_theta_true_beam, cos_theta);
            fillDisp(dCosDet,    cth,                 cos_theta_det);
        }
    }

    // ── summary ────────────────────────────────────────────────────────────
    std::cout << "\n---------- RESULTS ----------\n"
              << "  events              = " << nevents << "\n"
              << "  hits front (2 cath) = " << 100.*counts_forward /nevents << " %\n"
              << "  hits back  (2 cath) = " << 100.*counts_backward/nevents << " %\n"
              << "  coincidences        = " << 100.*coincidence    /nevents << " %\n"
              << "  missed              = " << 100.*missed         /nevents << " %\n"
              << "\n  single-plane test   = " << 100.*coinc_midplane/nevents << " %\n"
              << "  lost by requiring both cathodes = "
              << 100.*lost_2cathode/nevents << " %  ("
              << (coinc_midplane > 0
                  ? 100.*lost_2cathode/double(coinc_midplane) : 0.)
              << " % of single-plane coincidences)\n"
              << "\n  <theta_reco - theta_true> (beam) = "
              << hist_dtheta_vs_theta_beam->GetMean(2) << " deg,  RMS = "
              << hist_dtheta_vs_theta_beam->GetRMS(2)  << " deg\n"
              << "  <theta_reco - theta_true> (det)  = "
              << hist_dtheta_vs_theta_det->GetMean(2)  << " deg,  RMS = "
              << hist_dtheta_vs_theta_det->GetRMS(2)   << " deg\n"
              << "\n  <costheta_reco - costheta_true> (beam) = "
              << hist_dcostheta_vs_costheta_beam->GetMean(2) << ",  RMS = "
              << hist_dcostheta_vs_costheta_beam->GetRMS(2)  << "\n"
              << "  <costheta_reco - costheta_true> (det)  = "
              << hist_dcostheta_vs_costheta_det->GetMean(2)  << ",  RMS = "
              << hist_dcostheta_vs_costheta_det->GetRMS(2)   << "\n";

    if (strip_pitch <= 0.)
        std::cout << "\n  (strip_pitch = 0: the residuals above are the closure\n"
                  << "   test and must be zero to machine precision)\n";

    std::cout << "\n  note: `omega` in the CSV is a probability per emitted pair\n"
              << "  from forward-hemisphere sampling.  Multiply by 4*pi to get\n"
              << "  a solid angle in sr.\n" << std::endl;

    // ── CSV, same format as acceptance_coincidence.csv ─────────────────────
    std::ofstream csv("/Users/nico/Desktop/Tese/Analysis/acceptance_coincidence.csv");
    csv << "cos_theta_center,cos_theta_det_center,counts,omega\n";
    csv << std::fixed;
    for (int j = 0; j < nbins_beam; ++j) {
        const double c_beam = (j + 0.5) * dcos_beam;
        for (int k = 0; k < nbins_det; ++k) {
            const double c_det  = (k + 0.5) * dcos_det;
            const double counts = cell_counts[j][k];
            csv << std::setprecision(10) << c_beam << ','
                << std::setprecision(10) << c_det  << ','
                << std::setprecision(0)  << counts << ','
                << std::setprecision(10) << counts/double(nevents) << '\n';
        }
    }
    csv.close();
    std::cout << "wrote acceptance_coincidence.csv ("
              << nbins_beam*nbins_det << " rows)\n";

    // ── ROOT output ────────────────────────────────────────────────────────
    TFile* fout = new TFile("mc_acceptance.root", "RECREATE");

    hist_emitted->Write();
    eff_beam->Write();
    eff_det->Write();
    hist_x_front->Write();
    hist_x_back->Write();
    hist_y_back->Write();
    hist_cos_theta_det->Write();

    TH1D* acceptance = new TH1D("acceptance", ";cos#theta_{det};acceptance",
                                nbins_det, 0, 1);
    for (int k = 0; k < nbins_det; ++k) {
        double emitted = 0., hits = 0.;
        for (int j = 0; j < nbins_beam; ++j) {
            emitted += hist_emitted->GetBinContent(j+1, k+1);
            hits    += cell_counts[j][k];
        }
        const double a = (emitted > 0.) ? hits/emitted : 0.;
        acceptance->SetBinContent(k+1, a);
        acceptance->SetBinError  (k+1, (emitted > 0.)
                                       ? TMath::Sqrt(a*(1.-a)/emitted) : 0.);
    }
    acceptance->Write();

    // clone before dividing, so the raw coincidence map survives.  The "B"
    // option gives binomial errors, which is what a subset/total ratio needs;
    // the default assumes independent numerator and denominator.
    TH2D* hist_ratio = (TH2D*)hist_thetas->Clone("acceptance_2d");
    hist_ratio->Divide(hist_thetas, hist_emitted, 1., 1., "B");
    hist_ratio->Write();
    hist_thetas->Write();

    TH1D* h_cos_ratio = (TH1D*)hist_cos_theta->Clone("cos_theta_acceptance");
    h_cos_ratio->Divide(hist_cos_theta, hist_cos_theta_emitted, 1., 1., "B");
    h_cos_ratio->Write();
    hist_cos_theta->Write();
    hist_cos_theta_emitted->Write();

    // ── true-vs-reconstructed angle histograms ─────────────────────────────
    hist_theta_true_beam->Write();
    hist_theta_reco_beam->Write();
    hist_theta_true_det->Write();
    hist_theta_reco_det->Write();
    hist_dtheta_vs_theta_beam->Write();
    hist_dtheta_vs_theta_det->Write();

    hist_costheta_true_beam->Write();
    hist_costheta_reco_beam->Write();
    hist_costheta_true_det->Write();
    hist_costheta_reco_det->Write();
    hist_dcostheta_vs_costheta_beam->Write();
    hist_dcostheta_vs_costheta_det->Write();

    // migrations, so the pulls can be recomputed from the file
    for (Disp* d : {&dThetaBeam, &dThetaDet, &dCosBeam, &dCosDet})
        for (TH1D* h : {d->t, d->r, d->in, d->out}) h->Write();

    // ── figures: true vs reconstructed, with pulls, in theta and cos(theta) ─
    {
        using recofig::AxisSpec;
        const std::string sub = (strip_pitch > 0.)
            ? std::string(Form("strip pitch %g cm", strip_pitch))
            : std::string("strip pitch 0: closure test");

        // residual maps at the display binning (x rebinned; y unchanged)
        TH2D* mThB = (TH2D*)hist_dtheta_vs_theta_beam      ->RebinX(5, "map_theta_beam");
        TH2D* mThD = (TH2D*)hist_dtheta_vs_theta_det       ->RebinX(5, "map_theta_det");
        TH2D* mCB  = (TH2D*)hist_dcostheta_vs_costheta_beam->RebinX(2, "map_cos_beam");
        TH2D* mCD  = (TH2D*)hist_dcostheta_vs_costheta_det ->RebinX(2, "map_cos_det");
        for (TH2D* m : {mThB, mThD, mCB, mCD}) m->SetDirectory(nullptr);

        const double sTh = recofig::niceStep(dtheta_range / 2.5);
        const double sC  = recofig::niceStep(dcostheta_range / 2.5);
        const AxisSpec YmapTheta = {-dtheta_range,    dtheta_range,    0., sTh, recofig::minorsFor(sTh), false};
        const AxisSpec YmapCos   = {-dcostheta_range, dcostheta_range, 0., sC,  recofig::minorsFor(sC),  false};

        const recofig::Column colsTheta[2] = {
            {dThetaBeam.t, dThetaBeam.r, dThetaBeam.in, dThetaBeam.out, mThB,
             {0., 180., 0., 30., 3, true}, YmapTheta, "Beam frame",     "#theta_{beam} (#circ)"},
            {dThetaDet.t,  dThetaDet.r,  dThetaDet.in,  dThetaDet.out,  mThD,
             {0.,  90., 0., 15., 3, true}, YmapTheta, "Detector frame", "#theta_{det} (#circ)"}};
        recofig::drawFigure(colsTheta, "1#circ", "#theta_{reco} #minus #theta_{true} (#circ)",
                            sub, "angle_true_vs_reco");

        const recofig::Column colsCos[2] = {
            {dCosBeam.t, dCosBeam.r, dCosBeam.in, dCosBeam.out, mCB,
             {-1., 1., 0., 0.5, 5, true}, YmapCos, "Beam frame",     "cos#theta_{beam}"},
            {dCosDet.t,  dCosDet.r,  dCosDet.in,  dCosDet.out,  mCD,
             { 0., 1., 0., 0.2, 4, true}, YmapCos, "Detector frame", "cos#theta_{det}"}};
        recofig::drawFigure(colsCos, "0.02", "cos#theta_{reco} #minus cos#theta_{true}",
                            sub, "costheta_true_vs_reco");

        for (TH2D* m : {mThB, mThD, mCB, mCD}) delete m;
    }

    fout->Close();
    std::cout << "wrote mc_acceptance.root\n";
    std::cout << "wrote angle_true_vs_reco.png / .pdf\n";
    std::cout << "wrote costheta_true_vs_reco.png / .pdf\n";
}