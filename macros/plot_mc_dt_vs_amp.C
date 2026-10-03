// plot_mc_dt_vs_amp.C
//
// Simulación: tof1 - tof0 frente a amp0 + amp1.
//   mapa (kBird, escala log):  ppac0<6 && fOut_target_i!=0
//   contornos:                 ppac0<6 && fOut_target_i==0 && coinc_type==1
//                              && hasTrue_i && Z0>4 && Z1>4
//     línea gruesa: región que contiene el 68 % de estos eventos
//     línea fina:   región que contiene el 95 %
//     (densidad estimada con un kernel gaussiano, regla de Scott)
//
//   root -l
//   root [0] .x plot_mc_dt_vs_amp.C
//
// Guarda mc_dt_vs_amp.pdf y mc_dt_vs_amp.png.

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPaletteAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TList.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

namespace mcdt {

// ---- entrada, cortes y textos (editar aquí) -------------------------------------------
const char* kFile   = "/Users/nico/Desktop/Tese/Analysis/montecarlo/data/mc_complete.root";
const char* kTree   = "CoincTree";
const char* kX      = "amp0+amp1";
const char* kY      = "tof1-tof0";
const char* kBkgCut = "ppac0<6 && fOut_target_i!=0";
const char* kSelCut = "ppac0<6 && fOut_target_i==0 && coinc_type==1 && hasTrue_i && Z0>4 && Z1>4";

const char* kXTitle   = "#it{A}_{0} + #it{A}_{1}  (arb. units)";
const char* kYTitle   = "#it{t}_{1} #minus #it{t}_{0}  (ns)";
const char* kBkgLabel = "Outside the target";
const char* kSelLabel = "Target fission";

const double kLevel1 = 0.68;          // contorno interior (línea gruesa)
const double kLevel2 = 0.95;          // contorno exterior (línea fina)
const int    kLineColor = kBlack;     // color de los contornos

const int kNbX = 200, kNbY = 200;
// rango: si lo < hi se usa tal cual; si no, automático
const double kXlo = 0, kXhi = 0;
const double kYlo = 0, kYhi = 0;

// ---------------------------------------------------------------------------------------
void load(TTree* t, const char* cut, std::vector<double>& x, std::vector<double>& y)
{
    t->SetEstimate(t->GetEntries() + 1);
    const Long64_t n = std::max<Long64_t>(0, t->Draw(Form("%s:%s", kX, kY), cut, "goff"));
    x.assign(t->GetV1(), t->GetV1() + n);
    y.assign(t->GetV2(), t->GetV2() + n);
}

double quantile(std::vector<double> v, double q)
{
    const size_t k = (size_t)(q * (v.size() - 1));
    std::nth_element(v.begin(), v.begin() + k, v.end());
    return v[k];
}

// redondea [lo, hi] hacia fuera a múltiplos de un paso 1-2-5
void niceRange(double& lo, double& hi)
{
    const double raw = (hi - lo) / 5.0;
    const double e = std::pow(10.0, std::floor(std::log10(raw)));
    const double f = raw / e;
    const double step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * e;
    lo = std::floor(lo / step) * step;
    hi = std::ceil(hi / step) * step;
}

// densidad de los eventos del blanco: histograma fino suavizado con un kernel gaussiano
// (anchura por la regla de Scott en 2D, h = sigma * n^(-1/6))
TH2D* density(const std::vector<double>& x, const std::vector<double>& y,
              double xlo, double xhi, double ylo, double yhi)
{
    const int G = 200;
    TH2D* d = new TH2D("h_mcdt_dens", "", G, xlo, xhi, G, ylo, yhi);
    d->SetDirectory(nullptr);
    for (size_t i = 0; i < x.size(); ++i) d->Fill(x[i], y[i]);

    auto sigma = [](const std::vector<double>& v) {
        const double iqr = (quantile(v, 0.75) - quantile(v, 0.25)) / 1.349;
        double s = 0, ss = 0;
        for (double a : v) { s += a; ss += a * a; }
        const double sd = std::sqrt(std::max(0.0, ss / v.size() - (s / v.size()) * (s / v.size())));
        return iqr > 0 ? std::min(sd, iqr) : sd;
    };
    const double f = std::pow((double)x.size(), -1.0 / 6.0);
    const double sx = sigma(x) * f / d->GetXaxis()->GetBinWidth(1);   // en bins
    const double sy = sigma(y) * f / d->GetYaxis()->GetBinWidth(1);

    std::vector<double> a(G * G);
    for (int i = 0; i < G; ++i)
        for (int j = 0; j < G; ++j) a[i + G * j] = d->GetBinContent(i + 1, j + 1);
    auto pass = [&](double s, bool alongX) {
        const int R = std::max(1, (int)std::ceil(4 * s));
        std::vector<double> k(2 * R + 1), out(G * G, 0.0);
        for (int r = -R; r <= R; ++r) k[r + R] = std::exp(-0.5 * r * r / (s * s));
        for (int i = 0; i < G; ++i)
            for (int j = 0; j < G; ++j) {
                double v = 0, w = 0;
                for (int r = -R; r <= R; ++r) {
                    const int ii = alongX ? i + r : i, jj = alongX ? j : j + r;
                    if (ii < 0 || ii >= G || jj < 0 || jj >= G) continue;
                    v += k[r + R] * a[ii + G * jj];
                    w += k[r + R];
                }
                out[i + G * j] = v / w;
            }
        a.swap(out);
    };
    pass(sx, true);
    pass(sy, false);
    for (int i = 0; i < G; ++i)
        for (int j = 0; j < G; ++j) d->SetBinContent(i + 1, j + 1, a[i + G * j]);
    return d;
}

// nivel de densidad cuyo contorno encierra la fracción 'frac' de los eventos
double levelFor(TH2D* d, const std::vector<double>& x, const std::vector<double>& y, double frac)
{
    std::vector<double> v(x.size());
    for (size_t i = 0; i < x.size(); ++i) v[i] = d->Interpolate(x[i], y[i]);
    return quantile(v, 1.0 - frac);
}

// líneas de nivel del histograma 'd' (una lista de TGraph por nivel)
std::vector<std::vector<TGraph*>> contours(TH2D* d, double lvLow, double lvHigh)
{
    std::vector<std::vector<TGraph*>> out(2);
    double lv[2] = {lvLow, lvHigh};
    TH2D* h = (TH2D*)d->Clone("h_mcdt_cont");
    h->SetDirectory(nullptr);
    h->SetContour(2, lv);
    TCanvas tmp("c_mcdt_tmp", "", 400, 400);
    h->Draw("CONT LIST");
    tmp.Update();
    if (TObjArray* arr = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours"))
        for (int l = 0; l < 2 && l < arr->GetEntries(); ++l)
            if (TList* lst = (TList*)arr->At(l))
                for (TObject* o : *lst) out[l].push_back(new TGraph(*(TGraph*)o));
    delete h;
    return out;
}

}  // namespace mcdt

// ---------------------------------------------------------------------------------------
void plot_mc_dt_vs_amp()
{
    using namespace mcdt;

    TFile* f = TFile::Open(kFile, "READ");
    if (!f || f->IsZombie()) { printf("[ERROR] no se puede abrir %s\n", kFile); return; }
    TTree* t = (TTree*)f->Get(kTree);
    if (!t) { printf("[ERROR] no está el árbol %s\n", kTree); f->Close(); return; }

    std::vector<double> bx, by, sx, sy;
    load(t, kBkgCut, bx, by);
    load(t, kSelCut, sx, sy);
    f->Close();
    printf("Fondo (mapa):           N = %zu\nBlanco (contornos):     N = %zu\n", bx.size(), sx.size());
    if (bx.empty() || sx.size() < 20) { printf("[ERROR] no hay eventos suficientes\n"); return; }

    // ---- rango: 99.9 % del fondo y todos los eventos del blanco, redondeado ----
    double xlo = kXlo, xhi = kXhi, ylo = kYlo, yhi = kYhi;
    if (!(xlo < xhi)) {
        xlo = std::min(quantile(bx, 0.0005), *std::min_element(sx.begin(), sx.end()));
        xhi = std::max(quantile(bx, 0.9995), *std::max_element(sx.begin(), sx.end()));
        const bool pos = xlo >= 0;
        niceRange(xlo, xhi);
        if (pos && xlo < 0) xlo = 0;
    }
    if (!(ylo < yhi)) {
        ylo = std::min(quantile(by, 0.0005), *std::min_element(sy.begin(), sy.end()));
        yhi = std::max(quantile(by, 0.9995), *std::max_element(sy.begin(), sy.end()));
        const double m = std::max(std::fabs(ylo), std::fabs(yhi));   // simétrico en 0
        ylo = -m; yhi = m;
        niceRange(ylo, yhi);
    }

    if (TObject* o = gROOT->FindObject("h_mcdt")) delete o;
    TH2D* h = new TH2D("h_mcdt", "", kNbX, xlo, xhi, kNbY, ylo, yhi);
    h->SetDirectory(nullptr);
    for (size_t i = 0; i < bx.size(); ++i) h->Fill(bx[i], by[i]);
    std::vector<double>().swap(bx);
    std::vector<double>().swap(by);

    // ---- contornos de los eventos del blanco ----
    TH2D* dens = density(sx, sy, xlo, xhi, ylo, yhi);
    const double lv1 = levelFor(dens, sx, sy, kLevel1);
    const double lv2 = levelFor(dens, sx, sy, kLevel2);
    const Bool_t wasBatch = gROOT->IsBatch();
    gROOT->SetBatch(kTRUE);
    std::vector<std::vector<TGraph*>> iso = contours(dens, lv2, lv1);   // [0] = 95 %, [1] = 68 %

    // ---- estilo ----
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(255);

    if (TObject* o = gROOT->FindObject("c_mcdt")) delete o;
    TCanvas* c = new TCanvas("c_mcdt", "", 800, 700);
    c->SetCanvasSize(800, 700);
    c->SetFillColor(kWhite);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.17);
    c->SetBottomMargin(0.13);
    c->SetTopMargin(0.10);                    // hueco arriba para la leyenda
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetLogz(1);

    h->GetXaxis()->SetTitle(kXTitle);
    h->GetYaxis()->SetTitle(kYTitle);
    h->GetZaxis()->SetTitle("Counts");
    for (TAxis* a : {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()}) {
        a->SetTitleFont(43); a->SetTitleSize(26);
        a->SetLabelFont(43); a->SetLabelSize(22);
    }
    h->GetXaxis()->SetTitleOffset(1.25);
    h->GetYaxis()->SetTitleOffset(1.45);
    h->GetZaxis()->SetTitleOffset(1.25);
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetNdivisions(505);
    h->SetMinimum(0.9);
    h->Draw("COLZ");
    c->Update();
    if (TPaletteAxis* pal = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette")) {
        pal->SetX1NDC(0.845);
        pal->SetX2NDC(0.875);
        pal->SetY1NDC(0.13);
        pal->SetY2NDC(0.90);
    }

    for (int l = 0; l < 2; ++l)
        for (TGraph* g : iso[l]) {
            g->SetLineColor(kLineColor);
            g->SetLineWidth(l == 0 ? 2 : 3);      // 95 % fina, 68 % gruesa
            g->Draw("L");
        }

    // ---- leyenda (encima del marco, así nunca tapa los datos) ----
    TH1D* box = new TH1D("h_mcdt_box", "", 1, 0, 1);
    box->SetDirectory(nullptr);
    box->SetFillColor(gStyle->GetColorPalette(170));
    box->SetLineColor(gStyle->GetColorPalette(170));
    TGraph* ll = new TGraph(1);
    ll->SetLineColor(kLineColor);
    ll->SetLineWidth(3);

    TLegend* leg = new TLegend(0.14, 0.91, 0.83, 0.985);
    leg->SetNColumns(2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(43);
    leg->SetTextSize(22);
    leg->AddEntry(box, kBkgLabel, "f");
    leg->AddEntry(ll, Form("%s  (%.0f %%, %.0f %%)", kSelLabel, 100 * kLevel1, 100 * kLevel2), "l");
    leg->Draw();

    c->RedrawAxis();
    c->Modified();
    c->Update();
    c->Print("mc_dt_vs_amp.pdf", "pdf Portrait");   // "Portrait": página sin girar
    c->Print("mc_dt_vs_amp.png");

    gROOT->SetBatch(wasBatch);
    if (!wasBatch) c->DrawClone();
}