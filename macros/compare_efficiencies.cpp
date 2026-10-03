// plot_ce_data_vs_mc.C
//
// tof_1 - tof_0 del cerio: datos y simulación en un lienzo 2 x 1 (uno al lado
// del otro), mismo binning (60 bins entre -10 y 10 ns), mismo formato:
//   (a) datos       granate   (cos theta' > 0.8, 400 < E_n < 2000 MeV)
//   (b) simulación  azul
//
// Datos: cos theta' NO se toma de la rama cos_theta_det; se recalcula evento a
// evento a partir de x0, y0, x1, y1 con la corrección de offsets (medias de las
// posiciones), igual que compute_angles() / fillHistograms() del análisis.
// MC: se usa la rama cos_theta_det de la simulación.
//
//   root -l
//   root [0] .x plot_ce_data_vs_mc.C
//
// Guarda Ce_tof_data_vs_mc.pdf y .png. En la terminal se imprimen los offsets,
// la media y el rms de cada una (en su 90 % central) y el desplazamiento datos - MC.

#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TROOT.h"
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

// ---- entrada y cortes (editar aquí) ------------------------------------------------
// Datos: cortes SIN el de ángulo; cos theta' (recalculado) > kDataCosCut va aparte
static const char* kDataFile = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences_final.root";
static const char* kDataTree = "events_cerium";
static const char* kDataCut  =
    "neutron_energy>400 && neutron_energy<2000 && full_position && det0<6"
    " && abs(amp1-amp0)/(amp0+amp1)<0.5 && amp0+amp1>18e3";
static const char* kDataDt   = "tof1-tof0";
static const double kDataCosCut = 0.8;
// runs excluidos en el análisis (solo se aplica si el árbol tiene RunNumber)
static const char* kDataRunCut = "RunNumber!=118771 && RunNumber!=118789";

static const char* kMcFile   = "/Users/nico/Desktop/Tese/Analysis/montecarlo/data/mc_complete.root";
static const char* kMcTree   = "CoincTree";
static const char* kMcCut    = "cos_theta_det>0.8 && fOut_target_i==0 && ppac0<6 && hasTrue_i && coinc_type==1 && neutronE>400";
static const char* kMcDt     = "tof1-tof0";        // cambia si en la MC las ramas se llaman distinto

static const int    kNb = 60;
static const double kLo = -10., kHi = 10.;

static const int kFont = 43, kFontBold = 63;

// ------------------------------------------------------------------------------------
// ángulo del detector, igual que compute_angles() del análisis (solo cos_theta_det;
// los offsets en x no entran en la fórmula)
struct Offsets { double x0 = 0, x1 = 0, y0 = 0, y1 = 0; Long64_t n = 0; };

static double cosThetaDet(double x0, double y0, double x1, double y1, const Offsets& o)
{
    const double dx = (2*x1 + 2.5) - (2*x0 - 2.5);
    const double dy =  2*y1 - 2*y0 - (2*o.y1 - 2*o.y0);
    const double dz =  5.0;
    const double nd = std::sqrt(dx*dx + dy*dy + dz*dz);
    return nd > 0.0 ? dz / nd : -999.;
}

// relleno directo con TTree::Draw (MC, ángulo de la rama)
static TH1D* fillFrom(const char* file, const char* treeName, const char* expr,
                      const char* cut, const char* hname)
{
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) { printf("[ERROR] no se puede abrir %s\n", file); return nullptr; }
    TTree* t = (TTree*)f->Get(treeName);
    if (!t) { printf("[ERROR] no está el árbol %s en %s\n", treeName, file); f->Close(); return nullptr; }
    if (TObject* o = gROOT->FindObject(hname)) delete o;
    f->cd();
    TH1D* h = new TH1D(hname, "", kNb, kLo, kHi);
    h->Sumw2();
    t->Draw(Form("%s>>%s", expr, hname), cut, "goff");
    h->SetDirectory(nullptr);
    f->Close();
    return h;
}

// relleno de los datos con el ángulo recalculado:
//   1) offsets = medias de x0, x1, y0, y1 de los eventos que pasan 'cut'
//      (sin corte de ángulo), como el primer bucle de fillHistograms
//   2) eventos que pasan 'cut' y cos theta' recalculado > cosCut
static TH1D* fillManual(const char* file, const char* treeName, const char* dtExpr,
                        const char* cut, double cosCut, const char* hname, Offsets& off)
{
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) { printf("[ERROR] no se puede abrir %s\n", file); return nullptr; }
    TTree* t = (TTree*)f->Get(treeName);
    if (!t) { printf("[ERROR] no está el árbol %s en %s\n", treeName, file); f->Close(); return nullptr; }

    std::string fullCut = cut;
    if (t->GetBranch("RunNumber")) fullCut += std::string(" && ") + kDataRunCut;
    else printf("[AVISO] %s no tiene RunNumber: no se excluyen runs\n", file);

    if (TObject* o = gROOT->FindObject(hname)) delete o;
    TH1D* h = new TH1D(hname, "", kNb, kLo, kHi);
    h->SetDirectory(nullptr);
    h->Sumw2();

    bool ok = true;
    {   // las fórmulas deben destruirse antes de cerrar el fichero
        TTreeFormula fCut("fCut", fullCut.c_str(), t);
        TTreeFormula fX0 ("fX0", "x0", t), fX1("fX1", "x1", t);
        TTreeFormula fY0 ("fY0", "y0", t), fY1("fY1", "y1", t);
        TTreeFormula fDt ("fDt", dtExpr, t);
        for (TTreeFormula* fo : {&fCut, &fX0, &fX1, &fY0, &fY1, &fDt})
            if (fo->GetNdim() == 0) { printf("[ERROR] fórmula no válida en %s: %s\n", file, fo->GetTitle()); ok = false; }

        if (ok) {
            auto val = [](TTreeFormula& fo) { fo.GetNdata(); return fo.EvalInstance(0); };
            const Long64_t n = t->GetEntries();

            off = Offsets();
            for (Long64_t i = 0; i < n; ++i) {
                t->LoadTree(i);
                if (val(fCut) == 0) continue;
                off.x0 += val(fX0); off.x1 += val(fX1);
                off.y0 += val(fY0); off.y1 += val(fY1);
                ++off.n;
            }
            if (off.n > 0) { off.x0 /= off.n; off.x1 /= off.n; off.y0 /= off.n; off.y1 /= off.n; }

            for (Long64_t i = 0; i < n; ++i) {
                t->LoadTree(i);
                if (val(fCut) == 0) continue;
                const double c = cosThetaDet(val(fX0), val(fY0), val(fX1), val(fY1), off);
                if (c <= cosCut || c > 1.0) continue;
                h->Fill(val(fDt));
            }
        }
    }
    f->Close();
    if (!ok) { delete h; return nullptr; }
    return h;
}

// intervalo que contiene la fracción 'frac' central del histograma
static void centralInterval(TH1* h, double frac, int& b1, int& b2)
{
    const double tot = h->Integral();
    const double tail = 0.5 * (1.0 - frac) * tot;
    double s = 0; b1 = 1;
    for (int b = 1; b <= h->GetNbinsX(); ++b) { s += h->GetBinContent(b); if (s >= tail) { b1 = b; break; } }
    s = 0; b2 = h->GetNbinsX();
    for (int b = h->GetNbinsX(); b >= 1; --b) { s += h->GetBinContent(b); if (s >= tail) { b2 = b; break; } }
}

static void meanRms(TH1* h, int b1, int b2, double& m, double& r)
{
    double s = 0, sx = 0, sxx = 0;
    for (int b = b1; b <= b2; ++b) {
        const double w = h->GetBinContent(b), x = h->GetBinCenter(b);
        s += w; sx += w * x; sxx += w * x * x;
    }
    m = s > 0 ? sx / s : 0; r = s > 0 ? std::sqrt(std::max(0.0, sxx / s - m * m)) : 0;
}

// ------------------------------------------------------------------------------------
// un panel: relleno suave + contorno grueso + barras de error, todo en 'hex'
// ------------------------------------------------------------------------------------
static const int kW = 1700, kH = 680;          // lienzo 2 x 1
static const int kPadW = kW / 2;

static void drawPanel(TPad* p, TH1D* h, const char* hex, const char* tag,
                      const std::vector<std::string>& info)
{
    const int ink  = TColor::GetColor("#1B1B1B");
    const int grey = TColor::GetColor("#6E6E6E");
    const int col  = TColor::GetColor(hex);

    p->SetLeftMargin(130.0 / kPadW);
    p->SetRightMargin(30.0 / kPadW);
    p->SetTopMargin(40.0 / kH);
    p->SetBottomMargin(100.0 / kH);
    p->SetTickx(1);
    p->SetTicky(1);
    p->SetFillStyle(4000);
    p->SetBorderMode(0);
    p->Draw();
    p->cd();

    double ymax = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b)
        ymax = std::max(ymax, h->GetBinContent(b) + h->GetBinError(b));
    if (ymax <= 0) ymax = 1;
    const double yTop = 1.32 * ymax;                      // hueco arriba para el texto

    TH1F* fr = p->DrawFrame(kLo, 0.0, kHi, yTop);
    fr->GetXaxis()->SetTitle("tof_{1} #minus tof_{0}  (ns)");
    fr->GetYaxis()->SetTitle(Form("Events / %.2g ns", h->GetBinWidth(1)));
    TAxis* ax[2] = {fr->GetXaxis(), fr->GetYaxis()};
    for (TAxis* a : ax) {
        a->SetTitleFont(kFont); a->SetLabelFont(kFont);
        a->SetTitleSize(32);    a->SetLabelSize(28);
        a->SetAxisColor(ink);   a->SetLabelColor(ink); a->SetTitleColor(ink);
    }
    fr->GetXaxis()->SetNdivisions(510);
    fr->GetYaxis()->SetNdivisions(505);
    fr->GetXaxis()->SetTitleOffset(1.2);
    fr->GetYaxis()->SetTitleOffset(1.8);
    fr->GetXaxis()->SetTickLength(18.0 / kH);
    fr->GetYaxis()->SetTickLength(18.0 / kPadW);

    // referencia en 0
    TLine* l0 = new TLine(0.0, 0.0, 0.0, yTop);
    l0->SetLineColor(TColor::GetColor("#B8B8B8"));
    l0->SetLineStyle(3);
    l0->SetLineWidth(2);
    l0->Draw();

    TH1D* fill = (TH1D*)h->Clone(Form("%s_fill", h->GetName()));
    fill->SetDirectory(nullptr);
    fill->SetFillColorAlpha(col, 0.22);
    fill->SetFillStyle(1001);
    fill->SetLineWidth(0);
    fill->Draw("HIST SAME");

    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetLineWidth(3);
    h->SetMarkerSize(0);
    h->SetFillStyle(0);
    h->Draw("HIST ][ SAME");
    h->Draw("E1 SAME");

    // texto arriba a la izquierda: (a) + nombre en el color de la curva, resto en gris
    TLatex tx;
    tx.SetNDC();
    tx.SetTextAlign(13);
    const double x0 = 160.0 / kPadW;
    double y = 1.0 - 60.0 / kH;
    tx.SetTextFont(kFontBold); tx.SetTextSize(30); tx.SetTextColor(ink);
    tx.DrawLatex(x0, y, tag);
    tx.SetTextColor(col);
    tx.DrawLatex(x0 + 48.0 / kPadW, y, info[0].c_str());
    y -= 44.0 / kH;
    tx.SetTextFont(kFont); tx.SetTextSize(24); tx.SetTextColor(grey);
    for (size_t i = 1; i < info.size(); ++i) {
        tx.DrawLatex(x0, y, info[i].c_str());
        y -= 34.0 / kH;
    }

    p->RedrawAxis();
    p->Modified();
}

// ------------------------------------------------------------------------------------
void plot_ce_data_vs_mc()
{
    Offsets offD;
    TH1D* hD = fillManual(kDataFile, kDataTree, kDataDt, kDataCut, kDataCosCut, "h_ce_data", offD);
    TH1D* hM = fillFrom(kMcFile, kMcTree, kMcDt, kMcCut, "h_ce_mc");
    if (!hD || !hM) return;
    printf("Offsets datos (N = %lld): <x0> = %.3f  <x1> = %.3f  <y0> = %.3f  <y1> = %.3f\n",
           offD.n, offD.x0, offD.x1, offD.y0, offD.y1);
    const Long64_t nD = (Long64_t)hD->GetEntries(), nM = (Long64_t)hM->GetEntries();
    printf("Datos: %s && cos_theta_det(recalc) > %.2f\n  N = %lld\n", kDataCut, kDataCosCut, nD);
    printf("MC:    %s\n  N = %lld\n", kMcCut, nM);
    if (hD->Integral() <= 0 || hM->Integral() <= 0) { printf("[ERROR] histograma vacío\n"); return; }

    // media y rms de cada una en su 90 % central
    int d1, d2, m1, m2;
    centralInterval(hD, 0.90, d1, d2);
    centralInterval(hM, 0.90, m1, m2);
    double mD, rD, mM, rM;
    meanRms(hD, d1, d2, mD, rD);
    meanRms(hM, m1, m2, mM, rM);
    printf("  datos: media = %+.3f ns  rms = %.3f ns   (90 %% central)\n", mD, rD);
    printf("  MC:    media = %+.3f ns  rms = %.3f ns   (90 %% central)\n", mM, rM);
    printf("  desplazamiento datos - MC = %+.3f ns\n", mD - mM);

    // ---- lienzo 2 x 1 ------------------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetErrorX(0.);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetFrameLineColor(TColor::GetColor("#1B1B1B"));
    TGaxis::SetMaxDigits(3);

    if (TObject* o = gROOT->FindObject("c_ce_dvm")) delete o;
    TCanvas* c = new TCanvas("c_ce_dvm", "", kW, kH);
    c->SetFillColor(kWhite);
    c->SetBorderMode(0);

    c->cd();
    TPad* pL = new TPad("pL_ce_dvm", "", 0.0, 0.0, 0.5, 1.0);
    drawPanel(pL, hD, "#8C1C3E", "(a)",                      // granate
              {"Data",
               "^{nat}Ce,  cos#kern[0.15]{#theta'} > 0.8,  400 < #it{E}_{n} < 2000 MeV",
               Form("N = %lld,   #LTt#GT = %+.2f ns,   rms = %.2f ns", nD, mD, rD)});

    c->cd();
    TPad* pR = new TPad("pR_ce_dvm", "", 0.5, 0.0, 1.0, 1.0);
    drawPanel(pR, hM, "#1F5FA8", "(b)",                      // azul
              {"Simulation",
               "^{nat}Ce,  cos#kern[0.15]{#theta'} > 0.8,  #it{E}_{n} > 400 MeV",
               Form("N = %lld,   #LTt#GT = %+.2f ns,   rms = %.2f ns", nM, mM, rM)});

    c->cd();
    c->Modified();
    c->Update();
    c->SaveAs("Ce_tof_data_vs_mc.pdf");
    c->SaveAs("Ce_tof_data_vs_mc.png");
}