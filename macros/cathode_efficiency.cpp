// =============================================================================
//  cathode_efficiency_pretty.C
//
//  Análisis de eficiencia de cátodos + estudio del efecto del backing de Al
//  sobre la detección de fragmentos de fisión.
//
//  Uso:   root -l 'cathode_efficiency_pretty.C'
//         (o  .x cathode_efficiency_pretty.C  desde el prompt de ROOT)
//
//  Idea física central
//  -------------------
//  dt = tof1 - tof0.  dt > 0  ->  el fragmento que va al detector 1 es el LENTO
//  (el pesado).  Si el backing de Al está del lado del detector 1, el fragmento
//  pesado que lo atraviesa pierde energía y sale con menos ionización / más
//  straggling  ->  los cátodos X1,Y1 fallan preferentemente a dt > 0.
//  Por simetría del montaje, SIN backing se debería cumplir
//
//        eps_1(+dt)  ==  eps_0(-dt)
//
//  y toda la desviación de esa igualdad es el efecto del backing. Ese es el
//  "mirror test" (canvas 03) y su doble cociente R(dt) (panel inferior).
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

// =============================================================================
//  CONFIGURACIÓN  (lo único que deberías tocar)
// =============================================================================
namespace CFG {
    const char* kFile = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/"
                        "coincidences_final.root";
    const char* kTree = "events_cerium";          // <-- ojo: en tu macro la variable
                                                //     se llamaba events_uranium
    const char* kBase = "neutron_energy>400 && neutron_energy<2000 && amp0>8000 && amp1>8000 && abs(amp1-amp0)/(amp1+amp0)<0.5 && det0<6";

    const std::vector<int> kExcluded = {118771, 118789};

    const char* kDT   = "tof1-tof0";            // variable de asimetría temporal
    const double kDTmin = -10., kDTmax = 10.;
    const int    kDTbins = 80;

    const double kEnMin = 400., kEnMax = 2000.;    // MeV (o lo que uses)
    const int    kEnBins = 5;
    const bool   kEnLog  = true;

    const double kAmpMin = 16000., kAmpMax = 40000.;
    const int    kAmpBins = 40;

    const char* kOut = ".";                     // directorio de salida
    const bool  kSavePDF = true;
    const bool  kSavePNG = true;

    const int   kMinEntries2D = 10;             // stats mínimas por celda en 2D
}

// =============================================================================
//  ESTILO
// =============================================================================

// Paleta Okabe-Ito: bonita, imprimible en B/N y segura para daltonismo
Int_t NiceColor(int i) {
    static const char* hex[8] = {"#0072B2",  // azul
                                 "#D55E00",  // bermellón
                                 "#009E73",  // verde
                                 "#CC79A7",  // rosa
                                 "#E69F00",  // naranja
                                 "#56B4E9",  // celeste
                                 "#8C564B",  // marrón
                                 "#4C4C4C"}; // gris oscuro
    return TColor::GetColor(hex[i % 8]);
}

void SetPrettyStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.09);

    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetTitleFont(42, "");
    gStyle->SetLabelSize(0.045, "xyz");
    gStyle->SetTitleSize(0.050, "xyz");
    gStyle->SetTitleSize(0.055, "");
    gStyle->SetTitleOffset(1.15, "y");
    gStyle->SetTitleOffset(1.05, "x");
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(0.5);

    gStyle->SetHistLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(1.0);
    gStyle->SetEndErrorSize(0);

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendTextSize(0.038);

    gStyle->SetGridColor(kGray + 1);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);

    gStyle->SetPalette(kViridis);
    gStyle->SetNumberContours(255);
    TColor::InvertPalette();      // viridis con el claro abajo -> más legible en COLZ
    TColor::InvertPalette();      // (doble inversión = original; deja la línea para
                                  //  probar la orientación que más te guste)
}

// =============================================================================
//  UTILIDADES
// =============================================================================

bool HasVar(TTree* t, const char* v) {
    return (t->GetBranch(v) != nullptr) || (t->GetLeaf(v) != nullptr);
}

void Save(TCanvas* c, const char* base) {
    if (CFG::kSavePDF) c->SaveAs(Form("%s/%s.pdf", CFG::kOut, base));
    if (CFG::kSavePNG) c->SaveAs(Form("%s/%s.png", CFG::kOut, base));
}

TH1D* Proj1D(TTree* t, const char* name, const char* var, const char* cut,
             int nb, double lo, double hi) {
    if (auto* old = (TH1*)gDirectory->Get(name)) delete old;
    TH1D* h = new TH1D(name, "", nb, lo, hi);
    t->Project(name, var, cut);
    return h;
}

TH1D* Proj1DVar(TTree* t, const char* name, const char* var, const char* cut,
                const std::vector<double>& edges) {
    if (auto* old = (TH1*)gDirectory->Get(name)) delete old;
    TH1D* h = new TH1D(name, "", (int)edges.size() - 1, edges.data());
    t->Project(name, var, cut);
    return h;
}

std::vector<double> LogBins(double lo, double hi, int n) {
    std::vector<double> e(n + 1);
    const double s = std::log(hi / lo) / n;
    for (int i = 0; i <= n; ++i) e[i] = lo * std::exp(s * i);
    return e;
}

// Construye una TEfficiency con intervalos de Clopper-Pearson
TEfficiency* MakeEff(TTree* t, const char* tag, const char* var,
                     const char* cutTot, const char* cutPass,
                     int nb, double lo, double hi,
                     const std::vector<double>* edges = nullptr) {
    TString nP = TString::Format("hP_%s", tag);
    TString nT = TString::Format("hT_%s", tag);
    TString cP = TString::Format("(%s) && (%s)", cutTot, cutPass);

    TH1D *hp, *ht;
    if (edges) {
        hp = Proj1DVar(t, nP, var, cP, *edges);
        ht = Proj1DVar(t, nT, var, cutTot, *edges);
    } else {
        hp = Proj1D(t, nP, var, cP, nb, lo, hi);
        ht = Proj1D(t, nT, var, cutTot, nb, lo, hi);
    }

    if (!TEfficiency::CheckConsistency(*hp, *ht)) {
        std::cout << "[MakeEff] inconsistencia en " << tag << std::endl;
        return nullptr;
    }
    TEfficiency* e = new TEfficiency(*hp, *ht);
    e->SetName(Form("eff_%s", tag));
    e->SetStatisticOption(TEfficiency::kFCP);   // Clopper-Pearson
    e->SetConfidenceLevel(0.683);
    delete hp;
    delete ht;
    return e;
}

void StyleEff(TEfficiency* e, int ci, int mstyle = 20) {
    if (!e) return;
    e->SetLineColor(NiceColor(ci));
    e->SetMarkerColor(NiceColor(ci));
    e->SetMarkerStyle(mstyle);
    e->SetLineWidth(2);
    e->SetMarkerSize(1.1);
}

void FixEffAxis(TEfficiency* e, double ymin, double ymax) {
    gPad->Update();
    if (e && e->GetPaintedGraph())
        e->GetPaintedGraph()->GetYaxis()->SetRangeUser(ymin, ymax);
    gPad->Update();
}

// Cociente punto a punto de dos TEfficiency con el mismo binning.
// mirror = true  -> compara el bin i de 'a' con el bin (N+1-i) de 'b'
//                   (es decir eps_a(+dt) / eps_b(-dt))
TGraphAsymmErrors* EffRatio(TEfficiency* a, TEfficiency* b, bool mirror = false) {
    if (!a || !b) return nullptr;
    const TH1* h = a->GetTotalHistogram();
    const int nb = h->GetNbinsX();
    auto* g = new TGraphAsymmErrors();
    int k = 0;
    for (int i = 1; i <= nb; ++i) {
        const int j = mirror ? (nb + 1 - i) : i;
        const double ea = a->GetEfficiency(i);
        const double eb = b->GetEfficiency(j);
        if (ea <= 0 || eb <= 0) continue;
        const double sa = 0.5 * (a->GetEfficiencyErrorUp(i) + a->GetEfficiencyErrorLow(i));
        const double sb = 0.5 * (b->GetEfficiencyErrorUp(j) + b->GetEfficiencyErrorLow(j));
        const double r  = ea / eb;
        const double sr = r * std::sqrt((sa / ea) * (sa / ea) + (sb / eb) * (sb / eb));
        g->SetPoint(k, h->GetBinCenter(i), r);
        g->SetPointError(k, 0, 0, sr, sr);
        ++k;
    }
    return g;
}

void DrawHeader(const char* txt) {
    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.035);
    l.SetTextColor(kGray + 2);
    l.DrawLatex(0.14, 0.945, txt);
}

// =============================================================================
//  MACRO PRINCIPAL
// =============================================================================
void cathode_efficiency() {

    SetPrettyStyle();

    // ---------------------------------------------------------------- entrada
    TFile* fin = TFile::Open(CFG::kFile);
    if (!fin || fin->IsZombie()) { std::cout << "ERROR: no se abre el fichero\n"; return; }
    TTree* T = (TTree*)fin->Get(CFG::kTree);
    if (!T) { std::cout << "ERROR: no se encontró el TTree " << CFG::kTree << "\n"; return; }

    // -------------------------------------------------------------- selección
    TString runCut;
    for (size_t i = 0; i < CFG::kExcluded.size(); ++i) {
        if (i) runCut += " && ";
        runCut += TString::Format("RunNumber != %d", CFG::kExcluded[i]);
    }
    const TString BASE = TString::Format("(%s) && (%s)", CFG::kBase, runCut.Data());

    std::cout << "\nSelección base: " << BASE << "\n";
    const Long64_t Ntot = T->GetEntries(BASE);
    std::cout << "Eventos totales: " << Ntot << "\n";
    if (Ntot == 0) { std::cout << "Sin eventos, abortando.\n"; return; }

    // Cortes de "pasa"
    const char* pX0 = "hasX0";
    const char* pY0 = "hasY0";
    const char* pX1 = "hasX1";
    const char* pY1 = "hasY1";
    const char* pD0 = "hasX0 && hasY0";          // posición reconstruible en det 0
    const char* pD1 = "hasX1 && hasY1";          // posición reconstruible en det 1
    const char* pALL = "hasX0 && hasY0 && hasX1 && hasY1";

    // =========================================================================
    //  1) TABLA DE LAS 16 COMBINACIONES  (una sola pasada por el árbol)
    // =========================================================================
    // El patrón se codifica como  X0 + 2*Y0 + 4*X1 + 8*Y1  -> un histograma de
    // 16 bins en vez de 16 llamadas a GetEntries(). Mucho más rápido.
    TH1D* hPat = Proj1D(T, "hPat", "hasX0 + 2*hasY0 + 4*hasX1 + 8*hasY1",
                        BASE, 16, -0.5, 15.5);

    auto PatLabel = [](int p) {
        return TString::Format("%d%d%d%d", (p >> 0) & 1, (p >> 1) & 1,
                                            (p >> 2) & 1, (p >> 3) & 1);
    };

    std::vector<Long64_t> n(16, 0);
    for (int p = 0; p < 16; ++p) n[p] = (Long64_t)hPat->GetBinContent(p + 1);

    printf("\n=================== 16 COMBINACIONES  (X0 Y0 X1 Y1) ===================\n");
    printf("%-8s %12s %10s %12s\n", "Patrón", "Eventos", "Fracción", "Acumulada");
    printf("-----------------------------------------------------------------------\n");
    Long64_t acc = 0;
    // Orden descendente en estadística, que es como de verdad se lee la tabla
    std::vector<int> order(16);
    for (int i = 0; i < 16; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int a, int b) { return n[a] > n[b]; });
    for (int k = 0; k < 16; ++k) {
        const int p = order[k];
        acc += n[p];
        printf("%-8s %12lld %9.3f%% %11.3f%%\n", PatLabel(p).Data(), n[p],
               100. * n[p] / Ntot, 100. * acc / Ntot);
    }
    printf("-----------------------------------------------------------------------\n");
    printf("Suma: %lld / %lld  (debe ser 100%%: los 16 patrones son exhaustivos)\n",
           acc, Ntot);

    // -------------------------------------------------- eficiencias marginales
    Long64_t nX0 = 0, nY0 = 0, nX1 = 0, nY1 = 0, nD0 = 0, nD1 = 0, nALL = 0;
    for (int p = 0; p < 16; ++p) {
        const bool x0 = (p >> 0) & 1, y0 = (p >> 1) & 1;
        const bool x1 = (p >> 2) & 1, y1 = (p >> 3) & 1;
        if (x0) nX0 += n[p];
        if (y0) nY0 += n[p];
        if (x1) nX1 += n[p];
        if (y1) nY1 += n[p];
        if (x0 && y0) nD0 += n[p];
        if (x1 && y1) nD1 += n[p];
        if (x0 && y0 && x1 && y1) nALL += n[p];
    }
    auto binErr = [&](Long64_t k) { // error binomial simple, solo para la tabla
        const double e = (double)k / Ntot;
        return 100. * std::sqrt(e * (1 - e) / Ntot);
    };
    printf("\n======================= EFICIENCIAS MARGINALES ========================\n");
    printf("X0        %10lld   %6.2f %% +- %.2f\n", nX0, 100.*nX0/Ntot, binErr(nX0));
    printf("Y0        %10lld   %6.2f %% +- %.2f\n", nY0, 100.*nY0/Ntot, binErr(nY0));
    printf("X1        %10lld   %6.2f %% +- %.2f\n", nX1, 100.*nX1/Ntot, binErr(nX1));
    printf("Y1        %10lld   %6.2f %% +- %.2f\n", nY1, 100.*nY1/Ntot, binErr(nY1));
    printf("-----------------------------------------------------------------------\n");
    printf("Det0 (X0&Y0) %7lld   %6.2f %%\n", nD0, 100.*nD0/Ntot);
    printf("Det1 (X1&Y1) %7lld   %6.2f %%\n", nD1, 100.*nD1/Ntot);
    printf("Los 4        %7lld   %6.2f %%\n", nALL, 100.*nALL/Ntot);
    printf("\n>>> Asimetría det1/det0 = %.4f   (== 1 si el backing no afectara)\n",
           (double)nD1 / nD0);

    // =========================================================================
    //  CANVAS 01 - Poblaciones de los 16 patrones
    // =========================================================================
    TCanvas* c01 = new TCanvas("c01", "patrones", 1100, 700);
    c01->SetGridy();
    TH1D* hBar = new TH1D("hBar", ";X0 Y0 X1 Y1;fraction of events [%]",
                          16, 0, 16);
    for (int k = 0; k < 16; ++k) {
        const int p = order[k];
        hBar->SetBinContent(k + 1, 100. * n[p] / Ntot);
        hBar->GetXaxis()->SetBinLabel(k + 1, PatLabel(p));
    }
    hBar->SetFillColorAlpha(NiceColor(0), 0.75);
    hBar->SetLineColor(NiceColor(0));
    hBar->SetBarWidth(0.8);
    hBar->SetBarOffset(0.1);
    hBar->GetXaxis()->LabelsOption("v");
    hBar->GetXaxis()->SetLabelSize(0.040);
    hBar->SetMinimum(0);
    hBar->SetMaximum(1.25 * hBar->GetMaximum());
    hBar->Draw("BAR");
    {
        TLatex l; l.SetTextFont(42); l.SetTextSize(0.028); l.SetTextAlign(21);
        l.SetTextColor(kGray + 3);
        for (int k = 0; k < 16; ++k)
            if (hBar->GetBinContent(k + 1) > 0.05)
                l.DrawLatex(k + 0.5, hBar->GetBinContent(k + 1) * 1.03,
                            Form("%.1f", hBar->GetBinContent(k + 1)));
    }
    DrawHeader(Form("%s  |  N = %lld  |  runs excluidos: 118771, 118789",
                    CFG::kTree, Ntot));
    Save(c01, "01_patrones");

    // =========================================================================
    //  CANVAS 02 - eps(dt) para cada cátodo por separado
    //  Test directo de tu hipótesis: X1,Y1 deben caer a dt > 0
    //  (fragmento pesado hacia el detector 1, el que atraviesa el backing)
    // =========================================================================
    TEfficiency* eX0 = MakeEff(T, "X0", CFG::kDT, BASE, pX0, CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);
    TEfficiency* eY0 = MakeEff(T, "Y0", CFG::kDT, BASE, pY0, CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);
    TEfficiency* eX1 = MakeEff(T, "X1", CFG::kDT, BASE, pX1, CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);
    TEfficiency* eY1 = MakeEff(T, "Y1", CFG::kDT, BASE, pY1, CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);

    StyleEff(eX0, 0, 20); StyleEff(eY0, 5, 24);
    StyleEff(eX1, 1, 21); StyleEff(eY1, 4, 25);

    TCanvas* c02 = new TCanvas("c02", "eff vs dt", 1100, 750);
    c02->SetGridx(); c02->SetGridy();
    if (eX0) {
        eX0->SetTitle(Form(""
                           "#Deltat = tof_{1} - tof_{0}  [ns];#varepsilon"));
        eX0->Draw("AP");
        if (eY0) eY0->Draw("SAME P");
        if (eX1) eX1->Draw("SAME P");
        if (eY1) eY1->Draw("SAME P");
        FixEffAxis(eX0, 0., 1.05);
    }
    TLine* l0 = new TLine(0, 0, 0, 1.05);
    l0->SetLineStyle(2); l0->SetLineColor(kGray + 2); l0->Draw();
    {
        TLegend* lg = new TLegend(0.60, 0.18, 0.92, 0.42);
        lg->SetNColumns(2);
        if (eX0) lg->AddEntry(eX0, "X0", "lp");
        if (eY0) lg->AddEntry(eY0, "Y0", "lp");
        if (eX1) lg->AddEntry(eX1, "X1", "lp");
        if (eY1) lg->AddEntry(eY1, "Y1", "lp");
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.032);
        t.SetTextColor(kGray + 3);
        t.DrawLatex(0.17, 0.24, "#Deltat < 0: heavy #rightarrow det0");
        t.DrawLatex(0.17, 0.19, "#Deltat > 0: heavy #rightarrow det1 (backing)");
    }
    Save(c02, "02_eff_vs_dt");

    // =========================================================================
    //  CANVAS 03 - MIRROR TEST  (el plot clave)
    //  eps_det1(+dt)  vs  eps_det0(-dt).  Iguales <=> no hay efecto de backing.
    // =========================================================================
    TEfficiency* eD1 = MakeEff(T, "D1", CFG::kDT, BASE, pD1,
                               CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);
    TEfficiency* eD0m = MakeEff(T, "D0m", "tof0-tof1", BASE, pD0,   // eje espejado
                                CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);
    StyleEff(eD1, 1, 21);
    StyleEff(eD0m, 0, 20);

    TCanvas* c03 = new TCanvas("c03", "", 1100, 900);
    TPad* pTop = new TPad("pTop", "", 0, 0.32, 1, 1);
    TPad* pBot = new TPad("pBot", "", 0, 0.00, 1, 0.33);
    pTop->SetBottomMargin(0.02); pTop->SetGridx(); pTop->SetGridy();
    pBot->SetTopMargin(0.03); pBot->SetBottomMargin(0.32);
    pBot->SetGridx(); pBot->SetGridy();
    pTop->Draw(); pBot->Draw();

    pTop->cd();
    if (eD1) {
        eD1->Draw("AP");
        if (eD0m) eD0m->Draw("SAME P");
        FixEffAxis(eD1, 0., 1.05);
        if (eD1->GetPaintedGraph())
            eD1->GetPaintedGraph()->GetXaxis()->SetLabelSize(0);
    }
    {
        TLegend* lg = new TLegend(0.17, 0.16, 0.55, 0.34);
        if (eD1)  lg->AddEntry(eD1,  "#varepsilon_{det1}(+#Deltat)  (con backing)", "lp");
        if (eD0m) lg->AddEntry(eD0m, "#varepsilon_{det0}(-#Deltat)  (sin backing)", "lp");
        lg->Draw();
    }

    pBot->cd();
    TGraphAsymmErrors* gR = EffRatio(eD1, eD0m, /*mirror=*/false);
    if (gR) {
        gR->SetTitle(";#Deltat = tof_{1} - tof_{0}  [ns];"
                     "#varepsilon_{1}/#varepsilon_{0}");
        gR->SetMarkerStyle(20);
        gR->SetMarkerColor(NiceColor(2));
        gR->SetLineColor(NiceColor(2));
        gR->SetLineWidth(2);
        gR->GetYaxis()->SetNdivisions(505);
        gR->GetYaxis()->SetTitleSize(0.10);
        gR->GetYaxis()->SetLabelSize(0.085);
        gR->GetYaxis()->SetTitleOffset(0.55);
        gR->GetXaxis()->SetTitleSize(0.11);
        gR->GetXaxis()->SetLabelSize(0.095);
        gR->GetXaxis()->SetTitleOffset(1.20);
        gR->GetYaxis()->SetRangeUser(0.70, 1.30);
        gR->GetXaxis()->SetLimits(CFG::kDTmin, CFG::kDTmax);
        gR->Draw("AP");
        TLine* one = new TLine(CFG::kDTmin, 1., CFG::kDTmax, 1.);
        one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw();
    }
    Save(c03, "03_mirror_test_backing");

    // =========================================================================
    //  CANVAS 04 - Formas del espectro dt: total / seleccionado / perdido
    //  Cuantifica el SESGO que introduce exigir cátodos.
    // =========================================================================
    TH1D* hAll  = Proj1D(T, "hAll",  CFG::kDT, BASE, CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);
    TH1D* hSel  = Proj1D(T, "hSel",  CFG::kDT,
                         TString::Format("(%s)&&(%s)", BASE.Data(), pALL),
                         CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);
    TH1D* hLost = Proj1D(T, "hLost", CFG::kDT,
                         TString::Format("(%s)&&!(%s)", BASE.Data(), pALL),
                         CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);
    TH1D* hL1   = Proj1D(T, "hL1", CFG::kDT,
                         TString::Format("(%s)&&!(%s)&&(%s)", BASE.Data(), pD1, pD0),
                         CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);
    TH1D* hL0   = Proj1D(T, "hL0", CFG::kDT,
                         TString::Format("(%s)&&!(%s)&&(%s)", BASE.Data(), pD0, pD1),
                         CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);

    TCanvas* c04 = new TCanvas("c04", "formas dt", 1150, 750);
    c04->SetGridx(); c04->SetGridy();
    auto norm = [](TH1D* h) { if (h && h->Integral() > 0) h->Scale(1. / h->Integral()); };
    TH1D* hAllN  = (TH1D*)hAll->Clone("hAllN");   norm(hAllN);
    TH1D* hSelN  = (TH1D*)hSel->Clone("hSelN");   norm(hSelN);
    TH1D* hL1N   = (TH1D*)hL1->Clone("hL1N");     norm(hL1N);
    TH1D* hL0N   = (TH1D*)hL0->Clone("hL0N");     norm(hL0N);

    hAllN->SetTitle(""
                    "#Deltat = tof_{1} - tof_{0}  [ns];Density");
    hAllN->SetLineColor(kGray + 3);
    hAllN->SetFillColorAlpha(kGray + 1, 0.35);
    hAllN->SetLineWidth(2);
    hAllN->SetMaximum(1.45 * hAllN->GetMaximum());
    hAllN->Draw("HIST");

    hSelN->SetLineColor(NiceColor(2)); hSelN->SetLineWidth(3); hSelN->Draw("HIST SAME");
    hL1N->SetLineColor(NiceColor(1));  hL1N->SetLineWidth(3);
    hL1N->SetLineStyle(2);             hL1N->Draw("HIST SAME");
    hL0N->SetLineColor(NiceColor(0));  hL0N->SetLineWidth(3);
    hL0N->SetLineStyle(7);             hL0N->Draw("HIST SAME");
    {
        TLegend* lg = new TLegend(0.60, 0.68, 0.94, 0.90);
        lg->AddEntry(hAllN, "all events", "f");
        lg->AddEntry(hSelN, "with 4 validated cathodes", "l");
        lg->AddEntry(hL1N,  "det1 missing", "l");
        lg->AddEntry(hL0N,  "det0 missing", "l");
        lg->Draw();
    }
    Save(c04, "04_formas_dt");

    printf("\n=================== SESGO INTRODUCIDO POR LOS CÁTODOS =================\n");
    printf("<dt> todos            = %8.3f ns  (RMS %.3f)\n",
           hAll->GetMean(), hAll->GetRMS());
    printf("<dt> con 4 cátodos    = %8.3f ns  (RMS %.3f)\n",
           hSel->GetMean(), hSel->GetRMS());
    printf("<dt> perdidos         = %8.3f ns  (RMS %.3f)\n",
           hLost->GetMean(), hLost->GetRMS());
    printf("Desplazamiento del centroide: %.3f ns\n",
           hSel->GetMean() - hAll->GetMean());
    printf("Fracción perdida global: %.2f %%\n", 100. * hLost->Integral() / hAll->Integral());

    // =========================================================================
    //  CANVAS 05 - eps vs amplitud del ánodo (proxy de energía del fragmento)
    //  El "turn-on" del cátodo se ve directamente. Ojo a la correlación:
    //  amp bajo = fragmento pesado (o degradado por el backing).
    // =========================================================================
    TCanvas* c05 = new TCanvas("c05", "eff vs amp", 1300, 620);
    c05->Divide(2, 1);

    struct AmpCase { const char* var; const char* pass; const char* tag; const char* ttl; };
    AmpCase acases[2] = {
        {"amp0", pD0, "A0", "Detector 0 (no backing)"},
        {"amp1", pD1, "A1", "Detector 1 (backing)"}
    };
    for (int i = 0; i < 2; ++i) {
        c05->cd(i + 1); gPad->SetGridx(); gPad->SetGridy();
        TEfficiency* e = MakeEff(T, acases[i].tag, acases[i].var, BASE, acases[i].pass,
                                 CFG::kAmpBins, CFG::kAmpMin, CFG::kAmpMax);
        StyleEff(e, i, 20 + i);
        if (e) {
            e->SetTitle(Form("%s;%s  [canales];#varepsilon(X #upoint Y)",
                             acases[i].ttl, acases[i].var));
            e->Draw("AP");
            FixEffAxis(e, 0., 1.05);
        }
    }
    Save(c05, "05_eff_vs_amplitud");

    // =========================================================================
    //  CANVAS 06 - eps vs asimetría de amplitudes  (eje de asimetría de masa)
    // =========================================================================
    const char* asym = "(amp0-amp1)/(amp0+amp1)";
    TEfficiency* eAsD0 = MakeEff(T, "AsD0", asym, BASE, pD0, 30, -0.6, 0.6);
    TEfficiency* eAsD1 = MakeEff(T, "AsD1", asym, BASE, pD1, 30, -0.6, 0.6);
    StyleEff(eAsD0, 0, 20); StyleEff(eAsD1, 1, 21);

    TCanvas* c06 = new TCanvas("c06", "eff vs asimetria", 1100, 720);
    c06->SetGridx(); c06->SetGridy();
    if (eAsD0) {
        eAsD0->SetTitle(""
                        "(amp_{0}-amp_{1})/(amp_{0}+amp_{1});#varepsilon(X #upoint Y)");
        eAsD0->Draw("AP");
        if (eAsD1) eAsD1->Draw("SAME P");
        FixEffAxis(eAsD0, 0., 1.05);
        TLegend* lg = new TLegend(0.62, 0.18, 0.93, 0.34);
        lg->AddEntry(eAsD0, "det.0 (X0 #upoint Y0)", "lp");
        if (eAsD1) lg->AddEntry(eAsD1, "det.1 (X1 #upoint Y1)", "lp");
        lg->Draw();
    }
    Save(c06, "06_eff_vs_asimetria");

    // =========================================================================
    //  CANVAS 07 - Mapa 2D de eficiencia en (amp0, amp1)
    //  Dónde exactamente se pierden los eventos en el plano de energías.
    // =========================================================================
    auto Map2D = [&](const char* name, const char* pass, const char* ttl) -> TH2D* {
        TString nT = TString::Format("hT2_%s", name), nP = TString::Format("hP2_%s", name);
        if (auto* o = (TH2*)gDirectory->Get(nT)) delete o;
        if (auto* o = (TH2*)gDirectory->Get(nP)) delete o;
        TH2D* ht = new TH2D(nT, "", 40, CFG::kAmpMin, CFG::kAmpMax,
                                     40, CFG::kAmpMin, CFG::kAmpMax);
        TH2D* hp = new TH2D(nP, "", 40, CFG::kAmpMin, CFG::kAmpMax,
                                     40, CFG::kAmpMin, CFG::kAmpMax);
        T->Project(nT, "amp1:amp0", BASE);
        T->Project(nP, "amp1:amp0", TString::Format("(%s)&&(%s)", BASE.Data(), pass));
        TH2D* hr = (TH2D*)hp->Clone(Form("hR2_%s", name));
        hr->Divide(hp, ht, 1, 1, "B");
        for (int ix = 1; ix <= ht->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= ht->GetNbinsY(); ++iy)
                if (ht->GetBinContent(ix, iy) < CFG::kMinEntries2D)
                    hr->SetBinContent(ix, iy, 0);
        hr->SetTitle(Form("%s;amp_{0}  [arb. units];amp_{1}  [arb.units];#varepsilon", ttl));
        hr->SetMinimum(1e-6);
        hr->SetMaximum(1.0);
        hr->SetContour(255);
        return hr;
    };

    TCanvas* c07 = new TCanvas("c07", "2D map", 1400, 620);
    c07->Divide(2, 1);
    c07->cd(1); gPad->SetRightMargin(0.16);
    TH2D* m0 = Map2D("D0", pD0, "#varepsilon detector 0");
    m0->Draw("COLZ");
    c07->cd(2); gPad->SetRightMargin(0.16);
    TH2D* m1 = Map2D("D1", pD1, "#varepsilon detector 1 (backing)");
    m1->Draw("COLZ");
    Save(c07, "07_mapa_eff_amp0_amp1");

    // Diferencia de mapas: eps0 - eps1, el efecto neto del backing en el plano
    TCanvas* c07b = new TCanvas("c07b", "", 800, 650);
    gPad->SetRightMargin(0.16);
    TH2D* mdiff = (TH2D*)m0->Clone("mdiff");
    mdiff->Add(m1, -1);
    mdiff->SetTitle("#varepsilon_{det0} - #varepsilon_{det1};"
                    "amp_{0}  [arb. units];amp_{1}  [arb. units];#Delta#varepsilon");
    mdiff->SetMinimum(-0.5); mdiff->SetMaximum(0.5);
    gStyle->SetPalette(87);   // kLightTemperature / divergente
    mdiff->Draw("COLZ");
    Save(c07b, "07b_diferencia_mapas");
    gStyle->SetPalette(kViridis);

    // =========================================================================
    //  CANVAS 08 - eps vs energía de neutrón
    //  Si la anisotropía cambia con En, el camino efectivo en el Al cambia y la
    //  corrección deja de ser una constante -> impacto directo en sigma(En).
    // =========================================================================
    std::vector<double> enB = CFG::kEnLog
        ? LogBins(CFG::kEnMin, CFG::kEnMax, CFG::kEnBins)
        : [&] { std::vector<double> v(CFG::kEnBins + 1);
                for (int i = 0; i <= CFG::kEnBins; ++i)
                    v[i] = CFG::kEnMin + i * (CFG::kEnMax - CFG::kEnMin) / CFG::kEnBins;
                return v; }();

    TEfficiency* eEn0 = MakeEff(T, "En0", "neutron_energy", BASE, pD0, 0, 0, 0, &enB);
    TEfficiency* eEn1 = MakeEff(T, "En1", "neutron_energy", BASE, pD1, 0, 0, 0, &enB);
    TEfficiency* eEnA = MakeEff(T, "EnA", "neutron_energy", BASE, pALL, 0, 0, 0, &enB);
    StyleEff(eEn0, 0, 20); StyleEff(eEn1, 1, 21); StyleEff(eEnA, 2, 22);

    TCanvas* c08 = new TCanvas("c08", "eff vs En", 1100, 720);
    if (CFG::kEnLog) c08->SetLogx();
    c08->SetGridx(); c08->SetGridy();
    if (eEn0) {
        eEn0->SetTitle("Eficiencia vs energ#acute{i}a de neutr#acute{o}n;"
                       "E_{n}  [MeV];#varepsilon");
        eEn0->Draw("AP");
        if (eEn1) eEn1->Draw("SAME P");
        if (eEnA) eEnA->Draw("SAME P");
        FixEffAxis(eEn0, 0., 1.05);
        TLegend* lg = new TLegend(0.62, 0.16, 0.93, 0.36);
        lg->AddEntry(eEn0, "det.0 (X0 #upoint Y0)", "lp");
        if (eEn1) lg->AddEntry(eEn1, "det.1 (X1 #upoint Y1)", "lp");
        if (eEnA) lg->AddEntry(eEnA, "los 4 c#acute{a}todos", "lp");
        lg->Draw();
    }
    Save(c08, "08_eff_vs_energia_neutron");

    // =========================================================================
    //  CANVAS 09 - eps vs RunNumber (estabilidad; justifica las exclusiones)
    //  Aquí NO se excluye ningún run, precisamente para verlos.
    // =========================================================================
    const TString BASE_ALLRUNS = TString::Format("(%s)", CFG::kBase);
    TH1D* hR = Proj1D(T, "hRrange", "RunNumber", BASE_ALLRUNS, 1000, 0, 200000);
    int firstBin = hR->FindFirstBinAbove(0), lastBin = hR->FindLastBinAbove(0);
    double rlo = hR->GetBinLowEdge(firstBin) - 1;
    double rhi = hR->GetBinLowEdge(lastBin + 1) + 1;
    int nRuns = TMath::Max(10, (int)(rhi - rlo));
    if (nRuns > 400) nRuns = 400;

    TEfficiency* eR0 = MakeEff(T, "R0", "RunNumber", BASE_ALLRUNS, pD0, nRuns, rlo, rhi);
    TEfficiency* eR1 = MakeEff(T, "R1", "RunNumber", BASE_ALLRUNS, pD1, nRuns, rlo, rhi);
    StyleEff(eR0, 0, 20); StyleEff(eR1, 1, 21);

    TCanvas* c09 = new TCanvas("c09", "eff vs run", 1300, 650);
    c09->SetGridx(); c09->SetGridy();
    if (eR0) {
        eR0->SetTitle("Estabilidad por run (sin excluir nada);RunNumber;#varepsilon");
        eR0->Draw("AP");
        if (eR1) eR1->Draw("SAME P");
        FixEffAxis(eR0, 0., 1.05);
        TLatex t; t.SetTextFont(42); t.SetTextSize(0.030);
        t.SetTextColor(NiceColor(1)); t.SetTextAngle(90);
        for (int r : CFG::kExcluded) {
            TLine* lv = new TLine(r, 0, r, 1.05);
            lv->SetLineColor(NiceColor(1)); lv->SetLineStyle(2); lv->Draw();
            t.DrawLatex(r + 0.5, 0.10, Form("excluido %d", r));
        }
        TLegend* lg = new TLegend(0.62, 0.16, 0.93, 0.32);
        lg->AddEntry(eR0, "det.0", "lp");
        if (eR1) lg->AddEntry(eR1, "det.1", "lp");
        lg->Draw();
    }
    Save(c09, "09_eff_vs_run");

    // =========================================================================
    //  CANVAS 10 - dt de cada patrón (la versión bonita de tu 4x4)
    // =========================================================================
    TCanvas* c10 = new TCanvas("c10", "dt per combination of cathodes", 1500, 1150);
    c10->Divide(4, 4, 0.001, 0.001);
    for (int k = 0; k < 16; ++k) {
        const int p = order[k];
        c10->cd(k + 1);
        gPad->SetGridx(); gPad->SetGridy();
        gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.16);
        TString cut = TString::Format("(%s) && (hasX0 + 2*hasY0 + 4*hasX1 + 8*hasY1)==%d",
                                      BASE.Data(), p);
        TH1D* h = Proj1D(T, Form("hdt_%d", p), CFG::kDT, cut,
                         CFG::kDTbins, CFG::kDTmin, CFG::kDTmax);
        h->SetTitle(Form("%s   (N = %lld);#Deltat [ns];cuentas", PatLabel(p).Data(), n[p]));
        h->SetTitleSize(0.075);
        h->GetXaxis()->SetTitleSize(0.065); h->GetXaxis()->SetLabelSize(0.055);
        h->GetYaxis()->SetTitleSize(0.065); h->GetYaxis()->SetLabelSize(0.055);
        h->GetYaxis()->SetTitleOffset(1.25);
        // color continuo a lo largo de la paleta viridis, ordenado por población
        static const std::vector<int> kTurboPalette = {
    TColor::GetColor("#30123B"),  // 0
    TColor::GetColor("#3A2E8E"),  // 1
    TColor::GetColor("#4145AB"),  // 2
    TColor::GetColor("#4062D9"),  // 3
    TColor::GetColor("#4675ED"),  // 4
    TColor::GetColor("#3F8BF7"),  // 5
    TColor::GetColor("#39A2FC"),  // 6
    TColor::GetColor("#2AB9E5"),  // 7
    TColor::GetColor("#1BCFD4"),  // 8
    TColor::GetColor("#20DEBB"),  // 9
    TColor::GetColor("#26EDA2"),  // 10
    TColor::GetColor("#50F57E"),  // 11
    TColor::GetColor("#7BFC6B"),  // 12
    TColor::GetColor("#A6E54F"),  // 13
    TColor::GetColor("#D1E834"),  // 14
    TColor::GetColor("#F9BA38")   // 15
};
const int ci = kTurboPalette[k % kTurboPalette.size()];
        h->SetLineColor(ci);
        h->SetFillColorAlpha(ci, 0.45);
        h->SetLineWidth(2);
        if (h->GetEntries() > 0) {
            h->Draw("HIST");
        } else {
            TLatex t; t.SetNDC(); t.SetTextAlign(22); t.SetTextSize(0.12);
            t.SetTextColor(kGray + 1);
            t.DrawLatex(0.5, 0.5, "sin eventos");
        }
    }
    Save(c10, "10_dt_por_patron");

    // =========================================================================
    //  RESUMEN
    // =========================================================================
    printf("\n=======================================================================\n");
    printf("Ficheros guardados en %s :\n", CFG::kOut);
    const char* outs[11] = {"01_patrones", "02_eff_vs_dt", "03_mirror_test_backing",
                            "04_formas_dt", "05_eff_vs_amplitud", "06_eff_vs_asimetria",
                            "07_mapa_eff_amp0_amp1", "07b_diferencia_mapas",
                            "08_eff_vs_energia_neutron", "09_eff_vs_run",
                            "10_dt_por_patron"};
    for (auto o : outs) printf("   %s.pdf / .png\n", o);
    printf("=======================================================================\n");
    printf("\nLectura del canvas 03: si el cociente eps1/eps0 es plano y compatible\n"
           "con 1, el backing no sesga la deteccion. Una caida progresiva hacia\n"
           "dt > 0 es la firma del fragmento pesado frenado en el aluminio.\n\n");
}