// ---------------------------------------------------------------------------
// plot_costheta_slices.C
//
// Lee events_selection.root, aplica la seleccion nominal (getCuts) y dibuja
// las distribuciones de |cos_theta| normalizadas en rodajas de 20 MeV entre
// 500 y 800 MeV. Un canvas para uranio y otro para oro.
//
//   root -l -b -q plot_costheta_slices.C
//   root -l    'plot_costheta_slices.C("/otra/ruta.root", 500, 800, 20)'
// ---------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"

#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/config.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/cuts.h"

// ---------------------------------------------------------------------------
// Lector que acepta indistintamente ramas float o double.
// (amp0/amp1 aparecen como double en fillHistograms y float en fillRooDatasets)
// ---------------------------------------------------------------------------
struct DBranch {
    double d = 0.;
    float  f = 0.f;
    bool   isFloat = false;
    bool   ok = false;

    void connect(TTree* t, const char* name){
        TLeaf* l = t->GetLeaf(name);
        if(!l){
            std::cerr << "[ERROR] rama ausente: " << name << "\n";
            ok = false;
            return;
        }
        TString type = l->GetTypeName();          // "Float_t" / "Double_t"
        isFloat = type.Contains("Float");
        if(isFloat) t->SetBranchAddress(name, &f);
        else        t->SetBranchAddress(name, &d);
        ok = true;
    }
    inline double get() const { return isFloat ? double(f) : d; }
};

// ---------------------------------------------------------------------------
// Rellena un histograma de |cos_theta| por cada rodaja de energia.
// ---------------------------------------------------------------------------
static std::vector<TH1D*> buildSlices(TTree* tree,
                                      Sample sample,
                                      const char* tag,
                                      double e_lo, double e_hi, double e_step,
                                      int nbins_cos, bool use_abs_cos)
{
    const int nslices = int(std::lround((e_hi - e_lo) / e_step));

    const double cos_lo = use_abs_cos ? 0.0 : -1.0;
    const double cos_hi = 1.0;

    std::vector<TH1D*> h(nslices, nullptr);
    for(int s = 0; s < nslices; ++s){
        double a = e_lo + s * e_step;
        double b = a + e_step;
        h[s] = new TH1D(Form("h_%s_%d", tag, s),
                        Form("%s;%s;dN/d%s (norm.)",
                             tag,
                             use_abs_cos ? "|cos #theta|" : "cos #theta",
                             use_abs_cos ? "|cos #theta|" : "cos #theta"),
                        nbins_cos, cos_lo, cos_hi);
        h[s]->SetDirectory(nullptr);
        h[s]->Sumw2();
        h[s]->SetTitle(Form("%.0f #minus %.0f MeV", a, b));
    }

    DBranch tof0, tof1, amp0, amp1, en, cth, cthd;
    tree->SetBranchStatus("*", 1);
    tof0.connect(tree, "tof0");
    tof1.connect(tree, "tof1");
    amp0.connect(tree, "amp0");
    amp1.connect(tree, "amp1");
    en  .connect(tree, "neutron_energy");
    cth .connect(tree, "cos_theta");
    cthd.connect(tree, "cos_theta_det");

    if(!(tof0.ok && tof1.ok && amp0.ok && amp1.ok && en.ok && cth.ok && cthd.ok)){
        std::cerr << "[ERROR] faltan ramas en el arbol " << tree->GetName() << "\n";
        return h;
    }

    Long64_t nentries = tree->GetEntries();
    Long64_t nsel = 0;

    for(Long64_t i = 0; i < nentries; ++i){
        tree->GetEntry(i);

        const double E = en.get();
        if(E < e_lo || E >= e_hi) continue;

        const int s = int((E - e_lo) / e_step);
        if(s < 0 || s >= nslices) continue;

        const double ct  = cth.get();
        const double ctd = cthd.get();

        // cortes geometricos (identicos a fillHistograms)
        if(ctd < 0.0) continue;
        if(std::fabs(ctd) > 1.0 || std::fabs(ct) > 1.0) continue;

        // cortes de amplitud / ratio, dependientes de la energia
        EventCuts c = getCuts(sample, E);
        if(!passAmplitudeCut(float(amp0.get()), float(amp1.get()), c)) continue;

        // ventana de senal en tiempo de vuelo relativo
        const double dt = tof1.get() - tof0.get();
        if(dt < c.roi_min || dt > c.roi_max) continue;

        // por si acaso: el pico de U en la muestra de Au nunca solapa la ROI,
        // pero lo excluimos explicitamente
        if(inUraniumPeak(dt, c)) continue;

        h[s]->Fill(use_abs_cos ? std::fabs(ct) : ct);
        ++nsel;
    }

    std::cout << "\n=== " << tag << " ===  entradas: " << nentries
              << "   seleccionadas [" << e_lo << ", " << e_hi << ") MeV: "
              << nsel << "\n";
    for(int s = 0; s < nslices; ++s){
        std::cout << Form("  %4.0f - %4.0f MeV : %8.0f eventos\n",
                          e_lo + s*e_step, e_lo + (s+1)*e_step,
                          h[s]->Integral());
    }
    return h;
}

// ---------------------------------------------------------------------------
static void drawSlices(std::vector<TH1D*>& h,
                       const char* tag,
                       double e_lo, double e_step,
                       int min_counts)
{
    const int nslices = (int)h.size();

    TCanvas* c = new TCanvas(Form("c_%s", tag), Form("cos(theta) - %s", tag),
                             900, 700);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.05);
    c->SetTicks(1, 1);

    TLegend* leg = new TLegend(0.15, 0.62, 0.45, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.026);
    leg->SetNColumns(2);
    leg->SetHeader(tag);

    double ymax = 0.;
    int    ndrawn = 0;

    for(int s = 0; s < nslices; ++s){
        if(h[s]->Integral() < min_counts){
            std::cout << "  [skip] " << h[s]->GetTitle()
                      << " : estadistica insuficiente\n";
            continue;
        }
        // normalizacion a area unidad (con anchura de bin -> densidad)
        h[s]->Scale(1.0 / h[s]->Integral("width"));

        int col = TColor::GetColorPalette(
                      int(double(s) / std::max(1, nslices - 1) * 254));
        h[s]->SetLineColor(col);
        h[s]->SetMarkerColor(col);
        h[s]->SetLineWidth(2);
        h[s]->SetMarkerStyle(20);
        h[s]->SetMarkerSize(0.6);
        h[s]->SetStats(0);

        ymax = std::max(ymax, h[s]->GetMaximum() + h[s]->GetBinError(
                                  h[s]->GetMaximumBin()));
        leg->AddEntry(h[s], h[s]->GetTitle(), "lp");
        ++ndrawn;
    }

    if(ndrawn == 0){
        std::cerr << "[WARN] nada que dibujar para " << tag << "\n";
        return;
    }

    bool first = true;
    for(int s = 0; s < nslices; ++s){
        if(h[s]->Integral() <= 0) continue;
        h[s]->SetTitle("");
        h[s]->GetYaxis()->SetRangeUser(0., 1.35 * ymax);
        h[s]->GetYaxis()->SetTitleOffset(1.35);
        h[s]->Draw(first ? "E1 HIST P" : "E1 HIST P SAME");
        first = false;
    }
    leg->Draw();

    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.55, 0.92, Form("normalizado a area unidad, pasos de %.0f MeV",
                                  e_step));

    c->Update();
    c->SaveAs(Form("costheta_slices_%s.pdf", tag));
    c->SaveAs(Form("costheta_slices_%s.png", tag));
}

// ---------------------------------------------------------------------------
void angle_analysis(
    const char* filename =
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection.root",
    double e_lo   = 600.,
    double e_hi   = 1000.,
    double e_step = 50.,
    int    nbins_cos    = 20,
    bool   use_abs_cos  = true,
    int    min_counts   = 30)
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile* f = TFile::Open(filename, "READ");
    if(!f || f->IsZombie()){
        std::cerr << "[ERROR] no puedo abrir " << filename << "\n";
        return;
    }

    TTree* tU = (TTree*)f->Get("events_uranium");
    TTree* tA = (TTree*)f->Get("events_gold");
    if(!tU || !tA){
        std::cerr << "[ERROR] falta events_uranium o events_gold en el fichero\n";
        f->ls();
        return;
    }

    auto hU = buildSlices(tU, Sample::uranium, "uranium",
                          e_lo, e_hi, e_step, nbins_cos, use_abs_cos);
    auto hA = buildSlices(tA, Sample::gold,    "gold",
                          e_lo, e_hi, e_step, nbins_cos, use_abs_cos);

    drawSlices(hU, "uranium", e_lo, e_step, min_counts);
    drawSlices(hA, "gold",    e_lo, e_step, min_counts);

    f->Close();
}