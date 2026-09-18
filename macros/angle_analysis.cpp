// ---------------------------------------------------------------------------
// angle_analysis.C
//
// Lee events_selection.root, aplica la seleccion nominal (getCuts) y dibuja
// las distribuciones de |cos_theta| normalizadas en rodajas de energia,
// separadas en grupos de N runs consecutivos (por defecto 30).
//
// Salidas por muestra (uranium / gold):
//   costheta_groups_<tag>_slice<N>.pdf   grupos superpuestos, una rodaja
//   costheta_slices_<tag>_g<N>.pdf       rodajas superpuestas, un grupo  (opt)
//   meancos_stability_<tag>.pdf          <|cos theta|> vs grupo, todas las rodajas
//
//   root -l -b -q angle_analysis.C
//   root -l    'angle_analysis.C("/otra/ruta.root", 100, 800, 100)'
// ---------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <algorithm>
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
#include "TString.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

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
// Equivalente para ramas enteras (RunNumber puede ser Int_t o Long64_t).
// ---------------------------------------------------------------------------
struct IBranch {
    Int_t    i32 = 0;
    Long64_t i64 = 0;
    bool     isLong = false;
    bool     ok = false;

    void connect(TTree* t, const char* name){
        TLeaf* l = t->GetLeaf(name);
        if(!l){
            std::cerr << "[ERROR] rama ausente: " << name << "\n";
            ok = false;
            return;
        }
        TString type = l->GetTypeName();          // "Int_t" / "Long64_t" / "UInt_t"
        isLong = type.Contains("Long");
        if(isLong) t->SetBranchAddress(name, &i64);
        else       t->SetBranchAddress(name, &i32);
        ok = true;
    }
    inline long long get() const { return isLong ? (long long)i64 : (long long)i32; }
};

// ---------------------------------------------------------------------------
// Agrupacion por run number
// ---------------------------------------------------------------------------

// Recorre un arbol solo con la rama RunNumber activa y acumula los runs
// presentes. Barato: ROOT solo descomprime esa rama.
static bool collectRuns(TTree* tree, std::set<int>& runset)
{
    IBranch run;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("RunNumber", 1);
    run.connect(tree, "RunNumber");

    if(!run.ok){
        tree->SetBranchStatus("*", 1);
        tree->ResetBranchAddresses();
        return false;
    }

    for(Long64_t i = 0; i < tree->GetEntries(); ++i){
        tree->GetEntry(i);
        runset.insert((int)run.get());
    }

    // imprescindible: `run` es local y muere al salir, pero el arbol se vuelve
    // a leer despues en buildSlicesByGroup
    tree->SetBranchStatus("*", 1);
    tree->ResetBranchAddresses();
    return true;
}

// Trocea la lista ORDENADA de runs realmente presentes en bloques de
// group_size. Usar run/group_size en su lugar daria grupos de estadistica
// muy dispar cuando hay huecos (runs malos, paradas de haz).
static std::vector<std::vector<int>> makeRunGroups(const std::set<int>& runset,
                                                   int group_size)
{
    std::vector<int> runs(runset.begin(), runset.end());
    std::vector<std::vector<int>> groups;

    for(size_t i = 0; i < runs.size(); i += (size_t)group_size)
        groups.emplace_back(runs.begin() + i,
                            runs.begin() + std::min(i + (size_t)group_size,
                                                    runs.size()));

    std::cout << "[info] " << runs.size() << " runs -> "
              << groups.size() << " grupos de hasta " << group_size << "\n";
    for(size_t g = 0; g < groups.size(); ++g)
        std::cout << Form("       grupo %2zu : runs %d - %d  (%zu runs)\n",
                          g, groups[g].front(), groups[g].back(),
                          groups[g].size());
    return groups;
}

static std::map<int,int> runToGroup(const std::vector<std::vector<int>>& groups)
{
    std::map<int,int> m;
    for(size_t g = 0; g < groups.size(); ++g)
        for(int r : groups[g]) m[r] = (int)g;
    return m;
}

static std::string groupLabel(const std::vector<int>& g)
{
    return std::string(Form("runs %d#minus%d", g.front(), g.back()));
}

// ---------------------------------------------------------------------------
// Rellena un histograma de |cos_theta| por cada (grupo de runs, rodaja de E).
// Una sola pasada sobre el arbol para todos los grupos.
// ---------------------------------------------------------------------------
static std::vector<std::vector<TH1D*>> buildSlicesByGroup(
        TTree* tree,
        Sample sample,
        const char* tag,
        const std::map<int,int>& run_to_group,
        int ngroups,
        double e_lo, double e_hi, double e_step,
        int nbins_cos, bool use_abs_cos)
{
    const int nslices = int(std::lround((e_hi - e_lo) / e_step));

    const double cos_lo = use_abs_cos ? 0.0 : -1.0;
    const double cos_hi = 1.0;

    std::vector<std::vector<TH1D*>> h(ngroups,
                                      std::vector<TH1D*>(nslices, nullptr));
    for(int g = 0; g < ngroups; ++g){
        for(int s = 0; s < nslices; ++s){
            double a = e_lo + s * e_step;
            double b = a + e_step;
            h[g][s] = new TH1D(Form("h_%s_g%d_s%d", tag, g, s),
                               Form("%s;%s;dN/d%s (norm.)",
                                    tag,
                                    use_abs_cos ? "|cos #theta|" : "cos #theta",
                                    use_abs_cos ? "|cos #theta|" : "cos #theta"),
                               nbins_cos, cos_lo, cos_hi);
            h[g][s]->SetDirectory(nullptr);
            h[g][s]->Sumw2();
            h[g][s]->SetTitle(Form("%.0f #minus %.0f MeV", a, b));
        }
    }

    DBranch tof0, tof1, amp0, amp1, en, cth, cthd;
    IBranch run;
    tree->SetBranchStatus("*", 1);
    tof0.connect(tree, "tof0");
    tof1.connect(tree, "tof1");
    amp0.connect(tree, "amp0");
    amp1.connect(tree, "amp1");
    en  .connect(tree, "neutron_energy");
    cth .connect(tree, "cos_theta");
    cthd.connect(tree, "cos_theta_det");
    run .connect(tree, "RunNumber");

    if(!(tof0.ok && tof1.ok && amp0.ok && amp1.ok &&
         en.ok && cth.ok && cthd.ok && run.ok)){
        std::cerr << "[ERROR] faltan ramas en el arbol " << tree->GetName() << "\n";
        return h;
    }

    Long64_t nentries = tree->GetEntries();
    Long64_t nsel_tot = 0, norphan = 0;
    std::vector<Long64_t> nsel(ngroups, 0);

    for(Long64_t i = 0; i < nentries; ++i){
        tree->GetEntry(i);

        // corte mas barato primero: descarta el evento antes de tocar cortes
        auto it = run_to_group.find((int)run.get());
        if(it == run_to_group.end()){ ++norphan; continue; }
        const int g = it->second;

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

        h[g][s]->Fill(use_abs_cos ? std::fabs(ct) : ct);
        ++nsel[g];
        ++nsel_tot;
    }

    std::cout << "\n=== " << tag << " ===  entradas: " << nentries
              << "   seleccionadas [" << e_lo << ", " << e_hi << ") MeV: "
              << nsel_tot << "\n";
    if(norphan)
        std::cout << "[WARN] " << norphan
                  << " eventos con RunNumber fuera del mapa de grupos\n";
    for(int g = 0; g < ngroups; ++g){
        std::cout << Form("  grupo %2d : %8lld eventos  |", g, (long long)nsel[g]);
        for(int s = 0; s < nslices; ++s)
            std::cout << Form(" %6.0f", h[g][s]->Integral());
        std::cout << "\n";
    }
    return h;
}

// ---------------------------------------------------------------------------
// Una rodaja de energia, todos los grupos superpuestos.
// Es el grafico util para ver deriva entre periodos de toma de datos.
// ---------------------------------------------------------------------------
static void drawGroupsForSlice(std::vector<std::vector<TH1D*>>& h,
                               const std::vector<std::vector<int>>& groups,
                               const char* tag,
                               int s,
                               double e_lo, double e_step,
                               int min_counts)
{
    const int ngroups = (int)h.size();

    TCanvas* c = new TCanvas(Form("c_%s_s%d", tag, s),
                             Form("cos(theta) - %s - slice %d", tag, s),
                             900, 700);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.05);
    c->SetTicks(1, 1);

    TLegend* leg = new TLegend(0.15, 0.60, 0.48, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.024);
    leg->SetNColumns(2);
    leg->SetHeader(Form("%s   %.0f #minus %.0f MeV",
                        tag, e_lo + s*e_step, e_lo + (s+1)*e_step));

    // Solo se dibuja lo que pasa min_counts, y solo eso se normaliza.
    std::vector<TH1D*> keep;
    double ymax = 0.;

    for(int g = 0; g < ngroups; ++g){
        if(h[g][s]->Integral() < min_counts){
            std::cout << "  [skip] " << tag << " slice " << s
                      << " grupo " << g << " : estadistica insuficiente ("
                      << h[g][s]->Integral() << ")\n";
            continue;
        }
        // normalizacion a area unidad (con anchura de bin -> densidad)
        h[g][s]->Scale(1.0 / h[g][s]->Integral("width"));

        int col = TColor::GetColorPalette(
                      int(double(g) / std::max(1, ngroups - 1) * 254));
        h[g][s]->SetLineColor(col);
        h[g][s]->SetMarkerColor(col);
        h[g][s]->SetLineWidth(2);
        h[g][s]->SetMarkerStyle(20);
        h[g][s]->SetMarkerSize(0.6);
        h[g][s]->SetStats(0);
        h[g][s]->SetTitle("");

        ymax = std::max(ymax, h[g][s]->GetMaximum()
                            + h[g][s]->GetBinError(h[g][s]->GetMaximumBin()));
        leg->AddEntry(h[g][s], groupLabel(groups[g]).c_str(), "lp");
        keep.push_back(h[g][s]);
    }

    if(keep.empty()){
        std::cerr << "[WARN] nada que dibujar para " << tag
                  << " slice " << s << "\n";
        return;
    }

    for(size_t k = 0; k < keep.size(); ++k){
        keep[k]->GetYaxis()->SetRangeUser(0., 1.35 * ymax);
        keep[k]->GetYaxis()->SetTitleOffset(1.35);
        keep[k]->Draw(k == 0 ? "E1 HIST P" : "E1 HIST P SAME");
    }
    leg->Draw();

    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.028);
    tx.DrawLatex(0.55, 0.92, "normalizado a area unidad");

    c->Update();
    c->SaveAs(Form("costheta_groups_%s_slice%d.pdf", tag, s));
}

// ---------------------------------------------------------------------------
// Un grupo, todas las rodajas de energia superpuestas (el grafico original).
// ---------------------------------------------------------------------------
static void drawSlicesForGroup(std::vector<std::vector<TH1D*>>& h,
                               const std::vector<std::vector<int>>& groups,
                               const char* tag,
                               int g,
                               double e_lo, double e_step,
                               int min_counts)
{
    const int nslices = (int)h[g].size();

    TCanvas* c = new TCanvas(Form("c_%s_g%d", tag, g),
                             Form("cos(theta) - %s - group %d", tag, g),
                             900, 700);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.05);
    c->SetTicks(1, 1);

    TLegend* leg = new TLegend(0.15, 0.62, 0.45, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.026);
    leg->SetNColumns(2);
    leg->SetHeader(Form("%s  %s", tag, groupLabel(groups[g]).c_str()));

    std::vector<TH1D*> keep;
    double ymax = 0.;

    for(int s = 0; s < nslices; ++s){
        if(h[g][s]->Integral() < min_counts) continue;
        // puede haberse normalizado ya en drawGroupsForSlice; Integral("width")
        // sobre un histograma ya normalizado vale 1 y la operacion es inocua
        h[g][s]->Scale(1.0 / h[g][s]->Integral("width"));

        int col = TColor::GetColorPalette(
                      int(double(s) / std::max(1, nslices - 1) * 254));
        h[g][s]->SetLineColor(col);
        h[g][s]->SetMarkerColor(col);
        h[g][s]->SetLineWidth(2);
        h[g][s]->SetMarkerStyle(20);
        h[g][s]->SetMarkerSize(0.6);
        h[g][s]->SetStats(0);

        ymax = std::max(ymax, h[g][s]->GetMaximum()
                            + h[g][s]->GetBinError(h[g][s]->GetMaximumBin()));
        leg->AddEntry(h[g][s], Form("%.0f #minus %.0f MeV",
                                    e_lo + s*e_step, e_lo + (s+1)*e_step), "lp");
        keep.push_back(h[g][s]);
    }

    if(keep.empty()){
        std::cerr << "[WARN] nada que dibujar para " << tag
                  << " grupo " << g << "\n";
        return;
    }

    for(size_t k = 0; k < keep.size(); ++k){
        keep[k]->SetTitle("");
        keep[k]->GetYaxis()->SetRangeUser(0., 1.35 * ymax);
        keep[k]->GetYaxis()->SetTitleOffset(1.35);
        keep[k]->Draw(k == 0 ? "E1 HIST P" : "E1 HIST P SAME");
    }
    leg->Draw();

    c->Update();
    c->SaveAs(Form("costheta_slices_%s_g%d.pdf", tag, g));
}

// ---------------------------------------------------------------------------
// Resumen: <|cos theta|> vs indice de grupo, una curva por rodaja de energia.
// Un grupo que se desvia sistematicamente salta a la vista aqui mucho antes
// que comparando histogramas superpuestos.
// ---------------------------------------------------------------------------
static void drawMeanStability(std::vector<std::vector<TH1D*>>& h,
                              const std::vector<std::vector<int>>& groups,
                              const char* tag,
                              double e_lo, double e_step,
                              int min_counts)
{
    const int ngroups = (int)h.size();
    if(ngroups == 0) return;
    const int nslices = (int)h[0].size();

    TCanvas* c = new TCanvas(Form("c_mean_%s", tag),
                             Form("<|cos theta|> vs grupo - %s", tag), 950, 650);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.05);
    c->SetTicks(1, 1);

    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = new TLegend(0.72, 0.60, 0.94, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.026);
    leg->SetHeader(tag);

    for(int s = 0; s < nslices; ++s){
        std::vector<double> x, y, ex, ey;
        for(int g = 0; g < ngroups; ++g){
            // GetEntries no cambia con Scale(); Integral si. Usamos entradas.
            if(h[g][s]->GetEntries() < min_counts) continue;
            x .push_back(g);
            y .push_back(h[g][s]->GetMean());
            ex.push_back(0.0);
            ey.push_back(h[g][s]->GetMeanError());
        }
        if(x.size() < 2) continue;

        TGraphErrors* gr = new TGraphErrors((int)x.size(),
                                            x.data(), y.data(),
                                            ex.data(), ey.data());
        int col = TColor::GetColorPalette(
                      int(double(s) / std::max(1, nslices - 1) * 254));
        gr->SetLineColor(col);
        gr->SetMarkerColor(col);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.9);
        gr->SetLineWidth(2);
        mg->Add(gr, "LP");
        leg->AddEntry(gr, Form("%.0f #minus %.0f MeV",
                               e_lo + s*e_step, e_lo + (s+1)*e_step), "lp");
    }

    if(mg->GetListOfGraphs() == nullptr ||
       mg->GetListOfGraphs()->GetSize() == 0){
        std::cerr << "[WARN] sin puntos para la estabilidad de " << tag << "\n";
        return;
    }

    mg->Draw("A");
    mg->GetXaxis()->SetTitle("indice de grupo de runs");
    mg->GetYaxis()->SetTitle("#LT|cos #theta|#GT");
    mg->GetYaxis()->SetTitleOffset(1.35);
    leg->Draw();

    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.026);
    tx.DrawLatex(0.15, 0.92, Form("grupo 0 = runs %d#minus%d,  grupo %d = runs %d#minus%d",
                                  groups.front().front(), groups.front().back(),
                                  ngroups - 1,
                                  groups.back().front(), groups.back().back()));

    c->Update();
    c->SaveAs(Form("meancos_stability_%s.pdf", tag));
}

// ---------------------------------------------------------------------------
void angle_analysis(
    const char* filename =
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection_exp.root",
    double e_lo   = 100.,
    double e_hi   = 800.,
    double e_step = 100.,
    int    nbins_cos    = 20,
    bool   use_abs_cos  = true,
    int    min_counts   = 5,     // por grupo: mucho mas bajo que el valor global
    int    group_size   = 30,
    bool   per_group_slices = false)
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

    // --- grupos de runs, construidos a partir de LOS DOS arboles ------------
    // si el mapa se construyera solo con uranio, los runs presentes unicamente
    // en oro se descartarian en silencio
    std::set<int> runset;
    bool okU = collectRuns(tU, runset);
    bool okA = collectRuns(tA, runset);
    if(!okU || !okA){
        std::cerr << "[ERROR] no puedo leer RunNumber en alguno de los arboles\n";
        f->Close();
        return;
    }
    if(runset.empty()){
        std::cerr << "[ERROR] no hay runs en el fichero\n";
        f->Close();
        return;
    }

    auto groups = makeRunGroups(runset, group_size);
    auto r2g    = runToGroup(groups);
    const int ngroups = (int)groups.size();
    const int nslices = int(std::lround((e_hi - e_lo) / e_step));

    // --- llenado (una pasada por arbol) -------------------------------------
    auto hU = buildSlicesByGroup(tU, Sample::uranium, "uranium",
                                 r2g, ngroups,
                                 e_lo, e_hi, e_step, nbins_cos, use_abs_cos);
    auto hA = buildSlicesByGroup(tA, Sample::gold, "gold",
                                 r2g, ngroups,
                                 e_lo, e_hi, e_step, nbins_cos, use_abs_cos);

    // --- estabilidad: primero, porque usa las medias sin normalizar ---------
    drawMeanStability(hU, groups, "uranium", e_lo, e_step, min_counts);
    drawMeanStability(hA, groups, "gold",    e_lo, e_step, min_counts);

    // --- una rodaja, todos los grupos ---------------------------------------
    for(int s = 0; s < nslices; ++s){
        drawGroupsForSlice(hU, groups, "uranium", s, e_lo, e_step, min_counts);
        drawGroupsForSlice(hA, groups, "gold",    s, e_lo, e_step, min_counts);
    }

    // --- un grupo, todas las rodajas (opcional: ngroups x 2 ficheros) -------
    if(per_group_slices){
        for(int g = 0; g < ngroups; ++g){
            drawSlicesForGroup(hU, groups, "uranium", g, e_lo, e_step, min_counts);
            drawSlicesForGroup(hA, groups, "gold",    g, e_lo, e_step, min_counts);
        }
    }

    f->Close();
}