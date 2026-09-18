#include "../include/types.h"
#include "../include/constants.h"
#include "../include/config.h"
#include "../include/utils.h"
#include "../include/cuts.h"
#include "../include/acceptance.h"
#include "../include/efficiency.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

static inline int cosBin(double c, int nbins)
{
    double a = std::fabs(c);
    if (a < 0. || a > 1.) return -1;
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

// etiqueta corta para nombres de ficheros/histogramas: 0.05 -> "thr050"
static std::string thrTag(double thr)
{
    return std::string(Form("thr%03d", (int)std::lround(thr * 1000.)));
}

// ── MC: árbol CoincTree, mismo corte de selección que mc_analysis(), pero
//    el umbral de amplitud se barre. Se hace UNA sola pasada por el árbol y
//    cada evento se acumula en todos los umbrales que pasa (los conjuntos
//    son anidados: amp_min >= thr) ──────────────────────────────────────────
static bool fillCountsMC(TTree* t, std::vector<McCounts>& counts,
                          const std::vector<double>& thresholds,
                          const std::vector<double>& energy_bins_eff, int nbins_eff)
{
    Double_t neutronE, cos_theta_det, cos_theta;
    Int_t    target_i, Z0, Z1, A0, A1;
    Bool_t   hasTrue_i, hasTrue_j;
    Int_t    ppac0;
    Double_t amp0_c1, amp1_c1, amp0_c2, amp1_c2;

    t->SetBranchAddress("neutronE",      &neutronE);
    t->SetBranchAddress("cos_theta_det", &cos_theta_det);
    t->SetBranchAddress("cos_theta",     &cos_theta);
    t->SetBranchAddress("target_i",      &target_i);
    t->SetBranchAddress("Z0", &Z0);  t->SetBranchAddress("Z1", &Z1);
    t->SetBranchAddress("A0", &A0);  t->SetBranchAddress("A1", &A1);
    t->SetBranchAddress("hasTrue_i", &hasTrue_i);
    t->SetBranchAddress("hasTrue_j", &hasTrue_j);
    t->SetBranchAddress("ppac0",     &ppac0);
    t->SetBranchAddress("amp0_c1", &amp0_c1);
    t->SetBranchAddress("amp1_c1", &amp1_c1);
    t->SetBranchAddress("amp0_c2", &amp0_c2);
    t->SetBranchAddress("amp1_c2", &amp1_c2);

    const size_t nthr = thresholds.size();
    std::vector<Long64_t> nsel(nthr, 0);

    for (Long64_t k = 0, n = t->GetEntries(); k < n; ++k) {
        t->GetEntry(k);

        if (target_i != 2)                continue;   // TARGET_U
        if (!hasTrue_i || !hasTrue_j)     continue;
        if (Z0 <= 2 || Z1 <= 2)           continue;
        if (Z0 + Z1 < 80 || Z0 + Z1 > 92) continue;
        if (A0 + A1 < 200)                continue;
        if (ppac0 != 8)                   continue;

        int d = cosBin(cos_theta_det, nbins_det);
        if (d < 0) continue;

        int e_eff = findBin(energy_bins_eff, neutronE);
        if (e_eff < 0 || e_eff >= nbins_eff) continue;

        // OJO: hay que usar el bin de beam REAL (cos_theta), no un j=0
        // sintético -- si no, computeEfficiency solo usa dOmega_eff[0][i]
        // como normalización, en vez de la aceptancia real de cada evento.
        int j = cosBin(cos_theta, nbins_beam);
        if (j < 0) continue;

        // el corte original era amp*_c* >= 0.05 en las cuatro amplitudes:
        // equivale a exigir min(amp) >= umbral
        double amp_min = std::min({amp0_c1, amp1_c1, amp0_c2, amp1_c2});

        for (size_t s = 0; s < nthr; ++s) {
            if (amp_min < thresholds[s]) continue;
            counts[s].n[e_eff][j][d] += 1.;
            ++nsel[s];
        }
    }

    Long64_t tot = 0;
    for (size_t s = 0; s < nthr; ++s) {
        std::cout << "  seleccionados (MC, amp > " << thresholds[s] << "): "
                  << nsel[s] << "\n";
        tot += nsel[s];
    }
    return tot > 0;
}

// ── DATOS: árboles events_uranium / events_gold, ramas y cortes reales de
//    datos (tof0/tof1/amp0/amp1/x0/x1/cos_theta/cos_theta_det), igual que
//    fillHistograms() -- solo se acumulan los eventos dentro de la ROI
//    (señal), que es el equivalente a "hasTrue_i && hasTrue_j" en MC ───────
static bool fillCountsData(TTree* t, McCounts& counts, Sample sample,
                            const std::vector<double>& energy_bins_eff, int nbins_eff)
{
    Double_t cos_theta_det, cos_theta;
    Double_t tof1, tof0, neutron_energy;
    Double_t  amp0, amp1;
    Double_t x0, x1;

    t->SetBranchAddress("tof1",           &tof1);
    t->SetBranchAddress("tof0",           &tof0);
    t->SetBranchAddress("amp0",           &amp0);
    t->SetBranchAddress("amp1",           &amp1);
    t->SetBranchAddress("x1",             &x1);
    t->SetBranchAddress("x0",             &x0);
    t->SetBranchAddress("neutron_energy", &neutron_energy);
    t->SetBranchAddress("cos_theta",      &cos_theta);
    t->SetBranchAddress("cos_theta_det",  &cos_theta_det);

    Long64_t nsel = 0;
    for (Long64_t k = 0, n = t->GetEntries(); k < n; ++k) {
        t->GetEntry(k);

        if (neutron_energy > 1000)                      continue;
        if (cos_theta_det < 0.0)                        continue;
        if (std::fabs(cos_theta_det) > 1 || std::fabs(cos_theta) > 1) continue;

        EventCuts c = getCuts(sample, neutron_energy);
        if (!passAmplitudeCut(amp0, amp1, c))            continue;

        double dt = tof1 - tof0;
        if (dt < c.roi_min || dt > c.roi_max)            continue;   // solo señal (ROI)

        int e_eff = findBin(energy_bins_eff, neutron_energy);
        if (e_eff < 0 || e_eff >= nbins_eff)             continue;

        int d = cosBin(cos_theta_det, nbins_det);
        if (d < 0) continue;

        // bin de beam real, igual que en MC -- necesario para que
        // computeEfficiency normalice con la aceptancia correcta por evento
        int j = cosBin(cos_theta, nbins_beam);
        if (j < 0) continue;

        ++nsel;
        counts.n[e_eff][j][d] += 1.;
    }
    std::cout << "  seleccionados (DATA, ROI): " << nsel << "\n";
    return nsel > 0;
}

// ── eficiencia (método con overlap) para un McCounts ya lleno ─────────────
static void computeOverlapEff(McCounts& counts, const Vec2D& acceptance, int nbins_eff,
                               std::vector<EfficiencyResult>& eff)
{
    poissonErrors(counts);
    eff.resize(nbins_eff);
    for (int e = 0; e < nbins_eff; ++e)
        eff[e] = computeEfficiency(nbins_det - 1, nbins_beam, nbins_det,
                                    counts.n, counts.u_n, acceptance, e);
}

// ── grafica datos vs MC para cada umbral, y el cociente MC/datos ──────────
static void plotMCvsData(const std::vector<std::vector<EfficiencyResult>>& eff_mc,
                          const std::vector<double>& thresholds,
                          const std::vector<EfficiencyResult>& eff_data,
                          int nbins_eff, const std::string& outdir)
{
    const int colors[5] = { kRed + 1, kOrange + 7, kGreen + 2, kAzure + 2, kViolet + 1 };
    const size_t nthr = thresholds.size();

    for (int e = 0; e < nbins_eff; ++e) {

        // --- datos ---
        std::vector<double> xD, exD, yD, eyD;
        for (int i = 0; i < nbins_det; ++i) {
            double v = eff_data[e].eps[i];
            if (v <= 0.) continue;
            xD.push_back((i + 0.5) / nbins_det);
            exD.push_back(0.5 / nbins_det);
            yD.push_back(v);
            eyD.push_back(eff_data[e].u_eps[i]);
        }
        if (xD.empty()) continue;

        double ymin = 1e300, ymax = -1e300;
        for (size_t k = 0; k < yD.size(); ++k) {
            ymin = std::min(ymin, yD[k] - eyD[k]);
            ymax = std::max(ymax, yD[k] + eyD[k]);
        }

        // --- MC, un grafo por umbral ---
        std::vector<TGraphErrors*> gMC(nthr, nullptr);
        for (size_t s = 0; s < nthr; ++s) {
            if ((int)eff_mc[s].size() < nbins_eff) continue;
            std::vector<double> x, ex, y, ey;
            for (int i = 0; i < nbins_det; ++i) {
                double v = eff_mc[s][e].eps[i];
                if (v <= 0.) continue;
                x.push_back((i + 0.5) / nbins_det);
                ex.push_back(0.5 / nbins_det);
                y.push_back(v);
                ey.push_back(eff_mc[s][e].u_eps[i]);
            }
            if (x.empty()) continue;
            for (size_t k = 0; k < y.size(); ++k) {
                ymin = std::min(ymin, y[k] - ey[k]);
                ymax = std::max(ymax, y[k] + ey[k]);
            }
            gMC[s] = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
            gMC[s]->SetName(Form("eff_mc_%s_ebin%d", thrTag(thresholds[s]).c_str(), e));
            gMC[s]->SetMarkerStyle(24);
            gMC[s]->SetMarkerColor(colors[s % 5]);
            gMC[s]->SetLineColor(colors[s % 5]);
        }

        TGraphErrors* gD = new TGraphErrors(xD.size(), xD.data(), yD.data(), exD.data(), eyD.data());
        gD->SetName(Form("eff_data_ebin%d", e));
        gD->SetMarkerStyle(20);
        gD->SetMarkerColor(kBlack);
        gD->SetLineColor(kBlack);
        gD->SetLineWidth(2);

        // rango Y conjunto -- si no, al dibujar gD primero con "ALP" el eje se
        // autoescala solo a los datos y las curvas MC pueden quedar cortadas
        double pad = 0.1 * (ymax - ymin);
        gD->SetMinimum(ymin - pad);
        gD->SetMaximum(ymax + pad);

        TCanvas* cv = new TCanvas(Form("cv_eff_ebin%d", e),
                                   "efficiency MC thresholds vs data", 900, 600);
        gD->SetTitle(Form("Eficiencia (overlap) ebin %d;|cos#theta_{det}|;#varepsilon (rel.)", e));
        gD->Draw("ALP");
        for (size_t s = 0; s < nthr; ++s)
            if (gMC[s]) gMC[s]->Draw("LP same");

        TLegend* leg = new TLegend(0.60, 0.62, 0.89, 0.89);
        leg->AddEntry(gD, "datos", "lp");
        for (size_t s = 0; s < nthr; ++s)
            if (gMC[s]) leg->AddEntry(gMC[s], Form("MC amp > %.2f", thresholds[s]), "lp");
        leg->Draw();
        cv->SaveAs((outdir + "eff_mc_thresholds_vs_data_ebin" + std::to_string(e) + ".pdf").c_str());

        // --- cociente MC / datos ---
        std::vector<TGraphErrors*> gR(nthr, nullptr);
        double rmin = 1e300, rmax = -1e300;
        for (size_t s = 0; s < nthr; ++s) {
            if ((int)eff_mc[s].size() < nbins_eff) continue;
            std::vector<double> x, ex, y, ey;
            for (int i = 0; i < nbins_det; ++i) {
                double vD = eff_data[e].eps[i],  uD = eff_data[e].u_eps[i];
                double vM = eff_mc[s][e].eps[i], uM = eff_mc[s][e].u_eps[i];
                if (vD <= 0. || vM <= 0.) continue;
                double r = vM / vD;
                x.push_back((i + 0.5) / nbins_det);
                ex.push_back(0.5 / nbins_det);
                y.push_back(r);
                ey.push_back(r * std::sqrt(std::pow(uM / vM, 2) + std::pow(uD / vD, 2)));
            }
            if (x.empty()) continue;
            for (size_t k = 0; k < y.size(); ++k) {
                rmin = std::min(rmin, y[k] - ey[k]);
                rmax = std::max(rmax, y[k] + ey[k]);
            }
            gR[s] = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
            gR[s]->SetName(Form("eff_ratio_mc_%s_over_data_ebin%d", thrTag(thresholds[s]).c_str(), e));
            gR[s]->SetMarkerStyle(20);
            gR[s]->SetMarkerColor(colors[s % 5]);
            gR[s]->SetLineColor(colors[s % 5]);

            double sum = 0.;
            for (double v : y) sum += (v - 1.);
            std::cout << "ebin " << e << ", MC amp > " << thresholds[s]
                      << ": diferencia media (MC/datos - 1) = "
                      << sum / y.size() * 100. << " %\n";
        }

        TGraphErrors* gFirst = nullptr;
        for (size_t s = 0; s < nthr && !gFirst; ++s) gFirst = gR[s];
        if (!gFirst) continue;

        double rpad = 0.1 * (rmax - rmin);
        gFirst->SetMinimum(rmin - rpad);
        gFirst->SetMaximum(rmax + rpad);
        gFirst->SetTitle(Form("MC/Data ebin %d;|cos#theta_{det}|;MC / Data", e));

        TCanvas* cv2 = new TCanvas(Form("cv_ratio_ebin%d", e), "efficiency ratio", 900, 500);
        gFirst->Draw("ALP");
        for (size_t s = 0; s < nthr; ++s)
            if (gR[s] && gR[s] != gFirst) gR[s]->Draw("LP same");

        TLine* one = new TLine(0., 1., 1., 1.);
        one->SetLineStyle(2);
        one->SetLineColor(kGray + 2);
        one->Draw("same");

        TLegend* leg2 = new TLegend(0.60, 0.68, 0.89, 0.89);
        for (size_t s = 0; s < nthr; ++s)
            if (gR[s]) leg2->AddEntry(gR[s], Form("MC amp > %.2f", thresholds[s]), "lp");
        leg2->Draw();
        cv2->SaveAs((outdir + "eff_ratio_mc_over_data_ebin" + std::to_string(e) + ".pdf").c_str());
    }
}

// ── guarda un TH1D con eps/u_eps por bin de energía, si el vector no está
//    vacío (evita out-of-bounds si el fichero de origen no abrió) ─────────
static void saveEff(const std::vector<EfficiencyResult>& eff, int nbins_eff,
                     const std::string& tag)
{
    if ((int)eff.size() < nbins_eff) {
        std::cerr << "  [saveEff] " << tag << ": sin datos, se omite.\n";
        return;
    }
    for (int e = 0; e < nbins_eff; ++e) {
        TH1D* h = new TH1D(Form("heff_%s_ebin%d", tag.c_str(), e), "", nbins_det, 0, 1);
        for (int i = 0; i < nbins_det; ++i) {
            h->SetBinContent(i + 1, eff[e].eps[i]);
            h->SetBinError  (i + 1, eff[e].u_eps[i]);
        }
        h->Write();
    }
}

// ---------------------------------------------------------------------------
void compare_efficiencies()
{
    const std::string outdir = "/Users/nico/Desktop/Tese/Analysis/montecarlo/output/";

    const std::string mc_file   = "/Users/nico/Desktop/Tese/Analysis/montecarlo/data/mc_setup.root";
    const std::string mc_tree   = "CoincTree";

    const std::string data_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection_exp.root";
    const std::string data_tree = "events_uranium";   // target 2 = uranio

    const std::string acceptance_file =
        "/Users/nico/Desktop/Tese/Analysis/acceptance_coincidence.csv";

    const std::vector<double> energy_bins_eff = {1, 10, 100, 600, 1000};
    const int nbins_eff = (int)energy_bins_eff.size() - 1;

    // umbrales de amplitud aplicados al MC
    const std::vector<double> thresholds = {0.29, 0.3, 0.31, 0.32, 0.33, 0.34};
    const size_t nthr = thresholds.size();

    Vec2D acceptance, dOmega_fine;
    if (!loadAcceptanceCSV(acceptance_file, dOmega_fine)) {
        std::cerr << "Failed to load acceptance CSV\n";
        return;
    }
    acceptance = rebin(dOmega_fine);

    // ── MC: un juego de eficiencias por umbral ──────────────────────────
    std::vector<std::vector<EfficiencyResult>> eff_mc(nthr);
    {
        std::cout << "== MC (" << mc_file << ") ==\n";
        TFile* fin = TFile::Open(mc_file.c_str());
        if (!fin || fin->IsZombie()) {
            std::cerr << "cannot open " << mc_file << " -- se omite MC\n";
        } else {
            TTree* t = (TTree*)fin->Get(mc_tree.c_str());
            if (!t) {
                std::cerr << mc_tree << " not found in " << mc_file << " -- se omite MC\n";
            } else {
                std::vector<McCounts> counts;
                counts.reserve(nthr);
                for (size_t s = 0; s < nthr; ++s) counts.push_back(makeCounts(nbins_eff));

                fillCountsMC(t, counts, thresholds, energy_bins_eff, nbins_eff);
                for (size_t s = 0; s < nthr; ++s)
                    computeOverlapEff(counts[s], acceptance, nbins_eff, eff_mc[s]);
            }
            fin->Close();
        }
    }

    // ── Datos (un solo juego, cortes de datos sin cambios) ──────────────
    std::vector<EfficiencyResult> eff_data;
    {
        std::cout << "== DATA (" << data_file << ", " << data_tree << ") ==\n";
        TFile* fin = TFile::Open(data_file.c_str());
        if (!fin || fin->IsZombie()) {
            std::cerr << "cannot open " << data_file << " -- se omite DATA\n";
        } else {
            TTree* t = (TTree*)fin->Get(data_tree.c_str());
            if (!t) {
                std::cerr << data_tree << " not found in " << data_file << " -- se omite DATA\n";
            } else {
                McCounts counts = makeCounts(nbins_eff);
                fillCountsData(t, counts, Sample::uranium, energy_bins_eff, nbins_eff);
                computeOverlapEff(counts, acceptance, nbins_eff, eff_data);
            }
            fin->Close();
        }
    }

    // ── comparación MC (por umbral) vs datos ────────────────────────────
    bool have_mc = false;
    for (size_t s = 0; s < nthr; ++s) if (!eff_mc[s].empty()) have_mc = true;
    if (have_mc && !eff_data.empty())
        plotMCvsData(eff_mc, thresholds, eff_data, nbins_eff, outdir);
    else
        std::cerr << "No hay MC y/o datos suficientes para comparar.\n";

    // ── guardar todo lo que sí se pudo calcular ────────────────────────────
    TFile* fout = TFile::Open((outdir + "output_efficiency_threshold_comparison.root").c_str(), "RECREATE");
    for (size_t s = 0; s < nthr; ++s)
        saveEff(eff_mc[s], nbins_eff, "mc_" + thrTag(thresholds[s]));
    saveEff(eff_data, nbins_eff, "data");
    fout->Close();

    std::cout << "Guardado: " << outdir << "output_efficiency_threshold_comparison.root\n";
}