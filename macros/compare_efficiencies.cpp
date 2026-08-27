#include "../include/types.h"
#include "../include/constants.h"
#include "../include/config.h"
#include "../include/utils.h"
#include "../include/cuts.h"
#include "../include/acceptance.h"
#include "../include/efficiency.h"
#include "../include/efficiency_no_overlap.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <cmath>

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

// ── MC: árbol CoincTree, mismo corte de selección que mc_analysis() ───────
static bool fillCountsMC(TTree* t, McCounts& counts,
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

    Long64_t nsel = 0;
    for (Long64_t k = 0, n = t->GetEntries(); k < n; ++k) {
        t->GetEntry(k);

        if (target_i != 2)                continue;   // TARGET_U
        if (!hasTrue_i || !hasTrue_j)     continue;
        if (Z0 <= 2 || Z1 <= 2)           continue;
        if (Z0 + Z1 < 80 || Z0 + Z1 > 92) continue;
        if (A0 + A1 < 200)                continue;
        if (ppac0!=8 || amp0_c1<0.05 || amp1_c1<0.05 || amp0_c2<0.05 || amp1_c2<0.05) continue;
        ++nsel;

        int d = cosBin(cos_theta_det, nbins_det);
        if (d < 0) continue;

        int e_eff = findBin(energy_bins_eff, neutronE);
        if (e_eff < 0 || e_eff >= nbins_eff) continue;

        // OJO: hay que usar el bin de beam REAL (cos_theta), no un j=0
        // sintético -- si no, computeEfficiency solo usa dOmega_eff[0][i]
        // como normalización, en vez de la aceptancia real de cada evento.
        int j = cosBin(cos_theta, nbins_beam);
        if (j < 0) continue;

        counts.n[e_eff][j][d] += 1.;
    }
    std::cout << "  seleccionados (MC): " << nsel << "\n";
    return nsel > 0;
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
    Float_t  amp0, amp1;
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

        if (x1 > 5 || x0 < -5.5)                       continue;
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

// ── corre ambos métodos de eficiencia dado un McCounts ya lleno ───────────
static bool computeBoth(McCounts& counts, const Vec2D& acceptance, int nbins_eff,
                         std::vector<EfficiencyResult>& eff_overlap,
                         std::vector<EfficiencyResult>& eff_nooverlap)
{
    poissonErrors(counts);
    eff_overlap.resize(nbins_eff);
    eff_nooverlap.resize(nbins_eff);
    for (int e = 0; e < nbins_eff; ++e) {
        eff_overlap[e]   = computeEfficiency(nbins_det - 1, nbins_beam, nbins_det,
                                              counts.n, counts.u_n, acceptance, e);
        eff_nooverlap[e] = computeEfficiencyNoOverlap(nbins_det - 1, nbins_beam, nbins_det,
                                              counts.n, counts.u_n, acceptance, e);
    }
    return true;
}

// ── grafica overlap vs no-overlap y la diferencia relativa, para un dataset
static void plotComparison(const std::vector<EfficiencyResult>& eff_overlap,
                            const std::vector<EfficiencyResult>& eff_nooverlap,
                            int nbins_eff, const std::string& label,
                            const std::string& outdir)
{
    for (int e = 0; e < nbins_eff; ++e) {
        std::vector<double> x, ex, yO, eyO, yN, eyN, yRat, eyRat;
        for (int i = 0; i < nbins_det; ++i) {
            double eO = eff_overlap[e].eps[i],  ueO = eff_overlap[e].u_eps[i];
            double eN = eff_nooverlap[e].eps[i], ueN = eff_nooverlap[e].u_eps[i];
            if (eO <= 0. || eN <= 0.) continue;

            x.push_back((i + 0.5) / nbins_det);
            ex.push_back(0.5 / nbins_det);
            yO.push_back(eO); eyO.push_back(ueO);
            yN.push_back(eN); eyN.push_back(ueN);

            double r = eN / eO;
            yRat.push_back(r);
            eyRat.push_back(r * std::sqrt(std::pow(ueO/eO,2) + std::pow(ueN/eN,2)));
        }
        if (x.empty()) continue;

        TGraphErrors* gO = new TGraphErrors(x.size(), x.data(), yO.data(), ex.data(), eyO.data());
        TGraphErrors* gN = new TGraphErrors(x.size(), x.data(), yN.data(), ex.data(), eyN.data());
        gO->SetName(Form("eff_overlap_%s_ebin%d", label.c_str(), e));
        gN->SetName(Form("eff_nooverlap_%s_ebin%d", label.c_str(), e));
        gO->SetMarkerStyle(20); gO->SetMarkerColor(kBlue+1);  gO->SetLineColor(kBlue+1);
        gN->SetMarkerStyle(24); gN->SetMarkerColor(kRed+1);   gN->SetLineColor(kRed+1);

        // rango Y conjunto -- si no, al dibujar gO primero con "AP" el eje se
        // autoescala solo a gO y gN puede quedar cortado al superponerla
        double ymin = 1e300, ymax = -1e300;
        for (size_t k = 0; k < yO.size(); ++k) {
            ymin = std::min({ymin, yO[k]-eyO[k], yN[k]-eyN[k]});
            ymax = std::max({ymax, yO[k]+eyO[k], yN[k]+eyN[k]});
        }
        double pad = 0.1 * (ymax - ymin);
        gO->SetMinimum(ymin - pad);
        gO->SetMaximum(ymax + pad);

        TCanvas* cv = new TCanvas(Form("cv_eff_%s_%d", label.c_str(), e), "efficiency comparison", 900, 600);
        gO->SetTitle(Form("%s efficiency ebin %d;|cos#theta_{det}|;#varepsilon (rel.)", label.c_str(), e));
        gO->Draw("ALP");
        gN->Draw("LP same");
        TLegend* leg = new TLegend(0.6, 0.75, 0.89, 0.89);
        leg->AddEntry(gO, "con overlap", "lp");
        leg->AddEntry(gN, "sin overlap", "lp");
        leg->Draw();
        cv->SaveAs((outdir + "eff_compare_" + label + "_ebin" + std::to_string(e) + ".pdf").c_str());

        TGraphErrors* gRat = new TGraphErrors(x.size(), x.data(), yRat.data(), ex.data(), eyRat.data());
        gRat->SetName(Form("eff_ratio_%s_ebin%d", label.c_str(), e));
        gRat->SetTitle(Form("%s ratio (sin overlap / con overlap) ebin %d;|cos#theta_{det}|;ratio", label.c_str(), e));
        gRat->SetMarkerStyle(20);

        TCanvas* cv2 = new TCanvas(Form("cv_ratio_%s_%d", label.c_str(), e), "efficiency ratio", 900, 500);
        gRat->Draw("AP");
        cv2->SaveAs((outdir + "eff_ratio_" + label + "_ebin" + std::to_string(e) + ".pdf").c_str());

        std::cout << label << " ebin " << e << ": diferencia media (sin/con - 1) = ";
        double sum = 0.; int n = 0;
        for (double r : yRat) { sum += (r - 1.); ++n; }
        std::cout << (n ? sum / n * 100. : 0.) << " %\n";
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

    const std::string data_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences.root";
    const std::string data_tree = "events_uranium";   // target 2 = uranio

    const std::string acceptance_file =
        "/Users/nico/Desktop/Tese/Analysis/acceptance_coincidence.csv";

    const std::vector<double> energy_bins_eff = {1,10, 100, 600, 1000};
    const int nbins_eff = (int)energy_bins_eff.size() - 1;

    Vec2D acceptance, dOmega_fine;
    if (!loadAcceptanceCSV(acceptance_file, dOmega_fine)) {
        std::cerr << "Failed to load acceptance CSV\n";
        return;
    }
    acceptance = rebin(dOmega_fine);

    // ── MC ──────────────────────────────────────────────────────────────
    std::vector<EfficiencyResult> eff_mc_overlap, eff_mc_nooverlap;
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
                McCounts counts = makeCounts(nbins_eff);
                fillCountsMC(t, counts, energy_bins_eff, nbins_eff);
                computeBoth(counts, acceptance, nbins_eff, eff_mc_overlap, eff_mc_nooverlap);
            }
            fin->Close();
        }
    }
    if (!eff_mc_overlap.empty())
        plotComparison(eff_mc_overlap, eff_mc_nooverlap, nbins_eff, "mc", outdir);

    // ── Datos ───────────────────────────────────────────────────────────
    std::vector<EfficiencyResult> eff_data_overlap, eff_data_nooverlap;
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
                computeBoth(counts, acceptance, nbins_eff, eff_data_overlap, eff_data_nooverlap);
            }
            fin->Close();
        }
    }
    if (!eff_data_overlap.empty())
        plotComparison(eff_data_overlap, eff_data_nooverlap, nbins_eff, "data", outdir);

    // ── guardar todo lo que sí se pudo calcular ────────────────────────────
    TFile* fout = TFile::Open((outdir + "output_efficiency_overlap_comparison.root").c_str(), "RECREATE");
    saveEff(eff_mc_overlap,     nbins_eff, "mc_overlap");
    saveEff(eff_mc_nooverlap,   nbins_eff, "mc_nooverlap");
    saveEff(eff_data_overlap,   nbins_eff, "data_overlap");
    saveEff(eff_data_nooverlap, nbins_eff, "data_nooverlap");
    fout->Close();

    std::cout << "Guardado: " << outdir << "output_efficiency_overlap_comparison.root\n";
}