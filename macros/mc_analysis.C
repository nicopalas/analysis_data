#include "../include/types.h"
#include "../include/constants.h"
#include "../include/config.h"
#include "../include/utils.h"
#include "../include/acceptance.h"
#include "../include/efficiency.h"
#include "../include/anisotropy.h"
#include "../include/plotting.h"

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

// |cos| -> bin index in [0,1]
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

// ---------------------------------------------------------------------------
void mc_analysis()
{
    const std::string infile =
        "/Users/nico/Desktop/Tese/Analysis/montecarlo/data/mc_complete.root";
    const std::string outdir =
        "/Users/nico/Desktop/Tese/Analysis/montecarlo/output/";
    const std::string acceptance_file =
        "/Users/nico/Desktop/Tese/Analysis/cross_section/acceptance_coincidence.csv";

    const std::vector<double> energy_bins_eff = {10,880};
    const int nbins_eff = (int)energy_bins_eff.size() - 1;

    const int nbins_aniso = 9;
    std::vector<double> energy_bins_aniso = {40, 100, 200, 300, 400, 500, 600, 700, 800, 900 };

    // ── acceptance (same setup as data) ────────────────────────────────────
    Vec2D acceptance, dOmega_fine;
    if (!loadAcceptanceCSV(acceptance_file, dOmega_fine)) {
        std::cerr << "Failed to load acceptance CSV\n";
        return;
    }
    acceptance = rebin(dOmega_fine);

    // ── open tree ──────────────────────────────────────────────────────────
    TFile* fin = TFile::Open(infile.c_str());
    if (!fin || fin->IsZombie()) { std::cerr << "cannot open " << infile << "\n"; return; }
    TTree* t = (TTree*)fin->Get("CoincTree");
    if (!t) { std::cerr << "CoincTree not found\n"; fin->Close(); return; }
    std::cout << "Entries: " << t->GetEntries() << "\n";

    Double_t neutronE, cos_theta, cos_theta_det, cos_theta_cm_i;
    Double_t angle_axis_lab, angle_axis_cm;
    Int_t    coinc_type, target_i, Z0, Z1, A0, A1;
    Bool_t   hasCM, hasTrue_i, hasTrue_j;
    Int_t ppac0;
    Double_t amp0_c1, amp1_c1, amp0_c2, amp1_c2;

    t->SetBranchAddress("neutronE",       &neutronE);
    t->SetBranchAddress("cos_theta",      &cos_theta);
    t->SetBranchAddress("cos_theta_det",  &cos_theta_det);
    t->SetBranchAddress("cos_theta_cm_i", &cos_theta_cm_i);
    t->SetBranchAddress("angle_axis_lab", &angle_axis_lab);
    t->SetBranchAddress("angle_axis_cm",  &angle_axis_cm);
    t->SetBranchAddress("coinc_type",     &coinc_type);
    t->SetBranchAddress("target_i",       &target_i);
    t->SetBranchAddress("Z0", &Z0);  t->SetBranchAddress("Z1", &Z1);
    t->SetBranchAddress("A0", &A0);  t->SetBranchAddress("A1", &A1);
    t->SetBranchAddress("hasCM",     &hasCM);
    t->SetBranchAddress("ppac0",     &ppac0);
    t->SetBranchAddress("hasTrue_i", &hasTrue_i);
    t->SetBranchAddress("hasTrue_j", &hasTrue_j);
    t->SetBranchAddress("amp0_c1", &amp0_c1);
    t->SetBranchAddress("amp1_c1", &amp1_c1);
    t->SetBranchAddress("amp0_c2", &amp0_c2);
    t->SetBranchAddress("amp1_c2", &amp1_c2);

    // Efficiency: computed ONCE from the reconstructed angles, and reused for
    // both anisotropies (assumption: the detector efficiency does not depend on
    // which angular definition we histogram against).
    McCounts eff_counts = makeCounts(nbins_eff);

    // Anisotropy: two versions, same detector-frame binning, different beam angle
    McCounts ani_reco = makeCounts(nbins_aniso);   // beam angle = cos_theta   (reco)
    McCounts ani_cm   = makeCounts(nbins_aniso);   // beam angle = cos_theta_cm_i (CM truth)

    TH1D* h_ang_lab = new TH1D("h_ang_lab", "eje reco vs verdadero (lab);#Delta#theta (deg);counts", 100, 0, 10);
    TH1D* h_ang_cm  = new TH1D("h_ang_cm",  "eje reco vs verdadero (CM);#Delta#theta (deg);counts",  100, 0, 10);
    TH1D* h_dcos    = new TH1D("h_dcos", "cos#theta_{reco} - cos#theta*_{CM};#Deltacos;counts", 200, -0.2, 0.2);
    h_ang_lab->SetDirectory(0); h_ang_cm->SetDirectory(0); h_dcos->SetDirectory(0);

    // ── loop ───────────────────────────────────────────────────────────────
    Long64_t nsel = 0;
    for (Long64_t k = 0, n = t->GetEntries(); k < n; ++k) {
        t->GetEntry(k);

        if (target_i != 1)                continue;   // TARGET_U
        if (coinc_type!=1) continue;
        if (ppac0!=7 || amp0_c1<0.05 || amp1_c1<0.05 || amp0_c2<0.05 || amp1_c2<0.05) continue;
        ++nsel;

        if (angle_axis_lab >= 0.) h_ang_lab->Fill(angle_axis_lab);
        if (angle_axis_cm  >= 0.) h_ang_cm ->Fill(angle_axis_cm);
        h_dcos->Fill(std::fabs(cos_theta) - std::fabs(cos_theta_cm_i));

        // detector-frame bin: the SAME for both, always from the reconstructed
        // hit positions — this is what defines the acceptance cell
        int d = cosBin(cos_theta_det, nbins_det);
        if (d < 0) continue;

        int b_reco = cosBin(cos_theta,      nbins_beam);
        int b_cm   = cosBin(cos_theta_cm_i, nbins_beam);

        int e_eff = findBin(energy_bins_eff, neutronE);
        if (e_eff >= 0 && e_eff < nbins_eff && b_reco >= 0)
            eff_counts.n[e_eff][b_reco][d] += 1.;

        int e_ani = findBin(energy_bins_aniso, neutronE);
        if (e_ani >= 0 && e_ani < nbins_aniso) {
            if (b_reco >= 0) ani_reco.n[e_ani][b_reco][d] += 1.;
            if (b_cm   >= 0) ani_cm  .n[e_ani][b_cm  ][d] += 1.;
        }
    }
    fin->Close();

    std::cout << "Selected fission pairs: " << nsel << "\n"
              << "angle_axis_lab: mean=" << h_ang_lab->GetMean()
              << "  RMS=" << h_ang_lab->GetRMS() << " deg\n"
              << "angle_axis_cm : mean=" << h_ang_cm->GetMean()
              << "  RMS=" << h_ang_cm->GetRMS() << " deg\n"
              << "|cos_reco|-|cos_CM|: mean=" << h_dcos->GetMean()
              << "  RMS=" << h_dcos->GetRMS() << "\n";

    poissonErrors(eff_counts);
    poissonErrors(ani_reco);
    poissonErrors(ani_cm);

    // ── efficiency: computed once, shared by both ──────────────────────────
    std::vector<EfficiencyResult> eff(nbins_eff);
    for (int e = 0; e < nbins_eff; ++e)
        eff[e] = computeEfficiency(nbins_det-1, nbins_beam, nbins_det,
                                   eff_counts.n, eff_counts.u_n, acceptance, e);

    // ── anisotropy: two beam-angle definitions, same eff and acceptance ────
    AnalysisConfig cfg = makeUraniumConfig(energy_bins_aniso, "mc");

    std::vector<AnisotropyResult> ani_r(nbins_aniso), ani_c(nbins_aniso);
    for (int e = 0; e < nbins_aniso; ++e) {
        double Ec = std::sqrt(energy_bins_aniso[e] * energy_bins_aniso[e+1]);
        int e_eff = findBin(energy_bins_eff, Ec);
        if (e_eff < 0)          e_eff = 0;
        if (e_eff >= nbins_eff) e_eff = nbins_eff - 1;

        ani_r[e] = anisotropy(nbins_beam, nbins_det, ani_reco.n, ani_reco.u_n,
                              acceptance, e, eff[e_eff].eps, eff[e_eff].u_eps, cfg);
        ani_c[e] = anisotropy(nbins_beam, nbins_det, ani_cm.n, ani_cm.u_n,
                              acceptance, e, eff[e_eff].eps, eff[e_eff].u_eps, cfg);

        std::cout << "ebin " << e << "  E=" << Ec << " MeV"
                  << "   W0/W90  reco=" << ani_r[e].w[nbins_beam-1]
                  << " +/- "            << ani_r[e].u_w[nbins_beam-1]
                  << "   CM="           << ani_c[e].w[nbins_beam-1]
                  << " +/- "            << ani_c[e].u_w[nbins_beam-1] << "\n";
    }

    // ── W(0)/W(90) vs energy ───────────────────────────────────────────────
    std::vector<double> xE, exE, yR, eyR, yC, eyC, yRat, eyRat;
    for (int e = 0; e < nbins_aniso; ++e) {
        if (ani_c[e].w[nbins_beam-1] <= 0. || ani_r[e].w[nbins_beam-1] <= 0.) continue;
        double Ec = std::sqrt(energy_bins_aniso[e] * energy_bins_aniso[e+1]);
        double wr = ani_r[e].w[nbins_beam-1], uwr = ani_r[e].u_w[nbins_beam-1];
        double wc = ani_c[e].w[nbins_beam-1], uwc = ani_c[e].u_w[nbins_beam-1];

        xE .push_back(Ec);
        exE.push_back((energy_bins_aniso[e+1] - energy_bins_aniso[e]) / 2.);
        yR .push_back(wr);  eyR.push_back(uwr);
        yC .push_back(wc);  eyC.push_back(uwc);

        double r = wr / wc;
        yRat .push_back(r);
        eyRat.push_back(r * std::sqrt(std::pow(uwr/wr, 2) + std::pow(uwc/wc, 2)));
    }

    TGraphErrors* gR = new TGraphErrors(xE.size(), xE.data(), yR.data(), exE.data(), eyR.data());
    TGraphErrors* gC = new TGraphErrors(xE.size(), xE.data(), yC.data(), exE.data(), eyC.data());
    gR->SetName("aniso_reco"); gR->SetMarkerStyle(20); gR->SetMarkerColor(kRed+1);   gR->SetLineColor(kRed+1);
    gC->SetName("aniso_cm");   gC->SetMarkerStyle(22); gC->SetMarkerColor(kGreen+2); gC->SetLineColor(kGreen+2);

    TCanvas* cv = new TCanvas("cv", "anisotropy MC", 900, 600);
    cv->SetLogx();
    gC->SetTitle("W(0)/W(90);E_{n} (MeV);W(0)/W(90)");
    gC->Draw("AP");
    gR->Draw("P same");
    TLegend* leg = new TLegend(0.62, 0.75, 0.89, 0.89);
    leg->AddEntry(gC, "cos#theta* (CM)",  "lp");
    leg->AddEntry(gR, "cos#theta (reco)", "lp");
    leg->Draw();
    cv->SaveAs((outdir + "anisotropy_mc_compare.pdf").c_str());

    TGraphErrors* gRat = new TGraphErrors(xE.size(), xE.data(), yRat.data(), exE.data(), eyRat.data());
    gRat->SetName("aniso_ratio_reco_cm");
    gRat->SetTitle("W_{reco}/W_{CM};E_{n} (MeV);ratio");
    gRat->SetMarkerStyle(20);

    TCanvas* cv2 = new TCanvas("cv2", "bias", 900, 500);
    cv2->SetLogx();
    gRat->Draw("AP");
    cv2->SaveAs((outdir + "anisotropy_mc_bias.pdf").c_str());

    // ── save ───────────────────────────────────────────────────────────────
    TFile* fout = TFile::Open((outdir + "output_mc_analysis_au.root").c_str(), "RECREATE");
    for (int e = 0; e < nbins_eff; ++e) {
        TH1D* h = new TH1D(Form("heff_ebin%d", e), "", nbins_det, 0, 1);
        for (int i = 0; i < nbins_det; ++i) {
            h->SetBinContent(i+1, eff[e].eps[i]);
            h->SetBinError  (i+1, eff[e].u_eps[i]);
        }
        h->Write();
    }
    gR->Write(); gC->Write(); gRat->Write();
    h_ang_lab->Write(); h_ang_cm->Write(); h_dcos->Write();
    fout->Close();

    plotEfficiency(eff, nbins_eff, nbins_det, energy_bins_eff,
                   outdir + "efficiency_mc.pdf");
    plotAnisotropy(ani_c, nbins_aniso, nbins_beam, energy_bins_aniso,
                   outdir + "anisotropy_mc_cm.pdf");
    plotAnisotropy(ani_r, nbins_aniso, nbins_beam, energy_bins_aniso,
                   outdir + "anisotropy_mc_reco.pdf");
}