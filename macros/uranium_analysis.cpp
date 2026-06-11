#include "../include/types.h"
#include "../include/constants.h"
#include "../include/config.h"
#include "../include/cuts.h"
#include "../include/utils.h"
#include "../include/acceptance.h"
#include "../include/histograms.h"
#include "../include/background_subtraction.h"
#include "../include/signal_selection.h"
#include "../include/efficiency.h"
#include "../include/anisotropy.h"
#include "../include/plotting.h"
#include "../include/cross_section.h"
#include "TH2D.h"

void uranium_analysis(){

    const double emin = 1.;
    const double emax = 1000;

    const std::string outdir = "/Users/nico/Desktop/Tese/Analysis/cross_section/output/U-238/";

    // ================================================================
    // EFFICIENCY — coarse binning
    // ================================================================
    const std::vector<double> energy_bins_eff = {1, 10, 100, 500, 1000};
    const int nbins_eff = (int)energy_bins_eff.size() - 1;

    AnalysisConfig cfg_eff = makeUraniumConfig(energy_bins_eff, "eff");

    // --- open data ---
    TFile* fin = TFile::Open(cfg_eff.input_file.c_str());
    if(!fin || fin->IsZombie()){
        std::cerr << "Error opening " << cfg_eff.input_file << "\n";
        return;
    }
    TTree* tree = (TTree*)fin->Get(cfg_eff.tree_name.c_str());
    if(!tree){
        std::cerr << "Tree not found\n";
        fin->Close();
        return;
    }
    std::cout << "Entries: " << tree->GetEntries() << "\n";

    // --- histograms ---
    std::vector<TH1D*> hists_tof_eff(nbins_eff, nullptr);
    for(int i = 0; i < nbins_eff; ++i){
        hists_tof_eff[i] = new TH1D(Form("htof_eff_%d", i), "", 200, -30, 30);
        hists_tof_eff[i]->SetDirectory(0);
    }

    Vec3D counts_roi_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_bkg_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_upeak_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    fillHistograms(tree, cfg_eff, hists_tof_eff,
                   counts_roi_eff, counts_bkg_eff, counts_upeak_eff,
                   emin, emax);
    fin->Close();

    for(int i = 0; i < nbins_eff; ++i)
        std::cout << "Eff ebin " << i
                  << "  entries=" << hists_tof_eff[i]->GetEntries() << "\n";

    // --- fit background ---
    std::vector<BackgroundFit> bfs_eff(nbins_eff);
    for(int i = 0; i < nbins_eff; ++i){
        double Ec = std::sqrt(energy_bins_eff[i] * energy_bins_eff[i+1]);
        EventCuts c = getCuts(cfg_eff.sample, Ec);
        bfs_eff[i] = fitBackground(cfg_eff, hists_tof_eff[i], c.roi_min, c.roi_max, i);
        std::cout << "Eff ebin " << i << "  chi2/ndf=" << bfs_eff[i].chi2ndf << "\n";
        if(!bfs_eff[i].func || !bfs_eff[i].hist_subtracted){
            std::cerr << "Eff ebin " << i << " fit failed\n";
            return;
        }
    }

    std::vector<double> cs_eff(nbins_eff),            u_cs_eff(nbins_eff);
    std::vector<double> cs_upeak_eff(nbins_eff, 0.0), u_cs_upeak_eff(nbins_eff, 0.0);
    for(int i = 0; i < nbins_eff; ++i){
        cs_eff[i]   = bfs_eff[i].counts_subtract_bkg;
        u_cs_eff[i] = bfs_eff[i].u_counts_subtract_bkg;
    }

    Vec3D counts_signal_eff(nbins_eff,
    Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_counts_signal_eff(nbins_eff,
    Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    computeSignal(cfg_eff,
                  counts_roi_eff, counts_bkg_eff, counts_upeak_eff,
                  cs_eff,         u_cs_eff,
                  cs_upeak_eff,   u_cs_upeak_eff,
                  counts_signal_eff, u_counts_signal_eff);

    // --- acceptance ---
    TFile *solid_angle = TFile::Open("/Users/nico/Desktop/Tese/Analysis/mc_acceptance.root", "READ");
    TH2D* hist_theta_det = (TH2D*) solid_angle->Get("theta_det_beam");
    TH1D *acceptance_hist = (TH1D*) solid_angle->Get("acceptance");

    std::vector<std::vector<double>> dOmega_eff(nbins_beam, std::vector<double>(nbins_det, 0.0));
    for (int i = 1 ; i <= nbins_beam; i++){  
        for (int j = 1 ; j <= nbins_det; j++){  
            dOmega_eff[i-1][j-1] = hist_theta_det->GetBinContent(i, j);
        }
    }
    solid_angle->Close();
    // --- efficiency ---
    std::vector<EfficiencyResult> eff(nbins_eff);
    for(int e = 0; e < nbins_eff; ++e)
        eff[e] = computeEfficiency(
            nbins_beam,
            nbins_det,
            counts_signal_eff,
            u_counts_signal_eff,
            e);

    // --- save efficiency ---
    TFile* fout_eff = TFile::Open(
        (outdir + "output_efficiency_uranium.root").c_str(), "RECREATE");
    if(!fout_eff || fout_eff->IsZombie()){
        std::cerr << "Error creating efficiency output file\n";
        return;
    }
    for(int e = 0; e < nbins_eff; ++e){
        TH1D* heff = new TH1D(Form("heff_ebin%d", e), "", nbins_det, 0, 1);
        for(int i = 0; i < nbins_det; ++i){
            heff->SetBinContent(i+1, eff[e].eps[i]);
            heff->SetBinError(i+1,   eff[e].u_eps[i]);
        }
        heff->Write();
        hists_tof_eff[e]->Write();
        bfs_eff[e].hist_subtracted->Write();
    }
    fout_eff->Close();

    // --- plot efficiency ---
    std::vector<TF1*>  fits_eff(nbins_eff, nullptr);
    std::vector<TH1D*> hists_sub_eff(nbins_eff, nullptr);
    for(int i = 0; i < nbins_eff; ++i){
        fits_eff[i]      = bfs_eff[i].func;
        hists_sub_eff[i] = bfs_eff[i].hist_subtracted;
    }
    plotBackgroundFits(hists_tof_eff, hists_sub_eff, fits_eff, nbins_eff,
                       outdir + "background_subtraction_uranium_eff.pdf");
    plotEfficiency(eff, nbins_eff, nbins_det, energy_bins_eff,
                   outdir + "efficiency_uranium.pdf");

    // ================================================================
    // ANISOTROPY — fine logarithmic binning
    // ================================================================
    const int nbins_aniso = 45;
    std::vector<double> energy_bins_aniso = buildLogBins(nbins_aniso, 1.0, 1000.0);

    AnalysisConfig cfg_aniso = makeUraniumConfig(energy_bins_aniso, "aniso");

    TFile* fin2 = TFile::Open(cfg_aniso.input_file.c_str());
    if(!fin2 || fin2->IsZombie()){
        std::cerr << "Error opening file for anisotropy\n";
        return;
    }
    TTree* tree2 = (TTree*)fin2->Get(cfg_aniso.tree_name.c_str());

    std::vector<TH1D*> hists_tof_aniso(nbins_aniso, nullptr);
    for(int i = 0; i < nbins_aniso; ++i){
        hists_tof_aniso[i] = new TH1D(Form("htof_aniso_%d", i), "", 200, -20, 20);
        hists_tof_aniso[i]->SetDirectory(0);
    }

    Vec3D counts_roi_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_bkg_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_upeak_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    fillHistograms(tree2, cfg_aniso, hists_tof_aniso,
                   counts_roi_aniso, counts_bkg_aniso, counts_upeak_aniso,
                   emin, emax);
    fin2->Close();

    // --- fit background ---
    std::vector<BackgroundFit> bfs_aniso(nbins_aniso);
    for(int i = 0; i < nbins_aniso; ++i){
        double Ec = std::sqrt(energy_bins_aniso[i] * energy_bins_aniso[i+1]);
        EventCuts c = getCuts(cfg_aniso.sample, Ec);
        bfs_aniso[i] = fitBackground(cfg_aniso, hists_tof_aniso[i],
                                      c.roi_min, c.roi_max, i);
        std::cout << "Aniso ebin " << i << "  chi2/ndf=" << bfs_aniso[i].chi2ndf << "\n";
    }

    // ── diagnóstico señal y estadística ──────────────────────────────────
for(int i = 0; i < nbins_aniso; ++i){
    double Ec = std::sqrt(energy_bins_aniso[i] * energy_bins_aniso[i+1]);
    std::cout << "DIAG ebin " << i
              << "  E="          << Ec                              << " MeV"
              << "  entries="    << hists_tof_aniso[i]->GetEntries()
              << "  signal="     << bfs_aniso[i].counts_subtract_bkg
              << "  u_signal="   << bfs_aniso[i].u_counts_subtract_bkg
              << "  chi2/ndf="   << bfs_aniso[i].chi2ndf            << "\n";
}


    std::vector<double> cs_aniso(nbins_aniso),            u_cs_aniso(nbins_aniso);
    std::vector<double> cs_upeak_aniso(nbins_aniso, 0.0), u_cs_upeak_aniso(nbins_aniso, 0.0);
    for(int i = 0; i < nbins_aniso; ++i){
        cs_aniso[i]   = bfs_aniso[i].counts_subtract_bkg;
        u_cs_aniso[i] = bfs_aniso[i].u_counts_subtract_bkg;
    }
    Vec3D counts_signal_aniso(nbins_aniso,
    Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_counts_signal_aniso(nbins_aniso,
    Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    computeSignal(cfg_aniso,
                  counts_roi_aniso, counts_bkg_aniso, counts_upeak_aniso,
                  cs_aniso,         u_cs_aniso,
                  cs_upeak_aniso,   u_cs_upeak_aniso,
                  counts_signal_aniso, u_counts_signal_aniso);
    std::vector<AnisotropyResult> aniso(nbins_aniso);
    for(int e = 0; e < nbins_aniso; ++e){
        double Ec = std::sqrt(energy_bins_aniso[e] * energy_bins_aniso[e+1]);

        int e_eff = findBin(energy_bins_eff, Ec);
        if(e_eff < 0)          e_eff = 0;
        if(e_eff >= nbins_eff) e_eff = nbins_eff - 1;

        aniso[e] = anisotropy(
            nbins_beam, nbins_det,
            counts_signal_aniso,
            u_counts_signal_aniso,
            dOmega_eff,
            e,
            eff[e_eff].eps,
            eff[e_eff].u_eps,
            cfg_aniso);

        std::cout << "Aniso ebin " << e
                  << "  E=" << Ec << " MeV"
                  << "  W0/W90=" << aniso[e].w[nbins_beam-1]
                  << " +/- "     << aniso[e].u_w[nbins_beam-1] << "\n";
    }

    // --- save anisotropy ---
    TFile* fout_aniso = TFile::Open(
        (outdir + "output_anisotropy_uranium_complete.root").c_str(), "RECREATE");
    if(!fout_aniso || fout_aniso->IsZombie()){
        std::cerr << "Error creating anisotropy output file\n";
        return;
    }
    for(int e = 0; e < nbins_aniso; ++e){
        TGraphErrors* g = nullptr;

        if(e % 10 == 0){
            std::vector<double> x(nbins_beam), y(nbins_beam), ex(nbins_beam, 0.0);
            for(int i = 0; i < nbins_beam; ++i){
                x[i] = (i + 0.5) * dcos_beam;
                y[i] = aniso[e].w[i];
            }
            g = new TGraphErrors(
                nbins_beam,
                x.data(), y.data(),
                ex.data(), aniso[e].u_w.data());
            g->SetName(Form("anisotropy_ebin%d", e));
            g->SetTitle(Form(
                "W(#theta)/W(90) %.1f-%.1f MeV;"
                "cos(#theta_{beam});W(#theta)/W(90)",
                energy_bins_aniso[e], energy_bins_aniso[e+1]));
            g->SetMinimum(0.8);
            g->SetMaximum(2.5);
        }

        if(g) g->Write();
        hists_tof_aniso[e]->Write();
    }

    // --- plot anisotropy ---
    plotAnisotropy(aniso, nbins_aniso, nbins_beam, energy_bins_aniso,
                   outdir + "anisotropy_uranium_complete.pdf");
    plotAnisotropyRatio(aniso, nbins_aniso, nbins_beam, energy_bins_aniso,
                        outdir + "anisotropy_ratio_uranium_complete");
    fout_aniso->Close();


    // ================================================================
// CROSS SECTION ABSOLUTE — fine logarithmic binning
// ================================================================
const int nbins_cs = nbins_aniso;
const std::vector<double>& energy_bins_cs = energy_bins_aniso;

// E_low y E_high en MeV para readFlux
std::vector<double> E_low_cs(nbins_cs), E_high_cs(nbins_cs);
for (int e = 0; e < nbins_cs; ++e) {
    E_low_cs[e]  = energy_bins_cs[e];
    E_high_cs[e] = energy_bins_cs[e+1];
}

// reutilizamos counts_signal_aniso, u_counts_signal_aniso y eff ya calculados
std::vector<CrossSection> cs_abs(nbins_cs);
for (int e = 0; e < nbins_cs; ++e) {

    double Ec   = std::sqrt(energy_bins_cs[e] * energy_bins_cs[e+1]);
    int e_eff   = findBin(energy_bins_eff, Ec);
    if (e_eff < 0)          e_eff = 0;
    if (e_eff >= nbins_eff) e_eff = nbins_eff - 1;

    cs_abs[e] = cross_section_absolute(
        cfg_aniso.flux_file,
        cfg_aniso.flux_hist,
        nbins_beam, nbins_det,
        counts_signal_aniso,      u_counts_signal_aniso,
        eff[e_eff].eps,           eff[e_eff].u_eps,
        dOmega_eff,
        e,
        E_low_cs, E_high_cs,
        cfg_aniso.atoms,
        cfg_aniso);

}

// ── signal counts diagnostic ──────────────────────────────────────────
std::cout << "\n=== Signal counts per energy bin ===\n";
for (int e = 0; e < nbins_cs; ++e) {
    double Ec = std::sqrt(energy_bins_cs[e] * energy_bins_cs[e+1]);
    double total_signal = 0.;
    for (int j = 0; j < nbins_beam; ++j)
        for (int ii = 0; ii < nbins_det; ++ii)
            total_signal += counts_signal_aniso[e][j][ii];
    std::cout << "ebin " << e
              << "  E="          << Ec           << " MeV"
              << "  signal_sum=" << total_signal
              << "  sigma="      << cs_abs[e].sigma
              << " +/- "         << cs_abs[e].u_sigma << " barn\n";
}
for (int e = 0; e < nbins_eff; ++e){
    double Ec = std::sqrt(energy_bins_eff[e] * energy_bins_eff[e+1]);
    std::cout << "eff[" << e << "]  Ec=" << Ec << " MeV\n";
    for (int i = 0; i < nbins_det; ++i){
        std::cout << "  det " << i
                  << "  eps=" << eff[e].eps[i]
                  << " +/- "  << eff[e].u_eps[i] << "\n";
    }
}
// --- save cross section ---
TFile* fout_cs = TFile::Open(
    (outdir + "output_cross_section_uranium.root").c_str(), "RECREATE");
if (!fout_cs || fout_cs->IsZombie()) {
    std::cerr << "Error creating cross section output file\n";
    return;
}

// TGraphErrors con sigma vs energia
std::vector<double> x_cs, y_cs, ex_cs, ey_cs;
for (int e = 0; e < nbins_cs; ++e) {
    if (cs_abs[e].sigma <= 0.0) continue;
    double Ec = std::sqrt(energy_bins_cs[e] * energy_bins_cs[e+1]);
    x_cs.push_back(Ec);
    y_cs.push_back(cs_abs[e].sigma);
    ex_cs.push_back((E_high_cs[e] - E_low_cs[e]) / 2.0);
    ey_cs.push_back(cs_abs[e].u_sigma);
}

TGraphErrors* g_cs = new TGraphErrors(
    (int)x_cs.size(),
    x_cs.data(), y_cs.data(),
    ex_cs.data(), ey_cs.data());
g_cs->SetName("g_sigma_abs_U238");
g_cs->SetTitle("#sigma_{abs}(U-238);E_{n} (MeV);#sigma (barn)");
g_cs->SetMarkerStyle(20);
g_cs->SetMarkerColor(kRed+1);
g_cs->SetLineColor(kRed+1);
g_cs->SetLineWidth(2);
g_cs->Write();

// ── normalised cross section ──────────────────────────────────────────
TH1D* h_cs_raw = new TH1D("h_cs_raw", "", nbins_cs, energy_bins_cs.data());
for (int e = 0; e < nbins_cs; ++e){
    h_cs_raw->SetBinContent(e+1, cs_abs[e].sigma);
    h_cs_raw->SetBinError  (e+1, cs_abs[e].u_sigma);
}
TH1D* h_cs_norm = (TH1D*)h_cs_raw->Clone("h_cs_norm");
h_cs_norm->Reset();
normalised_xs(h_cs_raw, 8.0, 10.0, 2.006, h_cs_norm);

TCanvas* c_cs = new TCanvas("c_cs", "Cross section U-238", 900, 600);
c_cs->SetLogx();
h_cs_norm->SetLineColor(kBlue+1);
h_cs_norm->SetMarkerColor(kBlue+1);
h_cs_norm->SetMarkerStyle(20);
h_cs_norm->GetXaxis()->SetTitle("E_{n} (MeV)");
h_cs_norm->GetYaxis()->SetTitle("#sigma (barn)");
h_cs_norm->SetTitle("#sigma_{norm}(U-238)");
h_cs_norm->Draw("E");

fout_cs->cd();
g_cs->Write();
h_cs_norm->Write();
c_cs->Write();
fout_cs->Close();

for(int i = 0; i < nbins_eff;   ++i) delete hists_tof_eff[i];
for(int i = 0; i < nbins_aniso; ++i) delete hists_tof_aniso[i];
}