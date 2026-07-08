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



void gold_analysis()
{
    const std::string outdir = "/Users/nico/Desktop/Tese/Analysis/cross_section/output/Au-197/";

    // ================================================================
    // EFFICIENCY — coarse binning
    // ================================================================
    const double emin = 100.;
    const double emax = 1000.;
    const std::vector<double> energy_bins_eff = {80, 300, 600, 1000};
    const int nbins_eff = (int)energy_bins_eff.size() - 1;

    AnalysisConfig cfg_eff = makeGoldConfig(energy_bins_eff, "eff");

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
    std::cout << "Gold entries: " << tree->GetEntries() << "\n";

    // --- histograms ---
    std::vector<TH1D*> hists_tof_eff(nbins_eff, nullptr);
    for(int i = 0; i < nbins_eff; ++i){
        hists_tof_eff[i] = new TH1D(
            Form("htof_gold_eff_%d", i), "", 100, -20, 20);
        hists_tof_eff[i]->SetDirectory(0);
    }

    Vec3D counts_roi_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_bkg_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_upeak_eff(nbins_eff,                                  // NEW
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    fillHistograms(tree, cfg_eff, hists_tof_eff,
                   counts_roi_eff, counts_bkg_eff, counts_upeak_eff,  // NEW
                   emin, emax);
    fin->Close();

    for(int i = 0; i < nbins_eff; ++i)
        std::cout << "Gold eff ebin " << i
                  << "  entries=" << hists_tof_eff[i]->GetEntries() << "\n";

    // --- fit background ---
    std::vector<BackgroundFit> bfs_eff(nbins_eff);
    for(int i = 0; i < nbins_eff; ++i){
        double Ec = std::sqrt(energy_bins_eff[i] * energy_bins_eff[i+1]);
        EventCuts c = getCuts(cfg_eff.sample, Ec);
        bfs_eff[i] = fitBackground(cfg_eff, hists_tof_eff[i],
                                   c.roi_min, c.roi_max, i);
        std::cout << "Gold eff ebin " << i
                  << "  chi2/ndf=" << bfs_eff[i].chi2ndf << "\n";
        if(!bfs_eff[i].func || !bfs_eff[i].hist_subtracted){
            std::cerr << "Gold eff ebin " << i << " fit failed\n";
            return;
        }
    }

    // declara aquí, antes de computeSignal
    Vec3D counts_signal_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_counts_signal_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    std::vector<double> cs_bkg_eff(nbins_eff),    u_cs_bkg_eff(nbins_eff);
    std::vector<double> cs_upeak_eff(nbins_eff),   u_cs_upeak_eff(nbins_eff);
    for(int i = 0; i < nbins_eff; ++i){
        cs_bkg_eff[i]     = bfs_eff[i].counts_subtract_bkg;
        u_cs_bkg_eff[i]   = bfs_eff[i].u_counts_subtract_bkg;
        cs_upeak_eff[i]   = bfs_eff[i].counts_subtract_upeak;
        u_cs_upeak_eff[i] = bfs_eff[i].u_counts_subtract_upeak;
    }

    computeSignal(cfg_eff,
                  counts_roi_eff, counts_bkg_eff, counts_upeak_eff,
                  cs_bkg_eff,   u_cs_bkg_eff,
                  cs_upeak_eff, u_cs_upeak_eff,
                  counts_signal_eff, u_counts_signal_eff);

    // --- acceptance ---
    TFile *solid_angle = TFile::Open("/Users/nico/Desktop/Tese/Analysis/mc_acceptance.root", "READ");
    TH2D* hist_theta_det = (TH2D*) solid_angle->Get("theta_det_beam");
    std::vector<std::vector<double>> dOmega_eff(nbins_beam, std::vector<double>(nbins_det, 0.0));

    for (int i = 1; i <= nbins_beam_fine; i++) {
        for (int j = 1; j <= nbins_det_fine; j++) {
            int i_rebin = (i - 1) / 2;
            int j_rebin = (j - 1) / 2;
        
            if (i_rebin < nbins_beam && j_rebin < nbins_det) {
                dOmega_eff[i_rebin][j_rebin] += hist_theta_det->GetBinContent(i, j);
        }
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
        (outdir + "output_efficiency_gold.root").c_str(), "RECREATE");
    if(!fout_eff || fout_eff->IsZombie()){
        std::cerr << "Error creating gold efficiency file\n";
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
                       outdir + "background_subtraction_gold_eff.pdf");
    plotEfficiency(eff, nbins_eff, nbins_det, energy_bins_eff,
                   outdir + "efficiency_gold.pdf");

    // ================================================================
    // ANISOTROPY — 6 bins log entre 150 y 300 MeV
    // ================================================================
    const int nbins_aniso = 7;
    std::vector<double> energy_bins_aniso = buildLogBins(nbins_aniso, 100.0, 1000.0);

    AnalysisConfig cfg_aniso = makeGoldConfig(energy_bins_aniso, "aniso");

    // --- open data ---
    TFile* fin2 = TFile::Open(cfg_aniso.input_file.c_str());
    if(!fin2 || fin2->IsZombie()){
        std::cerr << "Error opening file for gold anisotropy\n";
        return;
    }
    TTree* tree2 = (TTree*)fin2->Get(cfg_aniso.tree_name.c_str());
    if(!tree2){ fin2->Close(); return; }

    // --- histograms ---
    std::vector<TH1D*> hists_tof_aniso(nbins_aniso, nullptr);
    for(int i = 0; i < nbins_aniso; ++i){
        hists_tof_aniso[i] = new TH1D(
            Form("htof_gold_aniso_%d", i), "", 100, -30, 30);
        hists_tof_aniso[i]->SetDirectory(0);
    }

    Vec3D counts_roi_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_bkg_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_upeak_aniso(nbins_aniso,                              // NEW
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    fillHistograms(tree2, cfg_aniso, hists_tof_aniso,
                   counts_roi_aniso, counts_bkg_aniso, counts_upeak_aniso,  // NEW
                   emin, emax);
    fin2->Close();

    // --- fit background ---
    std::vector<BackgroundFit> bfs_aniso(nbins_aniso);
    for(int i = 0; i < nbins_aniso; ++i){
        double Ec = std::sqrt(energy_bins_aniso[i] * energy_bins_aniso[i+1]);
        EventCuts c = getCuts(cfg_aniso.sample, Ec);
        bfs_aniso[i] = fitBackground(cfg_aniso, hists_tof_aniso[i],
                                     c.roi_min, c.roi_max, i);
        std::cout << "Gold aniso ebin " << i
                  << "  chi2/ndf=" << bfs_aniso[i].chi2ndf << "\n";
        if(!bfs_aniso[i].func || !bfs_aniso[i].hist_subtracted){
            std::cerr << "Gold aniso ebin " << i << " fit failed\n";
            return;
        }
    }

    // --- signal ---
    Vec3D counts_signal_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_counts_signal_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    // --- signal aniso ---
    std::vector<double> cs_bkg_aniso(nbins_aniso),   u_cs_bkg_aniso(nbins_aniso);
    std::vector<double> cs_upeak_aniso(nbins_aniso),  u_cs_upeak_aniso(nbins_aniso);
    for(int i = 0; i < nbins_aniso; ++i){
        cs_bkg_aniso[i]    = bfs_aniso[i].counts_subtract_bkg;
        u_cs_bkg_aniso[i]  = bfs_aniso[i].u_counts_subtract_bkg;
        cs_upeak_aniso[i]  = bfs_aniso[i].counts_subtract_upeak;
        u_cs_upeak_aniso[i]= bfs_aniso[i].u_counts_subtract_upeak;
    }
    computeSignal(cfg_aniso, counts_roi_aniso, counts_bkg_aniso, counts_upeak_aniso,
                  cs_bkg_aniso,   u_cs_bkg_aniso,
                  cs_upeak_aniso, u_cs_upeak_aniso,
                  counts_signal_aniso, u_counts_signal_aniso);


    // --- anisotropy ---
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

        std::cout << "Gold aniso ebin " << e
                  << "  E=" << Ec << " MeV"
                  << "  W0/W90=" << aniso[e].w[nbins_beam-1]
                  << " +/- "     << aniso[e].u_w[nbins_beam-1] << "\n";
    }

    // --- save anisotropy ---
    TFile* fout_aniso = TFile::Open(
        (outdir + "output_anisotropy_gold.root").c_str(), "RECREATE");
    if(!fout_aniso || fout_aniso->IsZombie()){
        std::cerr << "Error creating gold anisotropy file\n";
        return;
    }
    for(int e = 0; e < nbins_aniso; ++e){
        std::vector<double> x(nbins_beam), y(nbins_beam), ex(nbins_beam, 0.0);
        for(int i = 0; i < nbins_beam; ++i){
            x[i] = (i + 0.5) * dcos_beam;
            y[i] = aniso[e].w[i];
        }
        TGraphErrors* g = new TGraphErrors(
            nbins_beam,
            x.data(), y.data(),
            ex.data(), aniso[e].u_w.data());
        g->SetName(Form("anisotropy_gold_ebin%d", e));
        g->SetTitle(Form("W(#theta)/W(90) %.1f-%.1f MeV;"
                         "cos(#theta_{beam});W(#theta)/W(90)",
                         energy_bins_aniso[e], energy_bins_aniso[e+1]));
        g->Write();
        hists_tof_aniso[e]->Write();
    }
    fout_aniso->Close();

    // --- plot anisotropy ---
    plotAnisotropy(aniso, nbins_aniso, nbins_beam, energy_bins_aniso,
                   outdir + "anisotropy_gold.pdf");
    plotAnisotropyRatio(aniso, nbins_aniso, nbins_beam, energy_bins_aniso,
                        outdir + "anisotropy_ratio_gold.pdf");

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
        nbins_beam,
        nbins_det,
        counts_signal_aniso,
        u_counts_signal_aniso,
        eff[e_eff].eps,
        eff[e_eff].u_eps,
        dOmega_eff,
        e,
        E_low_cs,
        E_high_cs,
        cfg_aniso.atoms,                
        cfg_aniso);

    std::cout << "CS_abs ebin " << e
              << "  E=" << Ec << " MeV"
              << "  sigma=" << cs_abs[e].sigma
              << " +/- "    << cs_abs[e].u_sigma << " barn\n";
}

// --- save cross section ---
TFile* fout_cs = TFile::Open(
    (outdir + "output_cross_section_gold.root").c_str(), "RECREATE");
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
g_cs->SetName("g_sigma_abs_Au197");
g_cs->SetTitle("#sigma_{abs}(Au-197);E_{n} (MeV);#sigma (barn)");
g_cs->SetMarkerStyle(20);
g_cs->SetMarkerColor(kRed+1);
g_cs->SetLineColor(kRed+1);
g_cs->SetLineWidth(2);
g_cs->Write();

fout_cs->Close();


    for(int i = 0; i < nbins_eff;   ++i) delete hists_tof_eff[i];
    for(int i = 0; i < nbins_aniso; ++i) delete hists_tof_aniso[i];
}