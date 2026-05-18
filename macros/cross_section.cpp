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
#include "../include/cross_section.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"

void cross_section() {

    const std::string outdir           = "/Users/nico/Desktop/Tese/Analysis/U-238/output/";
    const std::string eff_file_uranium = "/Users/nico/Desktop/Tese/Analysis/U-238/output/output_efficiency_uranium.root";
    const std::string eff_file_gold    = "/Users/nico/Desktop/Tese/Analysis/U-238/output/output_efficiency_gold.root";

    // ================================================================
    // CONFIGS — each carries its own energy_bins_eff
    // ================================================================
    const double emin = 150.;
    const double emax = 1000.;
    const std::vector<double> energy_bins = buildLogBins(8, emin, emax);
    const int nbins = (int)energy_bins.size() - 1;

    AnalysisConfig cfg_u = makeUraniumConfig(energy_bins, "xs");
    AnalysisConfig cfg_g = makeGoldConfig(energy_bins, "xs");

    const int nbins_eff_u = (int)cfg_u.energy_bins_eff.size() - 1;
    const int nbins_eff_g = (int)cfg_g.energy_bins_eff.size() - 1;

    // ================================================================
    // LOAD EFFICIENCIES
    // ================================================================
    auto load_efficiency = [&](const std::string& path, int nbins_eff)
        -> std::vector<EfficiencyResult>
    {
        std::vector<EfficiencyResult> eff(nbins_eff);
        TFile* f = TFile::Open(path.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening efficiency file: " << path << "\n";
            return eff;
        }
        for (int e = 0; e < nbins_eff; ++e) {
            TH1D* h = (TH1D*)f->Get(Form("heff_ebin%d", e));
            if (!h) {
                std::cerr << "Missing heff_ebin" << e << " in " << path << "\n";
                continue;
            }
            eff[e].eps.resize(nbins_det, 0.0);
            eff[e].u_eps.resize(nbins_det, 0.0);
            for (int ii = 0; ii < nbins_det; ++ii) {
                eff[e].eps[ii]   = h->GetBinContent(ii + 1);
                eff[e].u_eps[ii] = h->GetBinError(ii + 1);
            }
        }
        f->Close();
        return eff;
    };

    std::vector<EfficiencyResult> eff_uranium = load_efficiency(eff_file_uranium, nbins_eff_u);
    std::vector<EfficiencyResult> eff_gold    = load_efficiency(eff_file_gold,    nbins_eff_g);

    // ================================================================
    // ACCEPTANCE
    // ================================================================
    Vec2D dOmega_fine;
    loadAcceptanceCSV(cfg_u.acceptance_file, dOmega_fine);
    Vec2D dOmega = rebin(dOmega_fine);

    // ================================================================
    // PIPELINE — fill histograms, fit background, compute signal
    // ================================================================
    auto run_pipeline = [&](
        AnalysisConfig& cfg,
        Vec3D& counts_signal,
        Vec3D& u_counts_signal)
    {
        const int nb = (int)cfg.energy_bins.size() - 1;

        TFile* fin = TFile::Open(cfg.input_file.c_str());
        if (!fin || fin->IsZombie()) {
            std::cerr << "Error opening " << cfg.input_file << "\n";
            return;
        }
        TTree* tree = (TTree*)fin->Get(cfg.tree_name.c_str());
        if (!tree) {
            std::cerr << "Tree not found\n";
            fin->Close();
            return;
        }

        std::string tag = (cfg.sample == Sample::uranium) ? "uranium" : "gold";

        std::vector<TH1D*> hists_tof(nb, nullptr);
        for (int i = 0; i < nb; ++i) {
            hists_tof[i] = new TH1D(Form("htof_%s_%d", tag.c_str(), i),
                                     "", 100, -30, 30);
            hists_tof[i]->SetDirectory(0);
        }

        Vec3D counts_roi(nb, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
        Vec3D counts_bkg(nb, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

        fillHistograms(tree, cfg, hists_tof, counts_roi, counts_bkg, emin, emax );
        fin->Close();

        std::vector<double> cs(nb, 0.0), u_cs(nb, 0.0);
        for (int i = 0; i < nb; ++i) {
            double Ec = std::sqrt(cfg.energy_bins[i] * cfg.energy_bins[i+1]);
            EventCuts c = getCuts(cfg.sample, Ec);
            BackgroundFit bf = fitBackground(cfg, hists_tof[i], c.roi_min, c.roi_max, i);
            if (!bf.func || !bf.hist_subtracted) {
                std::cerr << "Fit failed ebin " << i << "\n";
                continue;
            }
            cs[i]   = bf.counts_subtract;
            u_cs[i] = bf.u_counts_subtract;
        }

        counts_signal.assign(nb, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
        u_counts_signal.assign(nb, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
        computeSignal(cfg, counts_roi, counts_bkg, cs, u_cs,
                      counts_signal, u_counts_signal);
    };

    Vec3D counts_signal_uranium, u_counts_signal_uranium;
    Vec3D counts_signal_gold,    u_counts_signal_gold;

    run_pipeline(cfg_u, counts_signal_uranium, u_counts_signal_uranium);
    run_pipeline(cfg_g, counts_signal_gold,    u_counts_signal_gold);

    // ================================================================
    // CROSS SECTION RATIO
    // ================================================================
    std::vector<double> energy(nbins), ratio(nbins), u_ratio(nbins);
    std::vector<double> E_low(nbins), E_high(nbins);

    for (int e = 0; e < nbins; ++e) {
        energy[e] = std::sqrt(energy_bins[e] * energy_bins[e+1]);
        E_low[e]  = energy_bins[e];
        E_high[e] = energy_bins[e+1];

        double Ec = energy[e];

        int e_eff_u = findBin(cfg_u.energy_bins_eff, Ec);
        if (e_eff_u < 0)             e_eff_u = 0;
        if (e_eff_u >= nbins_eff_u)  e_eff_u = nbins_eff_u - 1;

        int e_eff_g = findBin(cfg_g.energy_bins_eff, Ec);
        if (e_eff_g < 0)             e_eff_g = 0;
        if (e_eff_g >= nbins_eff_g)  e_eff_g = nbins_eff_g - 1;

        CrossSection xs = cross_section_ratio(
            nbins_beam, nbins_det,
            counts_signal_gold,    u_counts_signal_gold,
            counts_signal_uranium, u_counts_signal_uranium,
            eff_gold[e_eff_g].eps,    eff_gold[e_eff_g].u_eps,
            eff_uranium[e_eff_u].eps, eff_uranium[e_eff_u].u_eps,
            dOmega, e, cfg_u);

        ratio[e]   = xs.sigma;
        u_ratio[e] = xs.u_sigma;

        std::cout << "Ebin " << e
                  << "  E=[" << E_low[e] << ", " << E_high[e] << "] MeV"
                  << "  e_eff_u=" << e_eff_u << "  e_eff_g=" << e_eff_g
                  << "  Au/U=" << ratio[e] << " +/- " << u_ratio[e] << "\n";
    }

    // ================================================================
    // SAVE AND PLOT
    // ================================================================
    std::vector<double> x, y, ex, ey;
    for (int e = 0; e < nbins; ++e) {
        if (ratio[e] <= 0.0) continue;
        x.push_back(energy[e]);
        y.push_back(ratio[e]);
        ex.push_back(0.0);
        ey.push_back(u_ratio[e]);
    }

    TGraphErrors* g_ratio = new TGraphErrors(
        (int)x.size(), x.data(), y.data(), ex.data(), ey.data());
    g_ratio->SetName("g_ratio_Au_over_U");
    g_ratio->SetTitle("#sigma(Au)/#sigma(U);E (MeV);Au/U");
    g_ratio->SetMarkerStyle(20);
    g_ratio->SetLineWidth(2);

    TCanvas* c1 = new TCanvas("c1", "Au/U cross section ratio", 800, 600);
    gStyle->SetOptStat(0);
    c1->SetLogx();
    g_ratio->Draw("AP");

    TFile* fout = TFile::Open((outdir + "cross_section_ratio.root").c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error creating output file\n";
        return;
    }
    g_ratio->Write();
    c1->Write();
    fout->Close();
}