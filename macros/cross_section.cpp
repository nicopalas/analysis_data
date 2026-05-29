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

    const std::string outdir           = "/Users/nico/Desktop/Tese/Analysis/cross_section/output/";
    const std::string eff_file_uranium = "/Users/nico/Desktop/Tese/Analysis/cross_section/output/U-238/output_efficiency_uranium.root";
    const std::string eff_file_gold    = "/Users/nico/Desktop/Tese/Analysis/cross_section/output/Au-197/output_efficiency_gold.root";

    // ================================================================
    // CONFIGS
    // ================================================================
    const double emin = 100.;
    const double emax = 1000.;
    const std::vector<double> energy_bins = buildLogBins(6, emin, emax);
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
    // PIPELINE
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

        Vec3D counts_roi(nb,
            Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
        Vec3D counts_bkg(nb,
            Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
        Vec3D counts_upeak(nb,
            Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

        fillHistograms(tree, cfg, hists_tof,
                       counts_roi, counts_bkg, counts_upeak,
                       emin, emax);
        fin->Close();

        std::vector<double> cs_bkg(nb, 0.0),   u_cs_bkg(nb, 0.0);
        std::vector<double> cs_upeak(nb, 0.0), u_cs_upeak(nb, 0.0);

        for (int i = 0; i < nb; ++i) {
            double Ec = std::sqrt(cfg.energy_bins[i] * cfg.energy_bins[i+1]);
            EventCuts c = getCuts(cfg.sample, Ec);
            BackgroundFit bf = fitBackground(cfg, hists_tof[i], c.roi_min, c.roi_max, i);
            if (!bf.func || !bf.hist_subtracted) {
                std::cerr << "Fit failed ebin " << i << "\n";
                continue;
            }
            cs_bkg[i]    = bf.counts_subtract_bkg;
            u_cs_bkg[i]  = bf.u_counts_subtract_bkg;
            cs_upeak[i]  = bf.counts_subtract_upeak;
            u_cs_upeak[i]= bf.u_counts_subtract_upeak;
        }

        counts_signal.assign(nb,
            Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
        u_counts_signal.assign(nb,
            Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

        computeSignal(cfg,
                      counts_roi, counts_bkg, counts_upeak,
                      cs_bkg,    u_cs_bkg,
                      cs_upeak,  u_cs_upeak,
                      counts_signal, u_counts_signal);

        for (int i = 0; i < nb; ++i) delete hists_tof[i];
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
        if (e_eff_u < 0)            e_eff_u = 0;
        if (e_eff_u >= nbins_eff_u) e_eff_u = nbins_eff_u - 1;

        int e_eff_g = findBin(cfg_g.energy_bins_eff, Ec);
        if (e_eff_g < 0)            e_eff_g = 0;
        if (e_eff_g >= nbins_eff_g) e_eff_g = nbins_eff_g - 1;

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
    // CROSS SECTION ABSOLUTE
    // ================================================================
    std::vector<double> sigma_u(nbins, 0.0),    u_sigma_u(nbins, 0.0);
    std::vector<double> sigma_au(nbins, 0.0),   u_sigma_au(nbins, 0.0);

    for (int e = 0; e < nbins; ++e) {
        double Ec = energy[e];

        int e_eff_u = findBin(cfg_u.energy_bins_eff, Ec);
        if (e_eff_u < 0)            e_eff_u = 0;
        if (e_eff_u >= nbins_eff_u) e_eff_u = nbins_eff_u - 1;

        int e_eff_g = findBin(cfg_g.energy_bins_eff, Ec);
        if (e_eff_g < 0)            e_eff_g = 0;
        if (e_eff_g >= nbins_eff_g) e_eff_g = nbins_eff_g - 1;

        CrossSection xs_u = cross_section_absolute(
            cfg_u.flux_file, cfg_u.flux_hist,
            nbins_beam, nbins_det,
            counts_signal_uranium, u_counts_signal_uranium,
            eff_uranium[e_eff_u].eps, eff_uranium[e_eff_u].u_eps,
            dOmega, e, E_low, E_high,
            cfg_u.atoms, cfg_u);

        CrossSection xs_g = cross_section_absolute(
            cfg_g.flux_file, cfg_g.flux_hist,
            nbins_beam, nbins_det,
            counts_signal_gold, u_counts_signal_gold,
            eff_gold[e_eff_g].eps, eff_gold[e_eff_g].u_eps,
            dOmega, e, E_low, E_high,
            cfg_g.atoms, cfg_g);

        sigma_u[e]   = xs_u.sigma;
        u_sigma_u[e] = xs_u.u_sigma;
        sigma_au[e]  = xs_g.sigma;
        u_sigma_au[e]= xs_g.u_sigma;

        std::cout << "Abs XS ebin " << e
                  << "  E=" << Ec << " MeV"
                  << "  sigma_U="  << sigma_u[e]  << " +/- " << u_sigma_u[e]
                  << "  sigma_Au=" << sigma_au[e] << " +/- " << u_sigma_au[e]
                  << " barn\n";
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

    TFile* fout = TFile::Open(
        (outdir + "cross_section_ratio.root").c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error creating output file\n";
        return;
    }
    g_ratio->Write();
    c1->Write();

    // graphs absolutos
    auto makeAbsGraph = [&](const std::vector<double>& sig,
                             const std::vector<double>& u_sig,
                             const char* name, const char* title,
                             int color) -> TGraphErrors*
    {
        std::vector<double> gx, gy, gex, gey;
        for (int e = 0; e < nbins; ++e) {
            if (sig[e] <= 0.0) continue;
            gx.push_back(energy[e]);
            gy.push_back(sig[e]);
            gex.push_back(0.0);
            gey.push_back(u_sig[e]);
        }
        TGraphErrors* g = new TGraphErrors(
            (int)gx.size(), gx.data(), gy.data(), gex.data(), gey.data());
        g->SetName(name);
        g->SetTitle(title);
        g->SetMarkerStyle(20);
        g->SetMarkerColor(color);
        g->SetLineColor(color);
        g->SetLineWidth(2);
        return g;
    };

    TGraphErrors* g_abs_u  = makeAbsGraph(sigma_u,  u_sigma_u,
        "g_sigma_abs_U238",  "#sigma_{abs}(U-238);E (MeV);#sigma (barn)",  kRed+1);
    TGraphErrors* g_abs_au = makeAbsGraph(sigma_au, u_sigma_au,
        "g_sigma_abs_Au197", "#sigma_{abs}(Au-197);E (MeV);#sigma (barn)", kBlue+1);

    TCanvas* c_abs = new TCanvas("c_abs", "Absolute cross sections", 800, 600);
    c_abs->SetLogx();
    c_abs->SetLogy();
    g_abs_au->Draw("AP");
    g_abs_u->Draw("P SAME");
    TLegend* leg_abs = new TLegend(0.6, 0.7, 0.88, 0.85);
    leg_abs->AddEntry(g_abs_au, "Au-197", "lp");
    leg_abs->AddEntry(g_abs_u,  "U-238",  "lp");
    leg_abs->Draw();
    g_ratio->Write();
    g_abs_u->Write();
    g_abs_au->Write();
    c1->Write();
    c_abs->Write();
    fout->Close();
}