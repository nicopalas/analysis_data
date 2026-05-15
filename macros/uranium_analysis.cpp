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
#include "../include/plotting.h"

void uranium_analysis(){
    const std::vector<double> energy_bins = {1, 10, 100, 500, 1000};
    const int nbins = (int)energy_bins.size() - 1;

    AnalysisConfig cfg = makeUraniumConfig(energy_bins, "nominal");

    // --- open data ---
    TFile* fin  = TFile::Open(cfg.input_file.c_str());
    TTree* tree = (TTree*)fin->Get(cfg.tree_name.c_str());

    // --- histograms ---
    std::vector<TH1D*> hists_tof(nbins, nullptr);
    for(int i = 0; i < nbins; ++i){
        hists_tof[i] = new TH1D(Form("htof_%d", i), "", 200, -30, 30);
        hists_tof[i]->SetDirectory(0);
    }

    Vec3D counts_roi(nbins, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_bkg(nbins, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    // --- fill ---
    fillHistograms(tree, cfg, hists_tof, counts_roi, counts_bkg);
    fin->Close();

    // --- fit background ---
    std::vector<BackgroundFit> bfs(nbins);
    for(int i = 0; i < nbins; ++i){
        double Ec = std::sqrt(energy_bins[i] * energy_bins[i+1]);
        EventCuts c = getCuts(cfg.sample, Ec);
        bfs[i] = fitBackground(cfg, hists_tof[i], c.roi_min, c.roi_max, i);
        std::cout << "Ebin " << i << "  chi2/ndf=" << bfs[i].chi2ndf << "\n";
    }

    // --- signal ---
    Vec3D counts_signal(nbins, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_counts_signal(nbins, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    std::vector<double> counts_subtract(nbins), u_counts_subtract(nbins);
    for(int i = 0; i < nbins; ++i){
        counts_subtract[i]   = bfs[i].counts_subtract;
        u_counts_subtract[i] = bfs[i].u_counts_subtract;
    }

    computeSignal(cfg, counts_roi, counts_bkg,
                  counts_subtract, u_counts_subtract,
                  counts_signal, u_counts_signal);

    // --- acceptance ---
    Vec2D dOmega_fine;
    loadAcceptanceCSV(cfg.acceptance_file, dOmega_fine);
    Vec2D dOmega_eff = rebin(dOmega_fine);

    // --- efficiency ---
    std::vector<EfficiencyResult> eff;
    for(int e = 0; e < nbins; ++e)
        eff.push_back(computeEfficiency(
            nbins_det - 1,   // ref_bin
            nbins_beam,
            nbins_det,
            counts_signal,
            u_counts_signal,
            dOmega_eff,
            e));

    // --- save ---
    TFile* fout = TFile::Open("/Users/nico/Desktop/Tese/Analysis/U-238/output/output_efficiency_uranium.root", "RECREATE");
    for(int e = 0; e < nbins; ++e){
        TH1D* heff = new TH1D(Form("heff_ebin%d", e), "", nbins_det, 0, 1);
        for(int i = 0; i < nbins_det; ++i){
            heff->SetBinContent(i+1, eff[e].eps[i]);
            heff->SetBinError(i+1,   eff[e].u_eps[i]);
        }
        heff->Write();
        hists_tof[e]->Write();
        bfs[e].hist_subtracted->Write();
    }
    fout->Close();

    // --- plot ---
    std::vector<TF1*>  fits(nbins, nullptr);
    std::vector<TH1D*> hists_sub(nbins, nullptr);
    for(int i = 0; i < nbins; ++i){
        fits[i]     = bfs[i].func;
        hists_sub[i] = bfs[i].hist_subtracted;
    }
    plotBackgroundFits(hists_tof, hists_sub, fits, nbins,
                       "/Users/nico/Desktop/Tese/Analysis/U-238/output/background_subtraction_uranium.pdf");
    plotEfficiency(eff, nbins, nbins_det, energy_bins,
                   "/Users/nico/Desktop/Tese/Analysis/U-238/output/efficiency_uranium.pdf");
}