#include "cuts.h"
#include "histograms.h"
#include "background.h"
#include "signal.h"
#include "efficiency.h"
#include "plotting.h"
#include "acceptance.h"   // loadAcceptanceCSV + rebin

void uranium_analysis(){
    const int nbins = 4;
    const std::vector<double> energy_bins = {1, 10, 100, 500, 1000};

    // --- open data ---
    TFile* fin  = TFile::Open("coincidences.root");
    TTree* tree = (TTree*)fin->Get("events_uranium");

    // --- initialize histograms ---
    TH1D* hists_tof[4];
    for(int i = 0; i < nbins; ++i){
        hists_tof[i] = new TH1D(Form("htof_%d", i), "", 200, -30, 30);
        hists_tof[i]->SetDirectory(0);
    }

    double counts_roi[4][10][10] = {};
    double counts_bkg[4][10][10] = {};

    // --- fill ---
    fillHistograms(tree, energy_bins, hists_tof, counts_roi, counts_bkg, nbins);
    fin->Close();

    // --- fit background ---
    double counts_subtract[4]   = {};
    double u_counts_subtract[4] = {};
    TF1*   fits[4]              = {};
    TH1D*  hists_sub[4]         = {};

    for(int i = 0; i < nbins; ++i){
        double Ec = std::sqrt(energy_bins[i] * energy_bins[i+1]);
        EventCuts c = getCuts(Ec);
        BackgroundFit bf = fitBackground(
            hists_tof[i], c.roi_min, c.roi_max, 500, i, "nominal");

        counts_subtract[i]   = bf.counts_subtract;
        u_counts_subtract[i] = bf.u_counts_subtract;
        fits[i]              = bf.func;
        hists_sub[i]         = bf.hist_subtracted;

        std::cout << "Ebin " << i << "  chi2/ndf=" << bf.chi2ndf << "\n";
    }

    // --- signal ---
    double counts_signal[4][10][10]   = {};
    double u_counts_signal[4][10][10] = {};
    computeSignal(nbins, 10, 10,
                  counts_roi, counts_bkg,
                  counts_subtract, u_counts_subtract,
                  counts_signal, u_counts_signal);

    // --- acceptance ---
    double dOmega_fine[100][20] = {};
    double dOmega_eff[10][10]   = {};
    loadAcceptanceCSV("acceptance_coincidence.csv", dOmega_fine);
    rebin(dOmega_fine, dOmega_eff);

    // --- efficiency ---
    EfficiencyResult eff[4];
    for(int e = 0; e < nbins; ++e)
        eff[e] = computeEfficiency(10, 10, 9,
                                   counts_signal, u_counts_signal,
                                   dOmega_eff, e);

    // --- save ---
    TFile* fout = TFile::Open("output_efficiency_uranium.root", "RECREATE");
    for(int e = 0; e < nbins; ++e){
        TH1D* heff = new TH1D(Form("heff_ebin%d", e), "", 10, 0, 1);
        for(int i = 0; i < 10; ++i){
            heff->SetBinContent(i+1, eff[e].eps[i]);
            heff->SetBinError(i+1,   eff[e].u_eps[i]);
        }
        heff->Write();
        hists_tof[e]->Write();
        hists_sub[e]->Write();
    }
    fout->Close();

    // --- plot ---
    plotBackgroundFits(hists_tof, hists_sub, fits, nbins,
                       "background_subtraction_uranium.pdf");
    plotEfficiency(eff, nbins, 10, energy_bins,
                   "efficiency_uranium.pdf");

    return 0;
}