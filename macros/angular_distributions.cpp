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

// Global variables for energy range
static double emin_global = 1.0;
static double emax_global = 1000.0;

void fill_histograms(TTree* tree, const AnalysisConfig& cfg, 
                     std::vector<TH1D*>& hists_bkg,
                     bool is_accidental = false,
                     bool debug = false) {
    
    if (!tree) {
        std::cerr << "ERROR: tree is null!\n";
        return;
    }
    
    std::cout << "Processing: " << (is_accidental ? "Accidental" : "Data") << "\n";
    std::cout << "Entries: " << tree->GetEntries() << "\n";
    std::cout << "Sample: " << (cfg.sample == Sample::uranium ? "Uranium" : "Gold") << "\n";

    int nbins = (int)cfg.energy_bins.size() - 1;
    
    // ONLY create histograms if the vector is empty
    if (hists_bkg.empty()) {
        for (int i = 0; i < nbins; i++) {
            TH1D* h = new TH1D(Form("bkg_%s_ebin_%d", 
                                    (is_accidental ? "acc" : "data"), i), 
                               Form("Background Angular Distribution E_bin_%d", i),
                               30, 0, 1.0);
            h->SetDirectory(0);
            hists_bkg.push_back(h);
        }
    }

    // Set branches
    double cos_theta_det, cos_theta;
    double tof1, tof0, neutron_energy;
    double amp0, amp1;
    double x0, x1, y0, y1;
    double phi, phi_det;
    int RunNumber, PSpulse;
    
    tree->SetBranchAddress("tof1", &tof1);
    tree->SetBranchAddress("tof0", &tof0);
    tree->SetBranchAddress("amp0", &amp0);
    tree->SetBranchAddress("amp1", &amp1);
    tree->SetBranchAddress("sumX1", &x1);
    tree->SetBranchAddress("sumX0", &x0);
    tree->SetBranchAddress("y0", &y0);
    tree->SetBranchAddress("sumY1", &y1);
    tree->SetBranchAddress("neutron_energy", &neutron_energy);
    tree->SetBranchAddress("cos_theta", &cos_theta);
    tree->SetBranchAddress("cos_theta_det", &cos_theta_det);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("phi_det", &phi_det);
    tree->SetBranchAddress("RunNumber", &RunNumber);
    tree->SetBranchAddress("PSpulse", &PSpulse);

    Long64_t nentries = tree->GetEntries();
    
    // Debug counters
    long long total_events = 0;
    long long pass_energy = 0;
    long long pass_amplitude = 0;
    long long pass_geometry = 0;
    long long pass_background = 0;
    
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        total_events++;
        
        // Use the global variables
        if (neutron_energy < emin_global || neutron_energy > emax_global) continue;
        pass_energy++;
        
        int e_bin = findBin(cfg.energy_bins, neutron_energy);
        if (e_bin < 0 || e_bin >= nbins) continue;
        if (neutron_energy > 2000) continue;
        
        EventCuts c = getCuts(cfg.sample, neutron_energy);
        
        if (!passAmplitudeCut(amp0, amp1, c)) continue;
        pass_amplitude++;
        
        if (cos_theta_det < 0.0) continue;
        if (std::fabs(cos_theta_det) > 1 || std::fabs(cos_theta) > 1) continue;
        pass_geometry++;
        
        double dt = tof1 - tof0;
        double ampsum = amp0 + amp1;
        
        // DEBUG: Print event info for first few events
        if (debug && i < 10) {
            std::cout << "Event " << i << ": E=" << neutron_energy 
                      << ", dt=" << dt << ", ampsum=" << ampsum 
                      << ", roi=[" << c.roi_min << "," << c.roi_max << "]";
            if (cfg.sample == Sample::uranium) {
                std::cout << ", bkg_min=" << c.bkg_min << ", bkg_amp_min=" << c.bkg_amp_min;
            }
            std::cout << "\n";
        }
        
        // Check if event is in background region using the NEW functions
        if (cfg.sample == Sample::uranium) {
            // Use the uranium background function from cuts.h
            if (inBackgroundUranium(dt, ampsum, c)) {
                pass_background++;
                hists_bkg[e_bin]->Fill(fabs(cos_theta));
                
                if (debug && pass_background < 5) {
                    std::cout << "  -> Uranium background event found! dt=" << dt 
                              << ", ampsum=" << ampsum << "\n";
                }
            }
        } 
        else if (cfg.sample == Sample::gold) {
            // Use the gold background function from cuts.h
            if (inBackgroundGold(dt, ampsum, c)) {
                pass_background++;
                hists_bkg[e_bin]->Fill(fabs(cos_theta));
                
                if (debug && pass_background < 5) {
                    std::cout << "  -> Gold background event found! dt=" << dt 
                              << ", ampsum=" << ampsum << "\n";
                }
            }
        }
    }
    
    // Print debug statistics
    std::cout << "\n=== Debug Statistics for " << (is_accidental ? "Accidental" : "Data") << " ===\n";
    std::cout << "Total events: " << total_events << "\n";
    std::cout << "Pass energy cut: " << pass_energy << "\n";
    std::cout << "Pass amplitude cut: " << pass_amplitude << "\n";
    std::cout << "Pass geometry cut: " << pass_geometry << "\n";
    std::cout << "Pass background: " << pass_background << "\n";
    std::cout << "================================\n\n";
}

void angular_distributions() {
    // Energy bins for analysis
    std::vector<double> energy_bins = {1,1000};
    AnalysisConfig cfg_uranium = makeUraniumConfig(energy_bins, "uranium");
    AnalysisConfig cfg_gold = makeGoldConfig(energy_bins, "gold");
    
    // Open data files
    TFile* fdata = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection_exp.root");
    TFile* facc = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/data/acc_selection.root");
    
    if (!fdata || fdata->IsZombie() || !facc || facc->IsZombie()) {
        std::cerr << "Error opening files\n";
        return;
    }
    
    // Get trees
    TTree* tree_data_uranium = (TTree*)fdata->Get("events_uranium");
    TTree* tree_data_gold = (TTree*)fdata->Get("events_gold");
    TTree* tree_acc_uranium = (TTree*)facc->Get("events_uranium");
    TTree* tree_acc_gold = (TTree*)facc->Get("events_gold");
    
    // Try alternative tree names if not found
    if (!tree_data_uranium) tree_data_uranium = (TTree*)fdata->Get("tree");
    if (!tree_data_gold) tree_data_gold = (TTree*)fdata->Get("tree");
    if (!tree_acc_uranium) tree_acc_uranium = (TTree*)facc->Get("tree");
    if (!tree_acc_gold) tree_acc_gold = (TTree*)facc->Get("tree");
    
    if (!tree_data_uranium || !tree_data_gold || !tree_acc_uranium || !tree_acc_gold) {
        std::cerr << "Some trees not found\n";
        std::cerr << "Data Uranium: " << (tree_data_uranium ? "FOUND" : "MISSING") << "\n";
        std::cerr << "Data Gold: " << (tree_data_gold ? "FOUND" : "MISSING") << "\n";
        std::cerr << "Acc Uranium: " << (tree_acc_uranium ? "FOUND" : "MISSING") << "\n";
        std::cerr << "Acc Gold: " << (tree_acc_gold ? "FOUND" : "MISSING") << "\n";
        fdata->Close();
        facc->Close();
        return;
    }
    
    // Separate vectors for data and accidental histograms
    std::vector<TH1D*> bkg_data_uranium;
    std::vector<TH1D*> bkg_acc_uranium;
    std::vector<TH1D*> bkg_data_gold;
    std::vector<TH1D*> bkg_acc_gold;
    
    // Fill histograms for uranium - data and accidental go to DIFFERENT vectors
    std::cout << "\n=== Processing Uranium ===\n";
    std::cout << "Data Uranium entries: " << tree_data_uranium->GetEntries() << "\n";
    std::cout << "Accidental Uranium entries: " << tree_acc_uranium->GetEntries() << "\n";
    fill_histograms(tree_data_uranium, cfg_uranium, bkg_data_uranium, false, true);
    fill_histograms(tree_acc_uranium, cfg_uranium, bkg_acc_uranium, true, false);
    
    // Fill histograms for gold - data and accidental go to DIFFERENT vectors
    std::cout << "\n=== Processing Gold ===\n";
    std::cout << "Data Gold entries: " << tree_data_gold->GetEntries() << "\n";
    std::cout << "Accidental Gold entries: " << tree_acc_gold->GetEntries() << "\n";
    fill_histograms(tree_data_gold, cfg_gold, bkg_data_gold, false, true);
    fill_histograms(tree_acc_gold, cfg_gold, bkg_acc_gold, true, false);
    
    // Create output file for comparison
    TFile* fout = TFile::Open("background_comparison.root", "RECREATE");
    fout->cd();
    
    // Compare uranium background
    std::cout << "\n=== Uranium Background Comparison ===\n";
    for (size_t i = 0; i < bkg_data_uranium.size(); i++) {
        TH1D* h_data = bkg_data_uranium[i];
        TH1D* h_acc = bkg_acc_uranium[i];
        
        if (!h_data || !h_acc) {
            std::cout << "Energy bin " << i << ": NULL histogram\n";
            continue;
        }
        
        double integral_data = h_data->Integral();
        double integral_acc = h_acc->Integral();
        
        std::cout << "Uranium bin " << i << ": Data=" << integral_data 
                  << ", Acc=" << integral_acc << "\n";
        
        if (integral_data > 0 && integral_acc > 0) {
            // Normalize
            h_data->Scale(1.0 / integral_data);
            h_acc->Scale(1.0 / integral_acc);
            
            h_data->SetLineColor(kBlue);
            h_data->SetLineWidth(2);
            h_data->SetStats(0);
            
            h_acc->SetLineColor(kRed);
            h_acc->SetLineWidth(2);
            h_acc->SetStats(0);
            
            h_data->Write();
            h_acc->Write();
            
            TH1D* h_ratio = (TH1D*)h_data->Clone(Form("ratio_uranium_%zu", i));
            h_ratio->SetTitle(Form("Data/Accidental Ratio - Uranium E_bin_%zu", i));
            h_ratio->SetLineColor(kGreen);
            h_ratio->SetLineWidth(2);
            h_ratio->Divide(h_acc);
            h_ratio->Write();
            
            // Calculate chi2
            double chi2 = 0;
            int ndf = 0;
            for (int bin = 1; bin <= h_data->GetNbinsX(); bin++) {
                double data_val = h_data->GetBinContent(bin);
                double acc_val = h_acc->GetBinContent(bin);
                double data_err = h_data->GetBinError(bin);
                double acc_err = h_acc->GetBinError(bin);
                
                if (data_val > 0 && acc_val > 0) {
                    double diff = data_val - acc_val;
                    double err = sqrt(data_err*data_err + acc_err*acc_err);
                    if (err > 0) {
                        chi2 += (diff*diff) / (err*err);
                        ndf++;
                    }
                }
            }
            
            if (ndf > 0) {
                std::cout << "  -> Chi2/NDF = " << chi2/ndf << " (NDF=" << ndf << ")\n";
            }
        } else {
            std::cout << "  -> Skipping bin " << i << " (zero entries)\n";
        }
    }
    
    // Compare gold background
    std::cout << "\n=== Gold Background Comparison ===\n";
    for (size_t i = 0; i < bkg_data_gold.size(); i++) {
        TH1D* h_data = bkg_data_gold[i];
        TH1D* h_acc = bkg_acc_gold[i];
        
        if (!h_data || !h_acc) {
            std::cout << "Energy bin " << i << ": NULL histogram\n";
            continue;
        }
        
        double integral_data = h_data->Integral();
        double integral_acc = h_acc->Integral();
        
        std::cout << "Gold bin " << i << ": Data=" << integral_data 
                  << ", Acc=" << integral_acc << "\n";
        
        if (integral_data > 0 && integral_acc > 0) {
            h_data->Scale(1.0 / integral_data);
            h_acc->Scale(1.0 / integral_acc);
            
            h_data->SetLineColor(kBlue);
            h_data->SetLineWidth(2);
            h_data->SetStats(0);
            
            h_acc->SetLineColor(kRed);
            h_acc->SetLineWidth(2);
            h_acc->SetStats(0);
            
            h_data->Write();
            h_acc->Write();
            
            TH1D* h_ratio = (TH1D*)h_data->Clone(Form("ratio_gold_%zu", i));
            h_ratio->SetTitle(Form("Data/Accidental Ratio - Gold E_bin_%zu", i));
            h_ratio->SetLineColor(kGreen);
            h_ratio->SetLineWidth(2);
            h_ratio->Divide(h_acc);
            h_ratio->Write();
            
            // Calculate chi2
            double chi2 = 0;
            int ndf = 0;
            for (int bin = 1; bin <= h_data->GetNbinsX(); bin++) {
                double data_val = h_data->GetBinContent(bin);
                double acc_val = h_acc->GetBinContent(bin);
                double data_err = h_data->GetBinError(bin);
                double acc_err = h_acc->GetBinError(bin);
                
                if (data_val > 0 && acc_val > 0) {
                    double diff = data_val - acc_val;
                    double err = sqrt(data_err*data_err + acc_err*acc_err);
                    if (err > 0) {
                        chi2 += (diff*diff) / (err*err);
                        ndf++;
                    }
                }
            }
            
            if (ndf > 0) {
                std::cout << "  -> Chi2/NDF = " << chi2/ndf << " (NDF=" << ndf << ")\n";
            }
        } else {
            std::cout << "  -> Skipping bin " << i << " (zero entries)\n";
        }
    }
    
    fout->Close();
    fdata->Close();
    facc->Close();
    
    std::cout << "\nComparison complete. Results saved to background_comparison.root\n";
}