#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/constants.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/types.h"

void efficiency_toy(){

    TFile *fin = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection.root", "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "Cannot open data file\n"; return; }
    TTree *tin = (TTree*)fin->Get("events_uranium");
    if (!tin) { std::cerr << "Tree not found\n"; return; }


    TFile *fcut1 = TFile::Open("/Users/nico/Desktop/Tese/Analysis/amps_sum_uranium.root", "READ");
    if (!fcut1 || fcut1->IsZombie()) { std::cerr << "Cannot open cut1 file\n"; return; }
    TCutG *cut1 = (TCutG*)fcut1->Get("amps_sum");
    if (!cut1) { std::cerr << "TCutG gold1 not found\n"; return; }

    TFile *fcut2 = TFile::Open("/Users/nico/Desktop/Tese/Analysis/amps_uranium.root", "READ");
    if (!fcut2 || fcut2->IsZombie()) { std::cerr << "Cannot open cut1 file\n"; return; }
    TCutG *cut2 = (TCutG*)fcut2->Get("amps");
    if (!cut2) { std::cerr << "TCutG gold1 not found\n"; return; }

    TFile *solid_angle = TFile::Open("/Users/nico/Desktop/Tese/Analysis/mc_acceptance.root", "READ");
    TH2D* hist_theta_det = (TH2D*) solid_angle->Get("theta_det_beam");
    TH1D *acceptance_hist = (TH1D*) solid_angle->Get("acceptance");

    std::vector<std::vector<double>> acceptance(nbins_beam, std::vector<double>(nbins_det, 0.0));
    for (int i = 1 ; i <= nbins_beam; i++){  
        for (int j = 1 ; j <= nbins_det; j++){  
            acceptance[i-1][j-1] = hist_theta_det->GetBinContent(i, j);
        }
    }


    const int nbins = 10;
    std::vector<double> energy_bins = buildLogBins(nbins, 1.0, 1000.0);
    std::vector<double> E_low(nbins), E_high(nbins);
    for (int e = 0; e < nbins; ++e) {
        E_low[e]  = energy_bins[e];
        E_high[e] = energy_bins[e+1];
    }

     // event loop 
    double tof1, tof0, neutron_energy;
    double  amp0, amp1;
    double x0,x1,y0,y1;
    double cos_theta, cos_theta_det;
    tin->SetBranchAddress("tof1",           &tof1);
    tin->SetBranchAddress("tof0",           &tof0);
    tin->SetBranchAddress("amp0",           &amp0);
    tin->SetBranchAddress("amp1",           &amp1);
    tin->SetBranchAddress("x1",           &x1);
    tin->SetBranchAddress("x0",           &x0);
    tin->SetBranchAddress("y0",           &y0);
    tin->SetBranchAddress("y1",           &y1);
    tin->SetBranchAddress("neutron_energy", &neutron_energy);
    tin->SetBranchAddress("cos_theta", &cos_theta);
    tin->SetBranchAddress("cos_theta_det", &cos_theta_det);
    Vec3D counts(nbins,
    Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Long64_t nentries = tin->GetEntries();
    for (int i = 0 ; i < nentries; i++){
        tin->GetEntry(i);
        double dt = tof1-tof0;
        if (std::fabs(cos_theta)>1 || std::fabs(cos_theta_det)>1) continue; // no cathodes available -> -9999
        if (neutron_energy < 1.0 || neutron_energy > 1000.0) continue;

        int bin = findBin(energy_bins, neutron_energy);
        if (bin < 0 || bin >= nbins) continue;
        int bin_beam = int(std::fabs(cos_theta)/dcos_beam);
        int bin_det = int(std::fabs(cos_theta_det)/dcos_det);
        if (cut1->IsInside(amp1+amp0,dt) && cut2->IsInside(amp1, amp0)){
            counts[bin][bin_beam][bin_det]++;
        }
    }
    std::cout << "Cheguei" << std::endl;
    fin->Close();
    fcut1->Close();
    fcut2->Close();

    std::vector<TH1D*> efficiencies (nbins, nullptr);
    for (int i = 0 ; i < nbins ; i++){
        TH1D* histo = new TH1D(Form("eff_ebin_%d", i), "", nbins_det, 0, 1);
        efficiencies[i] = histo;
        efficiencies[i]->SetDirectory(0);
    }

for (int ebin = 0; ebin < nbins; ebin++){
    int ref_bin = nbins_det - 1;
    
    // Set reference bin (highest cos(theta_det)) to 1
    efficiencies[ebin]->SetBinContent(ref_bin + 1, 1.0);
    efficiencies[ebin]->SetBinError(ref_bin + 1, 0.0);
    
    // Loop backwards from second-to-last bin to first bin (including bin 0)
    for (int bin_det = ref_bin - 1; bin_det >= 0; bin_det--){
        double counts_det = 0.0;
        double counts_ref = 0.0;
        int n_overlap = 0;  // Track number of valid beam bins
        
        // Accumulate counts over beam bins
        for (int bin_beam = 0; bin_beam < nbins_beam; bin_beam++){
            // Check if acceptance is valid for both detector bins
            if (acceptance[bin_beam][bin_det] <= 0.0 || acceptance[bin_beam][ref_bin] <= 0.0) 
                continue;
            
            // Check if counts are valid
            if (counts[ebin][bin_beam][ref_bin] <= 0.0 || counts[ebin][bin_beam][bin_det] <= 0.0) 
                continue;
            
            // Accumulate counts corrected by acceptance
            counts_det += counts[ebin][bin_beam][bin_det] / acceptance[bin_beam][bin_det];
            counts_ref += counts[ebin][bin_beam][ref_bin] / acceptance[bin_beam][ref_bin];
            n_overlap++;
        }
        
        // Check if we have valid counts
        if (counts_det <= 0.0 || counts_ref <= 0.0 || n_overlap == 0) {
            // Set to 0 instead of nan
            efficiencies[ebin]->SetBinContent(bin_det + 1, 0.0);
            efficiencies[ebin]->SetBinError(bin_det + 1, 0.0);
            
            // Print warning for debugging
            std::cout << "WARNING: Ebin " << ebin << ", bin_det " << bin_det 
                      << ": counts_det=" << counts_det 
                      << ", counts_ref=" << counts_ref 
                      << ", n_overlap=" << n_overlap << std::endl;
            
            ref_bin = bin_det;
            continue;
        }
        
        // Get reference efficiency and its uncertainty
        double eff_ref_bin = efficiencies[ebin]->GetBinContent(ref_bin + 1);
        double u_eff_ref_bin = efficiencies[ebin]->GetBinError(ref_bin + 1);
        
        // Check if reference efficiency is valid
        if (!std::isfinite(eff_ref_bin) || eff_ref_bin <= 0.0) {
            std::cout << "WARNING: Ebin " << ebin << ", bin_det " << bin_det 
                      << ": eff_ref_bin=" << eff_ref_bin << " is invalid!" << std::endl;
            
            efficiencies[ebin]->SetBinContent(bin_det + 1, 0.0);
            efficiencies[ebin]->SetBinError(bin_det + 1, 0.0);
            ref_bin = bin_det;
            continue;
        }
        
        // Calculate efficiency
        double eff = eff_ref_bin * (counts_det / counts_ref);
        
        // Calculate uncertainty with proper error propagation
        // Include Poisson uncertainties from counts
        double rel2 = 0.0;
        
        // Uncertainty from reference efficiency
        if (u_eff_ref_bin > 0.0 && eff_ref_bin > 0.0) {
            rel2 += pow(u_eff_ref_bin / eff_ref_bin, 2);
        }
        
        // Poisson uncertainty from counts_det (after acceptance correction)
        if (counts_det > 0.0) {
            // Need to properly propagate uncertainty through acceptance correction
            // For Poisson, variance = counts, but after dividing by acceptance:
            // var(counts_det) = sum(var(counts_raw)/acceptance^2)
            // For simplicity, we use the corrected count as an estimate
            rel2 += 1.0 / counts_det;
        }
        
        // Poisson uncertainty from counts_ref (after acceptance correction)
        if (counts_ref > 0.0) {
            rel2 += 1.0 / counts_ref;
        }
        
        double u_eff = 0.0;
        if (rel2 > 0.0 && eff > 0.0) {
            u_eff = eff * std::sqrt(rel2);
        } else {
            u_eff = 0.0;
        }
        
        // Final check for nan/inf
        if (!std::isfinite(eff) || !std::isfinite(u_eff)) {
            std::cout << "WARNING: Ebin " << ebin << ", bin_det " << bin_det 
                      << ": eff=" << eff << ", u_eff=" << u_eff << " are invalid!" << std::endl;
            
            efficiencies[ebin]->SetBinContent(bin_det + 1, 0.0);
            efficiencies[ebin]->SetBinError(bin_det + 1, 0.0);
        } else {
            // Store valid efficiency
            efficiencies[ebin]->SetBinContent(bin_det + 1, eff);
            efficiencies[ebin]->SetBinError(bin_det + 1, u_eff);
        }
        
        // Update reference bin for next iteration
        ref_bin = bin_det;
    }
}

// Print summary of efficiencies to check for nan/inf
std::cout << "\n=== EFFICIENCY SUMMARY ===\n";
for (int ebin = 0; ebin < nbins; ebin++) {
    std::cout << "\nEbin " << ebin << ":\n";
    int valid_bins = 0;
    for (int bin = 1; bin <= nbins_det; bin++) {
        double content = efficiencies[ebin]->GetBinContent(bin);
        double error = efficiencies[ebin]->GetBinError(bin);
        
        if (std::isfinite(content) && std::isfinite(error)) {
            if (content > 0.0 || error > 0.0) {
                valid_bins++;
                if (valid_bins <= 5) {  // Print first 5 valid bins
                    printf("  bin %d: eff=%.6f ± %.6f\n", bin, content, error);
                }
            }
        } else {
            std::cout << "  bin " << bin << ": INVALID (nan/inf)\n";
        }
    }
    std::cout << "  Total valid bins: " << valid_bins << " / " << nbins_det << std::endl;
}
    /*
    for (int ebin = 0; ebin < nbins; ebin++){
        for (int bin_det = 0; bin_det < nbins_det; bin_det++){
            double counts_det = 0.0;
            for (int bin_beam = 0; bin_beam < nbins_beam; bin_beam++){
                if (acceptance[bin_beam][bin_det] == 0) continue;
                counts_det += counts[ebin][bin_beam][bin_det];
            }
            efficiencies[ebin]->SetBinContent(bin_det+1, counts_det);
            efficiencies[ebin]->SetBinError(bin_det+1, std::sqrt(counts_det));
        }
    }
    */
    
    
    std::vector<TH1D*> efficiencies_beam (nbins, nullptr);
    for (int i = 0 ; i < nbins ; i++){
        // El eje X ahora es nbins_beam (ángulo del haz), no nbins_det
        TH1D* histo = new TH1D(Form("eff_ebin_beam_%d", i), "", nbins_beam, 0, 1);
        efficiencies_beam[i] = histo;
        efficiencies_beam[i]->SetDirectory(0);
    }

    for (int ebin = 0; ebin < nbins; ebin++){
        for (int bin_beam = 0; bin_beam < nbins_beam; bin_beam++){
            double counts_beam = 0.0;
            for (int bin_det = 0; bin_det < nbins_det; bin_det++){ 
                if (acceptance[bin_beam][bin_det] == 0) continue;
                counts_beam += counts[ebin][bin_beam][bin_det];
            }
            efficiencies_beam[ebin]->SetBinContent(bin_beam+1, counts_beam);
            efficiencies_beam[ebin]->SetBinError(bin_beam+1, std::sqrt(counts_beam));
        }
    }
auto counts_original = counts;

std::vector<TH1D*> h_eff_check(nbins, nullptr);
std::vector<TH1D*> h_avg_acceptance(nbins, nullptr);

for (int ebin = 0; ebin < nbins; ebin++) {
    h_eff_check[ebin] = new TH1D(Form("h_eff_check_ebin%d", ebin), 
                                  Form("Efficiency Check - Ebin %d;cos#theta_{det};N_{det}/N_{est}", ebin),
                                  nbins_det, 0, 1);
    h_avg_acceptance[ebin] = new TH1D(Form("h_avg_acceptance_ebin%d", ebin), 
                                       Form("Average Acceptance - Ebin %d;cos#theta_{det};<Acceptance>", ebin),
                                       nbins_det, 0, 1);
    
    for (int bin_det = 0; bin_det < nbins_det; bin_det++) {
        double N_detected = 0;
        double N_estimated = 0;
        double sum_acc = 0;
        int n_bins = 0;
        
        for (int bin_beam = 0; bin_beam < nbins_beam; bin_beam++) {
            double cnt = counts_original[ebin][bin_beam][bin_det];
            double acc = acceptance[bin_beam][bin_det];
            
            N_detected += cnt;
            if (acc > 0 && cnt > 0) {
                N_estimated += cnt / acc;
                sum_acc += acc;
                n_bins++;
            }
        }
        
        double eff_check = (N_estimated > 0) ? N_detected / N_estimated : 0;
        double avg_acc = (n_bins > 0) ? sum_acc / n_bins : 0;
        
        h_eff_check[ebin]->SetBinContent(bin_det + 1, eff_check);
        h_avg_acceptance[ebin]->SetBinContent(bin_det + 1, avg_acc);
    }
}
std::vector<TH1D*> h_eff_check_beam(nbins, nullptr);
std::vector<TH1D*> h_avg_acceptance_beam(nbins, nullptr);

for (int ebin = 0; ebin < nbins; ebin++) {
    h_eff_check_beam[ebin] = new TH1D(Form("h_eff_check_beam_ebin%d", ebin), 
                                       Form("Efficiency Check - Beam Ebin %d;cos#theta_{beam};N_{det}/N_{est}", ebin),
                                       nbins_beam, 0, 1);
    h_avg_acceptance_beam[ebin] = new TH1D(Form("h_avg_acceptance_beam_ebin%d", ebin), 
                                            Form("Average Acceptance - Beam Ebin %d;cos#theta_{beam};<Acceptance>", ebin),
                                            nbins_beam, 0, 1);
    
    for (int bin_beam = 0; bin_beam < nbins_beam; bin_beam++) {
        double N_detected = 0;
        double N_estimated = 0;
        double sum_acc = 0;
        int n_bins = 0;
        
        for (int bin_det = 0; bin_det < nbins_det; bin_det++) {
            double cnt = counts_original[ebin][bin_beam][bin_det];
            double acc = acceptance[bin_beam][bin_det];
            
            N_detected += cnt;
            if (acc > 0 && cnt > 0) {
                N_estimated += cnt / acc;
                sum_acc += acc;
                n_bins++;
            }
        }
        
        double eff_check = (N_estimated > 0) ? N_detected / N_estimated : 0;
        double avg_acc = (n_bins > 0) ? sum_acc / n_bins : 0;
        
        h_eff_check_beam[ebin]->SetBinContent(bin_beam + 1, eff_check);
        h_avg_acceptance_beam[ebin]->SetBinContent(bin_beam + 1, avg_acc);
    }
}

std::cout << "\n=== ACCEPTANCE CORRECTION CHECK ===\n";
for (int ebin = 0; ebin < nbins; ebin++) {
    std::cout << "\n--- Ebin " << ebin << " (E = " << E_low[ebin] << "-" << E_high[ebin] << " MeV) ---\n";
    std::cout << "bin_det  cos_theta_det    N_det  N_est  eff_check  avg_acc\n";
    std::cout << "-----------------------------------------------------------\n";
    
    for (int bin_det = 0; bin_det < nbins_det; bin_det++) {
        double N_det = 0, N_est = 0;
        
        for (int bin_beam = 0; bin_beam < nbins_beam; bin_beam++) {
            double cnt = counts_original[ebin][bin_beam][bin_det];
            double acc = acceptance[bin_beam][bin_det];
            N_det += cnt;
            if (acc > 0) N_est += cnt / acc;
        }
        
        double cos_det_lo = bin_det * dcos_det;
        
        if (N_det > 0) {
            printf("%d\t[%.3f,%.3f)\t%.0f\t%.0f\t%.4f\t%.4f\n",
                   bin_det, cos_det_lo, cos_det_lo + dcos_det,
                   N_det, N_est,
                   h_eff_check[ebin]->GetBinContent(bin_det+1),
                   h_avg_acceptance[ebin]->GetBinContent(bin_det+1));
        }
    }
}

std::cout << "Cheguei 3" << std::endl;

TFile *fout = new TFile("/Users/nico/Desktop/Tese/Analysis/cross_section/efficiencies_u_toy.root", "RECREATE");

for (int i = 0; i < nbins; i++) {
    efficiencies[i]->Write();
    efficiencies_beam[i]->Write();
}

// Guardar histogramas de chequeo
for (int i = 0; i < nbins; i++) {
    h_eff_check[i]->Write();
    h_avg_acceptance[i]->Write();
    h_eff_check_beam[i]->Write();
}

fout->Close();
solid_angle->Close();
}